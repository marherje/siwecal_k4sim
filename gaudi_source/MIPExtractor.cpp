#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>

#include <fstream>
#include <memory>
#include <vector>

// Extracts per-layer MIP peak values from raw SimCalorimeterHit data (in GeV).
// Accumulates per-layer energy histograms across all events, then fits them
// in finalize() to find the MIP peak.
//
// Output ROOT file contains:
//   - MIP_layer_N   TH1F: per-layer energy histogram with the fit function attached
//   - The fit TF1 is embedded in each histogram's function list (visible in TBrowser)
//
// Based on the CALICE ConversionProcessor (Marlin) by E. Brianne (DESY, 2015),
// adapted for the Gaudi/EDM4hep framework.
//
// Fit modes:
//   1 = Gaussian  — fit FWHM region around histogram maximum; returns mean
//   2 = Landau    — fit in [mean-RMS, mean+2*RMS]; returns MPV+0.22278*sigma
//   3 = LanGaus   — numerical Landau⊗Gaussian convolution (best for MIP peaks);
//                   returns MPV+0.22278*sigma (same convention as Landau)
class MIPExtractor : public Gaudi::Algorithm {
public:
  MIPExtractor(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      m_inputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(), Gaudi::DataHandle::Reader, this);

      const int nLayers = m_nLayers.value();
      if (nLayers <= 0) {
        error() << "[MIPExtractor] NLayers must be > 0." << endmsg;
        return StatusCode::FAILURE;
      }
      const int fitMode = m_fitMode.value();
      if (fitMode < 1 || fitMode > 3) {
        error() << "[MIPExtractor] FitMode must be 1 (Gauss), 2 (Landau) or 3 (LanGaus)." << endmsg;
        return StatusCode::FAILURE;
      }
      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());

      m_layerHistos.clear();
      m_layerHistos.reserve(nLayers);
      const double histMax = m_histMax.value();
      const int    histBins = m_histBins.value();
      for (int i = 0; i < nLayers; i++) {
        const std::string hname = "MIP_layer_" + std::to_string(i);
        const std::string htitle = "Energy per hit, layer " + std::to_string(i) + ";Energy [GeV];Hits";
        m_layerHistos.push_back(
            std::make_unique<TH1F>(hname.c_str(), htitle.c_str(), histBins, 0., histMax));
      }

      const char* modeStr[] = {"", "Gaussian", "Landau", "Landau#otimesGaussian"};
      info() << "[MIPExtractor] Initialized " << nLayers
             << " per-layer histograms. FitMode=" << fitMode
             << " (" << modeStr[fitMode] << ")." << endmsg;
      return sc;
    } catch (const std::exception& e) {
      error() << "[MIPExtractor] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[MIPExtractor] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input = m_inputHandle->get();
      for (const auto& hit : *input) {
        const int layer = static_cast<int>(m_decoder->get(hit.getCellID(), "layer"));
        if (layer < 0 || layer >= static_cast<int>(m_layerHistos.size())) continue;
        m_layerHistos[layer]->Fill(hit.getEnergy());
      }
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[MIPExtractor] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[MIPExtractor] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      const int fitMode    = m_fitMode.value();
      const int minEntries = m_minEntries.value();
      std::vector<double> mipValues;

      // Open ROOT output file.
      // Each histogram is written with its fit function embedded in GetListOfFunctions()
      // so that both histogram and fit are visible when opening in TBrowser/ROOT.
      TFile* rootOut = nullptr;
      if (!m_outputRootFile.value().empty()) {
        rootOut = TFile::Open(m_outputRootFile.value().c_str(), "RECREATE");
        if (!rootOut || rootOut->IsZombie()) {
          warning() << "[MIPExtractor] Could not open ROOT output file '"
                    << m_outputRootFile.value() << "'. Histograms will not be saved." << endmsg;
          rootOut = nullptr;
        }
      }

      for (int iz = 0; iz < static_cast<int>(m_layerHistos.size()); iz++) {
        TH1F* h = m_layerHistos[iz].get();
        double mip = 0.0;

        if (static_cast<int>(h->GetEntries()) >= minEntries) {
          if      (fitMode == 1) mip = fitGaus(h, iz);
          else if (fitMode == 2) mip = fitLandau(h, iz);
          else                   mip = fitLanGaus(h, iz);
        } else {
          warning() << "[MIPExtractor] Layer " << iz << " has only "
                    << h->GetEntries() << " entries (< MinEntries=" << minEntries
                    << ") — skipping fit, MIP set to 0." << endmsg;
        }
        mipValues.push_back(mip);

        // Write histogram to ROOT file.
        // The fit TF1 is already attached to h via Fit(...,"RQ") — it is saved automatically.
        if (rootOut) { rootOut->cd(); h->Write(); }
      }

      if (rootOut) {
        rootOut->Close();
        delete rootOut;
        info() << "[MIPExtractor] Histograms (with fits) saved to '"
               << m_outputRootFile.value() << "'." << endmsg;
      }

      // Print per-layer results
      info() << "[MIPExtractor] === Per-layer MIP peak values (GeV) ===" << endmsg;
      for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
        info() << "  Layer " << i << ": " << mipValues[i] << " GeV" << endmsg;
      }

      // Build the Python line ready to paste into job3_digitize.py
      std::string pyLine = "mip.MIPValues = [";
      for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
        pyLine += std::to_string(mipValues[i]);
        if (i + 1 < static_cast<int>(mipValues.size())) pyLine += ", ";
      }
      pyLine += "]";
      info() << "[MIPExtractor] Ready to paste into job3_digitize.py:" << endmsg;
      info() << "  " << pyLine << endmsg;

      // Write text file
      if (!m_outputTextFile.value().empty()) {
        std::ofstream ofs(m_outputTextFile.value());
        ofs << "# MIPExtractor per-layer MIP peak values [GeV]\n";
        ofs << "# FitMode=" << fitMode << "  NLayers=" << mipValues.size() << "\n";
        ofs << "# Layer  MIP_GeV\n";
        for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
          ofs << i << "  " << mipValues[i] << "\n";
        }
        ofs << "\n# Paste into job3_digitize.py:\n" << pyLine << "\n";
        info() << "[MIPExtractor] Results written to '" << m_outputTextFile.value() << "'." << endmsg;
      }

      m_inputHandle.reset();
      m_decoder.reset();
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[MIPExtractor] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[MIPExtractor] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // ── Fit mode 1: Gaussian ────────────────────────────────────────────────────
  // Fits the FWHM region around the histogram maximum.
  // Returns the Gaussian mean = MIP peak position.
  // The fitted TF1 is attached to the histogram (saved in ROOT file).
  double fitGaus(TH1F* h, int iz) const {
    const int    maxBin  = h->GetMaximumBin();
    const double halfMax = h->GetBinContent(maxBin) / 2.0;
    double minFit = h->GetXaxis()->GetXmin();
    double maxFit = h->GetXaxis()->GetXmax();
    for (int b = 1;           b < maxBin;          b++)
      if (h->GetBinContent(b) > halfMax) { minFit = h->GetBinCenter(b); break; }
    for (int b = maxBin + 1;  b <= h->GetNbinsX(); b++)
      if (h->GetBinContent(b) < halfMax) { maxFit = h->GetBinCenter(b); break; }

    const std::string fname = "fGaus_layer_" + std::to_string(iz);
    TF1 f(fname.c_str(), "gaus", minFit, maxFit);
    f.SetParameters(h->GetMaximum(), h->GetBinCenter(maxBin), (maxFit - minFit) / 2.0);
    // "RQ": R=use TF1 range, Q=quiet. No "0" flag → fit is attached to histogram and saved.
    h->Fit(&f, "RQ");
    return f.GetParameter(1);
  }

  // ── Fit mode 2: Landau ──────────────────────────────────────────────────────
  // Fits in [mean-RMS, mean+2*RMS].
  // Returns MPV + 0.22278*sigma (corrects for the Landau distribution shift).
  double fitLandau(TH1F* h, int iz) const {
    const double minFit = TMath::Max(h->GetMean() - h->GetRMS(), 0.0);
    const double maxFit = h->GetMean() + 2.0 * h->GetRMS();

    const std::string fname = "fLandau_layer_" + std::to_string(iz);
    TF1 f(fname.c_str(), "[0]*TMath::Landau(x,[1],[2],1)", minFit, maxFit);
    f.SetParameter(0, h->GetMaximum());
    f.SetParameter(1, h->GetMean());
    f.SetParameter(2, h->GetRMS() / 10.0);
    f.SetParLimits(1, 0.0,               h->GetMean() * 2.0);
    f.SetParLimits(2, h->GetRMS()/100.0, h->GetRMS());
    h->Fit(&f, "RQ");
    return f.GetParameter(1) + 0.22278298 * f.GetParameter(2);
  }

  // ── Fit mode 3: Landau ⊗ Gaussian convolution ──────────────────────────────
  // Numerical convolution following the CALICE/SiWECAL-TB-analysis implementation.
  // par[0] = Landau width (scale),  par[1] = Landau MPV,
  // par[2] = total area,            par[3] = Gaussian sigma.
  // Returns MPV + 0.22278*par[0] (same shift correction as Landau-only).
  //
  // Reference: E. Brianne, ConversionProcessor; originally from
  //   https://accserv.lepp.cornell.edu/svn/packages/root/tutorials/fit/langaus.C
  static double langaufun(double* x, double* par) {
    constexpr double invsq2pi = 0.3989422804014;  // (2π)^(-1/2)
    constexpr double mpshift  = -0.22278298;       // Landau maximum location
    constexpr double np       = 100.0;             // convolution steps
    constexpr double sc       = 5.0;               // ± sc·sigma Gaussian range

    const double mpc  = par[1] - mpshift * par[0];
    const double xlow = x[0] - sc * par[3];
    const double xupp = x[0] + sc * par[3];
    const double step = (xupp - xlow) / np;

    double sum = 0.0;
    for (double i = 1.0; i <= np / 2; i++) {
      double xx = xlow + (i - 0.5) * step;
      sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]);
      xx = xupp - (i - 0.5) * step;
      sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]);
    }
    return par[2] * step * sum * invsq2pi / par[3];
  }

  double fitLanGaus(TH1F* h, int iz) const {
    // Fit range and start values following the original ConversionProcessor.
    const double fr0 = h->GetMean() - 0.8 * h->GetRMS();
    const double fr1 = 1.5 * h->GetMean();

    double sv[4]   = { h->GetRMS() * 0.5,
                       h->GetMean() * 0.6,
                       h->Integral("width"),
                       h->GetRMS() / 5.0 };
    double pllo[4] = { h->GetRMS() / 20.0, 0.0,              0.005, 0.0           };
    double plhi[4] = { h->GetRMS(),         h->GetMean()*2.0, 3.0,   h->GetRMS()  };

    const std::string fname = "fLanGaus_layer_" + std::to_string(iz);
    TF1 f(fname.c_str(), langaufun, fr0, fr1, 4);
    f.SetParNames("LandauWidth", "MPV", "Area", "GausSigma");
    f.SetParameters(sv);
    for (int i = 0; i < 4; i++) f.SetParLimits(i, pllo[i], plhi[i]);

    h->Fit(&f, "RQ");
    // MPV correction: actual peak = par[1] + 0.22278 * par[0]
    return f.GetParameter(1) + 0.22278298 * f.GetParameter(0);
  }

  // ── Properties ──────────────────────────────────────────────────────────────
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiPadHitsWindowed",
      "Input SimCalorimeterHit collection (raw, in GeV)"};
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:4,x:9,y:9",
      "DD4hep BitField encoding string for CellID decoding"};
  Gaudi::Property<int> m_nLayers{
      this, "NLayers", 15,
      "Number of detector layers (one histogram per layer)"};
  Gaudi::Property<int> m_fitMode{
      this, "FitMode", 3,
      "MIP peak fit: 1=Gaussian, 2=Landau, 3=Landau⊗Gaussian (best)"};
  Gaudi::Property<int> m_minEntries{
      this, "MinEntries", 50,
      "Minimum histogram entries required to attempt a fit"};
  Gaudi::Property<int> m_histBins{
      this, "HistBins", 200,
      "Number of bins in each per-layer energy histogram"};
  Gaudi::Property<double> m_histMax{
      this, "HistMax", 0.001,
      "Upper edge of per-layer energy histograms [GeV]"};
  Gaudi::Property<std::string> m_outputRootFile{
      this, "OutputRootFile", "mip_extraction.root",
      "ROOT output file: contains per-layer TH1F histograms with fits embedded "
      "(set empty to skip)"};
  Gaudi::Property<std::string> m_outputTextFile{
      this, "OutputTextFile", "mip_values.txt",
      "Text file with per-layer MIP values and a ready-to-paste Python line "
      "(set empty to skip)"};

  // ── Internal state ───────────────────────────────────────────────────────────
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::vector<std::unique_ptr<TH1F>> m_layerHistos;
};

DECLARE_COMPONENT(MIPExtractor)
