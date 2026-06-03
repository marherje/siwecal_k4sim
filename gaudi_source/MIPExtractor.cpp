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
// in finalize() to find the MIP peak. Results are printed and saved to files.
//
// Based on the CALICE ConversionProcessor (Marlin) by E. Brianne (DESY, 2015),
// adapted for the Gaudi/EDM4hep framework.
//
// Fit modes:
//   1 = Gaussian (fit FWHM region around histogram maximum)
//   2 = Landau   (fit in [mean-RMS, mean+2*RMS], MPV corrected by +0.222*sigma)
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
      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());

      m_layerHistos.clear();
      m_layerHistos.reserve(nLayers);
      for (int i = 0; i < nLayers; i++) {
        const std::string hname = "MIP_layer_" + std::to_string(i);
        m_layerHistos.push_back(
            std::make_unique<TH1F>(hname.c_str(), hname.c_str(), 200, 0., 0.001));
      }

      info() << "[MIPExtractor] Initialized " << nLayers
             << " per-layer histograms (fit mode=" << m_fitMode.value() << ")." << endmsg;
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
      const int fitMode = m_fitMode.value();
      const int minEntries = m_minEntries.value();
      std::vector<double> mipValues;

      // Open ROOT output
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
          mip = (fitMode == 2) ? fitLandau(h) : fitGaus(h);
        } else {
          warning() << "[MIPExtractor] Layer " << iz << " has only "
                    << h->GetEntries() << " entries (< MinEntries=" << minEntries
                    << ") — skipping fit, MIP set to 0." << endmsg;
        }
        mipValues.push_back(mip);
        if (rootOut) { rootOut->cd(); h->Write(); }
      }

      if (rootOut) {
        rootOut->Write("", TObject::kOverwrite);
        rootOut->Close();
        delete rootOut;
      }

      // Print results
      info() << "[MIPExtractor] === Per-layer MIP values (GeV) ===" << endmsg;
      for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
        info() << "  Layer " << i << ": " << mipValues[i] << " GeV" << endmsg;
      }
      info() << "[MIPExtractor] MIPValues vector for job3_digitize.py:" << endmsg;
      std::string vec = "mip.MIPValues = [";
      for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
        vec += std::to_string(mipValues[i]);
        if (i + 1 < static_cast<int>(mipValues.size())) vec += ", ";
      }
      vec += "]";
      info() << "  " << vec << endmsg;

      // Write text file
      if (!m_outputTextFile.value().empty()) {
        std::ofstream ofs(m_outputTextFile.value());
        ofs << "# MIPExtractor per-layer MIP values [GeV]\n";
        ofs << "# FitMode=" << fitMode << "  NLayers=" << mipValues.size() << "\n";
        ofs << "# Layer  MIP_GeV\n";
        for (int i = 0; i < static_cast<int>(mipValues.size()); i++) {
          ofs << i << "  " << mipValues[i] << "\n";
        }
        ofs << "\n# Paste into job3_digitize.py:\n# " << vec << "\n";
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
  // Gaussian fit: fit around the FWHM region of the histogram maximum.
  // Returns the Gaussian mean (= MIP peak position).
  double fitGaus(TH1F* h) const {
    const int    maxBin = h->GetMaximumBin();
    const double halfMax = h->GetBinContent(maxBin) / 2.0;
    double minFit = h->GetXaxis()->GetXmin();
    double maxFit = h->GetXaxis()->GetXmax();

    for (int b = 1; b < maxBin; b++) {
      if (h->GetBinContent(b) > halfMax) { minFit = h->GetBinCenter(b); break; }
    }
    for (int b = maxBin + 1; b <= h->GetNbinsX(); b++) {
      if (h->GetBinContent(b) < halfMax) { maxFit = h->GetBinCenter(b); break; }
    }

    TF1 f("fGaus", "gaus", minFit, maxFit);
    f.SetParameters(h->GetMaximum(), h->GetBinCenter(maxBin), (maxFit - minFit) / 2.0);
    h->Fit(&f, "RQ0");
    return f.GetParameter(1);
  }

  // Landau fit: fit in [mean-RMS, mean+2*RMS].
  // Returns MPV corrected by +0.22278*sigma (Landau shift correction).
  double fitLandau(TH1F* h) const {
    const double minFit = TMath::Max(h->GetMean() - h->GetRMS(), 0.0);
    const double maxFit = h->GetMean() + 2.0 * h->GetRMS();

    TF1 f("fLandau", "[0]*TMath::Landau(x,[1],[2],1)", minFit, maxFit);
    f.SetParameter(0, h->GetMaximum());
    f.SetParameter(1, h->GetMean());
    f.SetParameter(2, h->GetRMS() / 10.0);
    f.SetParLimits(1, 0.0, h->GetMean() * 2.0);
    f.SetParLimits(2, h->GetRMS() / 100.0, h->GetRMS());
    h->Fit(&f, "RQ0");
    return f.GetParameter(1) + 0.22278298 * f.GetParameter(2);
  }

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
      this, "FitMode", 1,
      "Fit mode for MIP peak extraction: 1=Gaussian, 2=Landau"};
  Gaudi::Property<int> m_minEntries{
      this, "MinEntries", 50,
      "Minimum histogram entries required to attempt a fit"};
  Gaudi::Property<std::string> m_outputRootFile{
      this, "OutputRootFile", "mip_extraction.root",
      "ROOT output file for per-layer energy histograms (empty = skip)"};
  Gaudi::Property<std::string> m_outputTextFile{
      this, "OutputTextFile", "mip_values.txt",
      "Text file with per-layer MIP values ready to paste into job3 (empty = skip)"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::vector<std::unique_ptr<TH1F>> m_layerHistos;
};

DECLARE_COMPONENT(MIPExtractor)
