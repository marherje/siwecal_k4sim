#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "TFile.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

// Writes the digitised simulation as the test beam's "ecal" tree, in the same
// k4run as the digitisation, so the simulation reaches the reconstruction
// (EcalToEDM4hep + EcalPidTransformer) through the very tree format the data
// do, and the tree stays available for the comparisons.  Replaces the
// stand-alone analysis/sim_to_ecal_tree.py and writes the same branches, filled
// the same way (that script is kept as the reference it was validated against).
//
// Input: a hit collection with the test-beam CellID that ChannelMapper writes
// (SiPadHitsMapped on the simple chain, SiPadHitsRealAdc on the AdcDigitizer
// chain), plus the parallel UserDataCollections the chain writes next to it.
// Index alignment is the only thing tying those to the hits, so every one of
// them is optional (empty name = off, branch left at 0) and one whose length
// differs from the hits is ignored with a warning instead of being paired with
// the wrong hits.  With a KeptCollection (AdcDigitizer's chip-level hit
// selection) only the hits flagged 1 are written, in order.
class EcalTreeWriter : public Gaudi::Algorithm {
public:
  EcalTreeWriter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      m_hitsHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(), Gaudi::DataHandle::Reader, this);
      m_maskHandle = intHandle(m_maskName.value());
      m_keptHandle = intHandle(m_keptName.value());
      m_timeHandle = floatHandle(m_timeName.value());
      m_fastHandle = floatHandle(m_fastName.value());
      m_adcHighHandle = floatHandle(m_adcHighName.value());
      m_adcLowHandle = floatHandle(m_adcLowName.value());

      m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_bitField.value());
      for (const char* field : {"slab", "chip", "channel", "sca"}) {
        try {
          m_decoder->index(field);
        } catch (const std::exception&) {
          error() << "[EcalTreeWriter] BitField '" << m_bitField.value() << "' has no field '" << field
                  << "': the input must carry the test-beam CellID (ChannelMapper output)" << endmsg;
          return StatusCode::FAILURE;
        }
      }

      // Cumulative X0 in front of each layer, and the layer's own sampling
      // weight (its absorber / X0_W): the same tables as EcalPidTransformer's
      // WThicknesses and event_viewer/_metrics.py:hit_weights().
      const auto& w = m_wThicknesses.value();
      double cumulative = 0.;
      for (double t : w) {
        cumulative += t;
        m_layerX0.push_back(cumulative / m_x0W.value());
        m_layerWeight.push_back(t / m_x0W.value());
      }

      m_file.reset(TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
      if (!m_file || m_file->IsZombie()) {
        error() << "[EcalTreeWriter] cannot open output file " << m_outputFile.value() << endmsg;
        return StatusCode::FAILURE;
      }
      m_file->cd();
      m_tree = new TTree(m_treeName.value().c_str(), "ecal tree from k4sim simulation");  // owned by m_file
      bookBranches();

      info() << "[EcalTreeWriter] " << m_inputName.value() << " -> " << m_outputFile.value() << ":"
             << m_treeName.value() << "  (mask '" << m_maskName.value() << "', kept '" << m_keptName.value()
             << "', time '" << m_timeName.value() << "', fast '" << m_fastName.value() << "', ADC '"
             << m_adcHighName.value() << "'/'" << m_adcLowName.value() << "')" << endmsg;
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[EcalTreeWriter] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* hits = m_hitsHandle->get();
      const std::size_t nIn = hits->size();

      const auto* mask = parallel(m_maskHandle, nIn, m_maskWarned);
      const auto* kept = parallel(m_keptHandle, nIn, m_keptWarned);
      const auto* time = parallel(m_timeHandle, nIn, m_timeWarned);
      const auto* fast = parallel(m_fastHandle, nIn, m_fastWarned);
      const auto* adcHigh = parallel(m_adcHighHandle, nIn, m_adcWarned);
      const auto* adcLow = parallel(m_adcLowHandle, nIn, m_adcWarned);

      std::vector<std::size_t> indices;
      indices.reserve(nIn);
      for (std::size_t i = 0; i < nIn; ++i) {
        if (!kept || (*kept)[i] > 0) indices.push_back(i);
      }
      if (indices.size() > kMaxHits) {
        warning() << "[EcalTreeWriter] event " << m_event << ": " << indices.size() << " hits > " << kMaxHits
                  << ", truncating" << endmsg;
        indices.resize(kMaxHits);
      }

      std::lock_guard<std::mutex> lock(m_mutex);
      B& b = m_b;
      b.run = m_runNumber.value();
      b.event = m_event++;
      b.spill = 0;
      b.bcid = 0;
      b.nhit_chan = static_cast<int>(indices.size());

      // nhit_chip counts distinct chip NUMBERS, not (slab, chip) pairs -- as
      // sim_to_ecal_tree.py always did; kept so the tree is unchanged.
      std::set<int> slabs, chips;
      double sumEnergy = 0.;
      double sumHg = 0.;
      for (std::size_t j = 0; j < indices.size(); ++j) {
        const std::size_t i = indices[j];
        const auto hit = (*hits)[i];
        const auto cid = hit.getCellID();
        const int slab = static_cast<int>(m_decoder->get(cid, "slab"));
        const int chip = static_cast<int>(m_decoder->get(cid, "chip"));
        const float energy = hit.getEnergy();
        const auto pos = hit.getPosition();
        const bool inGeometry = slab >= 0 && static_cast<std::size_t>(slab) < m_layerX0.size();

        b.hit_slab[j] = slab;
        b.hit_chip[j] = chip;
        b.hit_chan[j] = static_cast<int>(m_decoder->get(cid, "channel"));
        b.hit_sca[j] = static_cast<int>(m_decoder->get(cid, "sca"));
        b.hit_ismasked[j] = mask ? (*mask)[i] : 0;
        b.hit_energy[j] = energy;
        b.hit_hg[j] = adcHigh ? (*adcHigh)[i] : 0.f;
        b.hit_lg[j] = adcLow ? (*adcLow)[i] : 0.f;
        b.hit_x[j] = pos.x;
        b.hit_y[j] = pos.y;
        b.hit_z[j] = pos.z;
        b.hit_X0[j] = inGeometry ? static_cast<float>(m_layerX0[slab]) : 0.f;
        b.hit_w_energy[j] = inGeometry ? static_cast<float>(energy * m_layerWeight[slab]) : 0.f;
        b.hit_time[j] = time ? (*time)[i] : 0.f;
        b.hit_fast[j] = fast ? (*fast)[i] : 0.f;

        slabs.insert(slab);
        chips.insert(chip);
        sumEnergy += energy;
        sumHg += b.hit_hg[j];
      }
      b.nhit_slab = static_cast<int>(slabs.size());
      b.nhit_chip = static_cast<int>(chips.size());
      b.sum_energy = static_cast<float>(sumEnergy);
      b.sum_hg = static_cast<float>(sumHg);
      m_tree->Fill();
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EcalTreeWriter] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    if (m_file) {
      m_file->cd();
      m_tree->Write();
      info() << "[EcalTreeWriter] wrote " << m_tree->GetEntries() << " events to " << m_outputFile.value()
             << endmsg;
      m_file->Close();
      m_file.reset();
    }
    m_hitsHandle.reset();
    m_maskHandle.reset();
    m_keptHandle.reset();
    m_timeHandle.reset();
    m_fastHandle.reset();
    m_adcHighHandle.reset();
    m_adcLowHandle.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
  static constexpr std::size_t kMaxHits = 4096;

  template <typename T>
  using UDHandle = std::unique_ptr<k4FWCore::DataHandle<podio::UserDataCollection<T>>>;

  UDHandle<int> intHandle(const std::string& name) {
    if (name.empty()) return nullptr;
    return std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<int>>>(name, Gaudi::DataHandle::Reader,
                                                                                  this);
  }
  UDHandle<float> floatHandle(const std::string& name) {
    if (name.empty()) return nullptr;
    return std::make_unique<k4FWCore::DataHandle<podio::UserDataCollection<float>>>(
        name, Gaudi::DataHandle::Reader, this);
  }

  /// The parallel collection, or nullptr when it is off or not aligned with the hits.
  template <typename T>
  const podio::UserDataCollection<T>* parallel(const UDHandle<T>& handle, std::size_t nHits, bool& warned) const {
    if (!handle) return nullptr;
    const auto* col = handle->get();
    if (col->size() != nHits) {
      if (!warned) {
        warning() << "[EcalTreeWriter] '" << handle->objKey() << "' has " << col->size() << " entries but '"
                  << m_inputName.value() << "' has " << nHits << "; its branch is set to 0" << endmsg;
        warned = true;
      }
      return nullptr;
    }
    return col;
  }

  void bookBranches() {
    B& b = m_b;
    m_tree->Branch("run", &b.run, "run/I");
    m_tree->Branch("event", &b.event, "event/I");
    m_tree->Branch("spill", &b.spill, "spill/I");
    m_tree->Branch("bcid", &b.bcid, "bcid/I");
    m_tree->Branch("nhit_chan", &b.nhit_chan, "nhit_chan/I");
    m_tree->Branch("nhit_slab", &b.nhit_slab, "nhit_slab/I");
    m_tree->Branch("nhit_chip", &b.nhit_chip, "nhit_chip/I");
    m_tree->Branch("sum_energy", &b.sum_energy, "sum_energy/F");
    m_tree->Branch("sum_hg", &b.sum_hg, "sum_hg/F");
    m_tree->Branch("hit_slab", b.hit_slab.data(), "hit_slab[nhit_chan]/I");
    m_tree->Branch("hit_chip", b.hit_chip.data(), "hit_chip[nhit_chan]/I");
    m_tree->Branch("hit_chan", b.hit_chan.data(), "hit_chan[nhit_chan]/I");
    m_tree->Branch("hit_sca", b.hit_sca.data(), "hit_sca[nhit_chan]/I");
    m_tree->Branch("hit_ismasked", b.hit_ismasked.data(), "hit_ismasked[nhit_chan]/I");
    m_tree->Branch("hit_energy", b.hit_energy.data(), "hit_energy[nhit_chan]/F");
    m_tree->Branch("hit_hg", b.hit_hg.data(), "hit_hg[nhit_chan]/F");
    m_tree->Branch("hit_lg", b.hit_lg.data(), "hit_lg[nhit_chan]/F");
    m_tree->Branch("hit_x", b.hit_x.data(), "hit_x[nhit_chan]/F");
    m_tree->Branch("hit_y", b.hit_y.data(), "hit_y[nhit_chan]/F");
    m_tree->Branch("hit_z", b.hit_z.data(), "hit_z[nhit_chan]/F");
    m_tree->Branch("hit_X0", b.hit_X0.data(), "hit_X0[nhit_chan]/F");
    m_tree->Branch("hit_w_energy", b.hit_w_energy.data(), "hit_w_energy[nhit_chan]/F");
    // Fast-channel trigger time [ns] and peak amplitude [MIP] (RealDigitizer chain only).
    m_tree->Branch("hit_time", b.hit_time.data(), "hit_time[nhit_chan]/F");
    m_tree->Branch("hit_fast", b.hit_fast.data(), "hit_fast[nhit_chan]/F");
  }

  // Branch buffers.
  struct B {
    int run = 0, event = 0, spill = 0, bcid = 0, nhit_chan = 0, nhit_slab = 0, nhit_chip = 0;
    float sum_energy = 0.f, sum_hg = 0.f;
    std::array<int, kMaxHits> hit_slab{}, hit_chip{}, hit_chan{}, hit_sca{}, hit_ismasked{};
    std::array<float, kMaxHits> hit_energy{}, hit_hg{}, hit_lg{}, hit_x{}, hit_y{}, hit_z{}, hit_X0{},
        hit_w_energy{}, hit_time{}, hit_fast{};
  };

  Gaudi::Property<std::string> m_inputName{this, "InputCollection", "SiPadHitsRealAdc",
                                           "Hits with the test-beam CellID (ChannelMapper or AdcDigitizer output)"};
  Gaudi::Property<std::string> m_maskName{this, "MaskingCollection", "",
                                          "Parallel UserDataCollection<int> of masking flags -> hit_ismasked"};
  Gaudi::Property<std::string> m_keptName{this, "KeptCollection", "",
                                          "Parallel UserDataCollection<int>: only hits flagged 1 are written"};
  Gaudi::Property<std::string> m_timeName{this, "TimeCollection", "",
                                          "Parallel UserDataCollection<float>, trigger time [ns] -> hit_time"};
  Gaudi::Property<std::string> m_fastName{this, "FastCollection", "",
                                          "Parallel UserDataCollection<float>, fast-channel peak [MIP] -> hit_fast"};
  Gaudi::Property<std::string> m_adcHighName{this, "AdcHighCollection", "",
                                             "Parallel UserDataCollection<float>, high-gain ADC -> hit_hg"};
  Gaudi::Property<std::string> m_adcLowName{this, "AdcLowCollection", "",
                                            "Parallel UserDataCollection<float>, low-gain ADC -> hit_lg"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "ecal.root", "Output ROOT file"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "ecal", "Output TTree name"};
  Gaudi::Property<int> m_runNumber{this, "RunNumber", 0, "Value of the run branch"};
  Gaudi::Property<std::string> m_bitField{this, "BitField", "system:8,slab:8,chip:16,channel:8,sca:8",
                                          "Test-beam CellID bitfield"};
  Gaudi::Property<std::vector<double>> m_wThicknesses{
      this, "WThicknesses", {2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6},
      "W absorber thickness [mm] in front of each layer (EcalPidTransformer.WThicknesses)"};
  Gaudi::Property<double> m_x0W{this, "X0W", 3.5, "Radiation length of tungsten [mm]"};

  std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_hitsHandle;
  UDHandle<int> m_maskHandle, m_keptHandle;
  UDHandle<float> m_timeHandle, m_fastHandle, m_adcHighHandle, m_adcLowHandle;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  std::vector<double> m_layerX0, m_layerWeight;
  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  mutable B m_b;
  mutable int m_event = 0;
  mutable std::mutex m_mutex;
  mutable bool m_maskWarned = false, m_keptWarned = false, m_timeWarned = false, m_fastWarned = false,
               m_adcWarned = false;
};

DECLARE_COMPONENT(EcalTreeWriter)
