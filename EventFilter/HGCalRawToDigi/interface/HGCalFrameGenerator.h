/****************************************************************************
 *
 * This is a part of HGCAL offline software.
 * Authors:
 *   Laurent Forthomme
 *
 ****************************************************************************/

#ifndef EventFilter_HGCalRawToDigi_HGCalFrameGenerator_h
#define EventFilter_HGCalRawToDigi_HGCalFrameGenerator_h

#include "DataFormats/HGCalDigi/interface/HGCalRawDataEmulatorInfo.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalECONDEmulatorParameters.h"
#include "EventFilter/HGCalRawToDigi/interface/SlinkTypes.h"

#include <cstdint>
#include <vector>

namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm
namespace CLHEP {
  class HepRandomEngine;
}

namespace hgcal {
  namespace econd {
    class Emulator;
  }
  /// A S-link/ECON-D payload generator helper
  class HGCalFrameGenerator {
  public:
    explicit HGCalFrameGenerator(const edm::ParameterSet&);

    static edm::ParameterSetDescription description();

    void setRandomEngine(CLHEP::HepRandomEngine& rng);
    void setEmulator(econd::Emulator&);

    /// Produce a S-link event from an emulated event
    std::vector<uint64_t> produceSlinkEvent(uint32_t fed_id) const;
    const HGCalSlinkEmulatorInfo& lastSlinkEmulatedInfo() const { return last_slink_emul_info_; }

    /// Produce a ECON-D event from an emulated event
    std::vector<uint32_t> produceECONEvent(unsigned int) const;
    const HGCalECONDEmulatorInfo& lastECONDEmulatedInfo() const { return last_econd_emul_info_; }
    const econd::ECONDInput& lastECONDEmulatedInput() const { return last_emul_event_; }

    struct SlinkParameters {
      std::vector<unsigned int> active_econds{};
      unsigned int boe_marker{0}, eoe_marker{0}, format_version{0}, num_capture_blocks{1};
      bool store_header_trailer{true};
    };
    const SlinkParameters& slinkParams() const { return slink_; }
    const econd::EmulatorParameters& econdParams(unsigned int) const;

  private:
    std::vector<bool> generateEnabledChannels(unsigned int) const;
    std::vector<uint32_t> generateERxData(unsigned int, const econd::ERxInput&) const;

    static constexpr size_t max_num_econds_ = 12;

    struct HeaderBits {
      bool bitO, bitB, bitE, bitT, bitH, bitS;
    };
    HeaderBits generateStatusBits(unsigned int) const;
    /// 32bit CRC
    uint32_t computeCRC(const std::vector<uint32_t>&) const;

    SlinkParameters slink_;
    std::map<unsigned int, econd::EmulatorParameters> econd_params_;

    CLHEP::HepRandomEngine* rng_{nullptr};    // NOT owning
    mutable econd::Emulator* emul_{nullptr};  // NOT owning

    mutable HGCalECONDEmulatorInfo last_econd_emul_info_;
    mutable HGCalSlinkEmulatorInfo last_slink_emul_info_;
    mutable econd::ECONDInput last_emul_event_;
  };
}  // namespace hgcal

#endif
