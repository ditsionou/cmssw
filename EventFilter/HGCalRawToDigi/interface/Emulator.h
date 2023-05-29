#ifndef EventFilter_HGCalRawToDigi_Emulator_h
#define EventFilter_HGCalRawToDigi_Emulator_h

#include <cstddef>

#include "EventFilter/HGCalRawToDigi/interface/HGCalECONDEmulatorParameters.h"
#include "EventFilter/HGCalRawToDigi/interface/SlinkTypes.h"

namespace hgcal::econd {
  /**
     @class Emulator
     @short Pure virtual base class for a ECON-D event emulator implementation
  */
  class Emulator {
  public:
    explicit Emulator(const EmulatorParameters& params) : params_(params) {}
    virtual ~Emulator() = default;
    /// Fetch the next ECON-D event
    virtual ECONDInput next() = 0;

  protected:
    const EmulatorParameters params_;
  };

  /**
     @class EmptyEmulator
     @short an empty ECON-D payloads emulator
  */
  class EmptyEmulator : public Emulator {
  public:
    using Emulator::Emulator;

    ECONDInput next() override {
      EventId evt_id{event_id_++, bx_id_++, orbit_id_++};
      ERxData dummy_data{.cm0 = 0, .cm1 = 0, .tctp = {}, .adc = {}, .adcm = {}, .toa = {}, .tot = {}, .meta = {}};
      ERxInput empty_evt = {{ERxId_t{0 /*chip*/, 0 /*half*/}, dummy_data}};  // map<ERx_t, ERxData>
      return ECONDInput{evt_id, empty_evt};
    }

  private:
    uint32_t event_id_{1}, bx_id_{2}, orbit_id_{3};
  };

  /**
     @class TrivialEmulator
     @short A "trivial" ECON-D emulator emulating non-empty ECON-D events
  */
  class TrivialEmulator : public Emulator {
  public:
    using Emulator::Emulator;

    ECONDInput next() override;

  private:
    uint32_t event_id_{1}, bx_id_{2}, orbit_id_{3};
  };

}  // namespace hgcal::econd

#endif
