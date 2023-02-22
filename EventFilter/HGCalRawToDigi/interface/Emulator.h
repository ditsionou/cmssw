#ifndef EventFilter_HGCalRawToDigi_Emulator_h
#define EventFilter_HGCalRawToDigi_Emulator_h

#include "EventFilter/HGCalRawToDigi/interface/SlinkTypes.h"

namespace hgcal::econd {
  class Emulator {
  public:
    explicit Emulator(size_t num_channels) : num_channels_(num_channels) {}
    virtual ~Emulator() = default;

    virtual ECONDEvent next() = 0;

  protected:
    const size_t num_channels_;
  };

  class TrivialEmulator : public Emulator {
  public:
    using Emulator::Emulator;

    ECONDEvent next() override {
      EventId evt_id{event_id_++, bx_id_++, orbit_id_++};
      ERxData dummy_data{.adc = {}, .adcm = {}, .toa = {}, .tot = {}, .tctp = {}, .cm0 = 0, .cm1 = 0};
      ERxEvent evt = {// map<ERx_t, ERxData>
                      {ERx_t{
                           0,  // chip
                           0   // half
                       },
                       dummy_data}};
      return std::make_pair(evt_id, evt);
    }

  private:
    uint32_t event_id_{1}, bx_id_{2}, orbit_id_{3};
  };
}  // namespace hgcal::econd

#endif
