#ifndef EventFilter_HGCalRawToDigi_HGCalRawDataBaseEmulator_h
#define EventFilter_HGCalRawToDigi_HGCalRawDataBaseEmulator_h

#include <cstddef>
#include <vector>
#include "EventFilter/HGCalRawToDigi/interface/HGCalECONDEmulatorParameters.h"
#include "EventFilter/HGCalRawToDigi/interface/ECONDTypes.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

namespace hgcal {

  /**
     @class template class for raw data emulation
  */
  template<class T, class C>
  class Emulator {
  public:
    explicit Emulator(const C& params) : params_(params) {}
    virtual ~Emulator() = default;

    /// Fetch the next ECON-D event
    virtual T next() = 0;

    /// get the configuration parameters
    const C &config() { return params_; }
    
  protected:
    const C params_;
  };


  typedef Emulator<hgcal::econd::ECONDInput,hgcal::econd::EmulatorParameters> ECONDEmulatorBase;
  typedef Emulator<FEDRawDataCollection, edm::ParameterSet> SlinkEmulatorBase;
  
}  // namespace hgcal

#endif
