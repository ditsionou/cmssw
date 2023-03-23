#ifndef _geometry_hgcalmapping_hgcalmappingrcd_h_
#define _geometry_hgcalmapping_hgcalmappingrcd_h_

#include "FWCore/Utilities/interface/mplVector.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "CondFormats/DataRecord/interface/HGCalMappingElectronicsRcd.h"

class HGCalMappingRcd
    : public edm::eventsetup::DependentRecordImplementation<HGcalMappingRcd,
                                                            edm::mpl::Vector<HGCalMappingElectronicsRcd> > {};

#endif
