#include "HGCalElectronicsMappingBuilder.h"
//#include "DataFormats/HGCalDetId/interface/EEDetId.h"   //fixme add Si,Scintillator/LogicalId
//#include "DataFormats/HGCalDetId/interface/HGCalElectronicsId.h"
//#include "DataFormats/HGCalDetId/interface/HGCalTriggerElectronicsId.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "CondFormats/DataRecord/interface/HGCalMappingElectronicsRcd.h" //FIXME needs to be created

HGCalElectronicsMappingBuilder::HGCalElectronicsMappingBuilder(const edm::ParameterSet&)
    : token_{
             setWhatProduced(this).consumesFrom<HGCalMappingElectronics, HGCalMappingElectronicsRcd>(edm::ESInputTag{})} {}

//
HGCalElectronicsMappingBuilder::ReturnType HGCalElectronicsMappingBuilder::produce(const HGCalMappingRcd& iRecord) {

  auto prod = std::make_unique<HGCalElectronicsMapping>();

  const auto& item = iRecord.get(token_);
  const std::vector<HGCalMappingElement> &coll = item.items();

  //loop over valid DetIds
  //  if det2elec do something
  //  if elec2det do the opposite
  //  prod->assign(rawId1,rawId2);
  
  return prod;
}
