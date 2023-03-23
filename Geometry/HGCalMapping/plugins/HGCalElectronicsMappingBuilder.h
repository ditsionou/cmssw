#ifndef Geometry_HGCalMapping_HGCalElectronicsMappingBuilder
#define Geometry_HGCalMapping_HGCalElectronicsMappingBuilder

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "Geometry/HGCalMapping/interface/HGCalMappingRcd.h"
#include "Geometry/HGCalMapping/interface/HGCalSiCellLocator.h"
#include "CondFormats/HGCalObjects/interface/HGCalCondObjectContainer.h"


class HGCalElectronicsMappingBuilder : public edm::ESProducer {
public:
  HGCalElectronicsMappingBuilder(const edm::ParameterSet&);

  using ReturnType = std::unique_ptr<HGCalCondUInt32Container>;

  ReturnType produce(const HGCalMappingRcd&);

private:

  edm::ESGetToken<HGCalCondUInt32Container, HGCalMappingElectronicsRcd> token_;
};
#endif
