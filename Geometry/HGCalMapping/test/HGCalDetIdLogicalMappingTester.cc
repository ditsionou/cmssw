#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/HGCalCondSerializableModuleInfoRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableModuleInfo.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalMapping/interface/HGCalSiPMCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalSiCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalModuleLocator.h"

#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include <iostream>

/**
   @short a tester for the logical mapping as ESsource
 */
class HGCalDetIdLogicalMappingTester : public edm::one::EDAnalyzer<> {

public:
  
  explicit HGCalDetIdLogicalMappingTester(const edm::ParameterSet& iConfig)
    : geomToken_(esConsumes<CaloGeometry,CaloGeometryRecord>()),
      modFilename_(iConfig.getParameter<std::string>("modFile")),
      sipmFilename_(iConfig.getParameter<std::string>("sipmFile")),
      siFilename_(iConfig.getParameter<std::string>("siFile")) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("label", {});
    desc.add<std::string>("modFile", "");
    desc.add<std::string>("sipmFile", "");
    desc.add<std::string>("siFile", "");
    descriptions.addWithDefaultLabel(desc);
  }

private:

  /**
     @short the main analysis method to check what has been parsed from the files
   */
  void analyze(const edm::Event&, const edm::EventSetup& iSetup) override;

  //tokens and record watches
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  //edm::ESGetToken<HGCalCondSerializableModuleInfo, HGCalCondSerializableModuleInfoRcd> moduleInfoToken_;
  std::string modFilename_;
  std::string sipmFilename_;
  std::string siFilename_;

};


//
void HGCalDetIdLogicalMappingTester::analyze(const edm::Event&, const edm::EventSetup& iSetup) {

  HGCalSiPMCellLocator sipmCellLocator;
  sipmCellLocator.buildLocatorFrom(sipmFilename_);
  HGCalSiCellLocator siCellLocator;
  siCellLocator.buildLocatorFrom(siFilename_, false, true);
  HGCalModuleLocator modLocator;
  modLocator.buildLocatorFrom(modFilename_, true);

  std::map<DetId, HGCalElectronicsId> mappedIds;

  const CaloGeometry *caloGeom = &iSetup.getData(geomToken_);
  const HGCalGeometry *sipmGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));
  const HGCalGeometry *siGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty));

  bool isSiPM;
  //for (const auto& baseCellId : siGeom->getValidDetIds()) {
  //  HGCSiliconDetId siDetId(baseCellId.rawId());
  //  isSiPM = false;
  //  int modU = siDetId.waferU();
  //  int modV = siDetId.waferV();
  //  int layer = siDetId.layer();
  //  int econdIdx = modLocator.getEcondIdx(layer, modU, modV, isSiPM);
  //  int captureblock = modLocator.getCaptureBlockIdx(layer, modU, modV, isSiPM);
  //  int fedId = modLocator.getFedId(layer, modU, modV, isSiPM);

  //  HGCalSiCellChannelInfo cellInfo = siCellLocator.locateCellByGeom(siDetId.cellU(), siDetId.cellV(), siDetId.type(), false); //isHD
  //  uint16_t rocPin = cellInfo.rocpin;
  //  uint8_t rocHalf = cellInfo.half; 
  //  int rocNumber = HGCSiliconDetIdToROC().getROCNumber(siDetId.triggerCellU(), siDetId.triggerCellV(), siDetId.type())-1;
  //  int econdErx = rocNumber*2+rocHalf;
  //  HGCalElectronicsId siEleId(fedId, captureblock, econdIdx, econdErx, rocPin);
  //  if (mappedIds.count(siDetId)==0) {
  //    mappedIds[siDetId]=siEleId;
  //  }
  //  else {
  //     edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a duplicate detID in silicon HGCAL geom.");
  //     throw e;
  //  }
  //}

  int goodCells = 0;
  int missingCells = 0;
  for (const auto& baseCellId : sipmGeom->getValidDetIds()) {
    isSiPM = true;
    HGCScintillatorDetId sipmDetId(baseCellId.rawId());
    // Cell locator returns <layer, iring, iphi>
    std::tuple<int, int, int> modLoc = sipmCellLocator.getModuleLocation(sipmDetId);
    int layer = std::get<0>(modLoc);
    int modiring = std::get<1>(modLoc);
    int modiphi = std::get<2>(modLoc);
    
    // For module locator, iu->iring, iv->iphi
    int econdIdx = modLocator.getEcondIdx(layer, modiring, modiphi, isSiPM);
    int captureBlock = modLocator.getCaptureBlockIdx(layer, modiring, modiphi, isSiPM);
    int fedId = modLocator.getFedId(layer, modiring, modiphi, isSiPM);
    HGCalSiPMTileInfo cellInfo;
    try{ 
      cellInfo = sipmCellLocator.getCellByGeom(layer,sipmDetId.ring(),(sipmDetId.iphi()-1)%8);
    }
    catch(...) {
      missingCells++;
      continue;
    }
    // Math here requires seq<36 and halfRoc 0 or 1. econdErx = ROC*2 + halfRoc
    int econdErx = 2*int(cellInfo.sipmcell/72) + int((cellInfo.sipmcell%72)/36);
    int halfRocCh = cellInfo.sipmcell%36;
    goodCells++;
    HGCalElectronicsId sipmEleId(fedId, captureBlock, econdIdx, econdErx, halfRocCh);
    if (mappedIds.count(sipmDetId)==0) {
      mappedIds[sipmDetId]=sipmEleId;
    }
    else {
       edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a duplicate detID in scintillator HGCAL geom.");
       throw e;
    }
  }

  std::cout << "Good sipm cells: " << goodCells << ", bad sipm cells: " << missingCells << std::endl;
  
  /*for (const auto& eleIdIt : mappedIds) {
    HGCalElectronicsId thisId = eleIdIt.second;
    HGCalModuleInfo modInfo = modLocator.getModule(thisId.econdIdx(),thisId.captureBlock(), thisId.fedId());
    if(modInfo.isSiPM){
       HGCScintillatorDetId sipmDetId = sipmCellLocator.getDetId(thisId, thisId.halfrocChannel(), modInfo.zside, modInfo.plane, modInfo.u, modInfo.v);
       HGCScintillatorDetId originalId(eleIdIt.first.rawId());
       //if(sipmDetId!=originalId) {
       //   std::cout << "siPM det id layer: " << originalId.layer() << ", ring: " << originalId.ring() << ", iphi: " << originalId.iphi() << ", type: " << originalId.type() << ", sipm: " << originalId.sipm() << ", trigger: " << originalId.trigger() << ", zside: " << originalId.zside() << std::endl;
       //   std::cout << "log id det id layer: " << sipmDetId.layer() << ", ring: " << sipmDetId.ring() << ", iphi: " << sipmDetId.iphi() << ", type: " << sipmDetId.type() << ", sipm: " << sipmDetId.sipm() << ", trigger: " << originalId.trigger() << ", zside: " << sipmDetId.zside() << std::endl;
       //   edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a scintillator logical ID that does not map back to its detID.");
       //   throw e;
       //}
    }
    else {
       HGCSiliconDetId siDetId = siCellLocator.getDetId(thisId, modInfo.zside, modInfo.plane, modInfo.u, modInfo.v);
       if(siDetId!=eleIdIt.first) {
          edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a scintillator logical ID that does not map back to its detID.");
          throw e; 
       }
    }
  }*/

}


DEFINE_FWK_MODULE(HGCalDetIdLogicalMappingTester);
