#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/HGCalCondSerializableModuleInfoRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableModuleInfo.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalMapping/interface/HGCalSiPMCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalSiCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalModuleLocator.h"

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
  std::vector<DetId::Detector> siDets = {DetId::HGCalEE, DetId::HGCalHSi};

  bool isSiPM;
  int foundSipmCells = 0, foundSiCells = 0, missingSipmCells = 0, missingSiCells = 0, missingModules = 0;
  for (const auto& d : siDets) {
   
    const HGCalGeometry *siGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(d, ForwardSubdetector::ForwardEmpty));
    continue;
    for (const auto& baseCellId : siGeom->getValidDetIds()) {
      HGCSiliconDetId siDetId(baseCellId.rawId());
      isSiPM = false;
      int modU = siDetId.waferU();
      int modV = siDetId.waferV();
      int layer = siDetId.layer();
      HGCalModuleInfo modInfo;
      try{
         modInfo = modLocator.getModuleFromGeom(layer, modU, modV, isSiPM);
      }
      catch(...) {
        missingModules++;
        continue;
      }

      HGCalParameters::waferInfo wafInfo = siGeom->topology().dddConstants().waferInfo(layer,modU,modV);
      HGCalSiCellChannelInfo cellInfo;
      try{ 
        cellInfo = siCellLocator.locateCellByGeom(siDetId.cellU(), siDetId.cellV(), wafInfo.part, modInfo.isHD);
      }
      catch(...) {
        missingSiCells++;
        continue;
      }
      foundSiCells++; 
      int econdErx = cellInfo.chip*2+cellInfo.half;
 
      // The electronics ID only uses 4 bits for the captureblock, so we need to skip values over 15.  
      if(modInfo.captureblock>15){continue;}

      HGCalElectronicsId siEleId(modInfo.fedid, modInfo.captureblock, modInfo.econdidx, econdErx, cellInfo.seq);
      if (mappedIds.count(siDetId)==0) {
        mappedIds[siDetId]=siEleId;
      }
      else {
         edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a duplicate detID in silicon HGCAL geom.");
         throw e;
      }
    }
    std::string detType("HGCalEE");
    if(d==DetId::HGCalHSi){detType="HGCalHSi";}
    edm::LogInfo("HGCalMapping") << "Found " << foundSiCells << " " << detType << " cells in both geometry and si cell mapper, " << missingSiCells << " included in geometry but missing in map." << std::endl;
  }

  const HGCalGeometry *sipmGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));

  for (const auto& baseCellId : sipmGeom->getValidDetIds()) {
    isSiPM = true;
    HGCScintillatorDetId sipmDetId(baseCellId.rawId());
    // Cell locator returns <layer, iring, iphi>
    std::tuple<int, int, int> modLoc = sipmCellLocator.getModuleLocation(sipmDetId);
    int layer = std::get<0>(modLoc);
    int modiring = std::get<1>(modLoc);
    int modiphi = std::get<2>(modLoc);
    
    // For module locator, iu->iring, iv->iphi
    HGCalModuleInfo modInfo;
    try{
      modInfo = modLocator.getModuleFromGeom(layer, modiring, modiphi, isSiPM);
    }
    catch(...) {
      missingModules++;
      continue;
    }
    HGCalSiPMTileInfo cellInfo;
    try{ 
      cellInfo = sipmCellLocator.getCellByGeom(layer,sipmDetId.ring(),(sipmDetId.iphi()-1)%8);
    }
    catch(...) {
      missingSipmCells++;
      continue;
    }
    foundSipmCells++;

    // Math here requires seq<36 and halfRoc 0 or 1. econdErx = ROC*2 + halfRoc
    int econdErx = 2*int(cellInfo.sipmcell/72) + int((cellInfo.sipmcell%72)/36);
    int seqHalfRocCh = cellInfo.sipmcell%36;

    // The electronics ID only uses 4 bits for the captureblock, so we need to skip values over 15.
    if(modInfo.captureblock>15){continue;}

    HGCalElectronicsId sipmEleId(modInfo.fedid, modInfo.captureblock, modInfo.econdidx, econdErx, seqHalfRocCh);
    if (mappedIds.count(sipmDetId)==0) {
      mappedIds[sipmDetId]=sipmEleId;
    }
    else {
       edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a duplicate detID in scintillator HGCAL geom.");
       throw e;
    }
  }

  edm::LogInfo("HGCalMapping") << "Found " << foundSipmCells << " in geom and siCellLocator, missing " << missingSipmCells << " in geometry but not cell map." << std::endl;
  edm::LogInfo("HGCalMapping") << "Number of modules in geometry but not found in module map: " << missingModules << std::endl;

  const HGCalGeometry *siGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty));
  const HGCalGeometry *eeGeom = static_cast<const HGCalGeometry *>(caloGeom->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty));

 
  for (const auto& eleIdIt : mappedIds) {
    HGCalElectronicsId thisId = eleIdIt.second;
    HGCalModuleInfo modInfo = modLocator.getModule(thisId.econdIdx(),thisId.captureBlock(), thisId.fedId());
    if(modInfo.isSiPM){
       HGCScintillatorDetId originalId(eleIdIt.first.rawId());
       HGCalParameters::tileInfo tileInfo = sipmGeom->topology().dddConstants().tileInfo(originalId.zside(),modInfo.plane,modInfo.u);
       HGCScintillatorDetId sipmDetId = sipmCellLocator.getDetId(thisId, thisId.halfrocChannel(), originalId.zside(), modInfo.plane, modInfo.u, modInfo.v, originalId.type(), tileInfo.sipm);

       //SipmDetId constructor doesn't have a zside method. Instead, manually set the bit in the id.
       uint32_t newRawId = sipmDetId.rawId();
       //For now, all modules have zside = -1 in the mapper, so we cheat by using the original detId
       if(originalId.zside()==1){newRawId = newRawId&0xfcffffff;}
       newRawId = (newRawId&0xff7fffff)|(originalId.sipm()<<23);
       HGCScintillatorDetId modifiedSipmDetId(newRawId); 
       if(modifiedSipmDetId!=originalId) {
          edm::LogError("HGCalMapping") << "original raw ID: "  << std::hex << originalId.rawId() << ", new raw ID: " << modifiedSipmDetId.rawId() << std::dec << std::endl;  
          edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a scintillator logical ID that does not map back to its detID.");
          throw e;
       }
    }
    else {
       HGCSiliconDetId originalId(eleIdIt.first.rawId());

       uint8_t roc = uint8_t(thisId.econdeRx()/2.);
       uint8_t rocHalf = uint8_t(thisId.econdeRx()%2);
       HGCalParameters::waferInfo wafInfo;

       if(originalId.isHE()){ 
          wafInfo = siGeom->topology().dddConstants().waferInfo(modInfo.plane,modInfo.u,modInfo.v);
       }
       else{
         wafInfo = eeGeom->topology().dddConstants().waferInfo(modInfo.plane,modInfo.u,modInfo.v);
       }
       HGCalSiCellChannelInfo cellInfo = siCellLocator.locateCellByChannel(roc, rocHalf, thisId.halfrocChannel(), wafInfo.part, modInfo.isHD);
       //For now, all modules have zside = -1 in the mapper, so we cheat by using the original detId
       //HGCSiliconDetId siDetId(DetId::HGCalHSi, modInfo.zside, wafInfo.type, modInfo.plane, modInfo.u, modInfo.v, cellInfo.iu, cellInfo.iv);
       
       // !! Need to figure out whether module is EE or HE. For now, we cheat and use the original.
       DetId::Detector detType=DetId::HGCalEE;
       if(originalId.isHE()){detType=DetId::HGCalHSi;}

       HGCSiliconDetId siDetId(detType, originalId.zside(), wafInfo.type, modInfo.plane, modInfo.u, modInfo.v, cellInfo.iu, cellInfo.iv);

       if(siDetId!=originalId) {
          edm::LogError("HGCalMapping") << "original raw ID: "  << std::hex << originalId.rawId() << ", new raw ID: " << siDetId.rawId() << std::dec << std::endl;  
          edm::Exception e(edm::errors::NotFound,"HGCalDetIdLogicalMappingTester::Found a silicon logical ID that does not map back to its detID.");
          throw e; 
       }
    }
  }

}


DEFINE_FWK_MODULE(HGCalDetIdLogicalMappingTester);
