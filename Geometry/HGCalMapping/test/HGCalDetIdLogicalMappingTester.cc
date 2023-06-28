#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalMapping/interface/HGCalSiPMCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalSiCellLocator.h"
#include "Geometry/HGCalMapping/interface/HGCalModuleLocator.h"

/**
   @short a tester for the logical mapping as ESsource
 */
class HGCalDetIdLogicalMappingTester : public edm::one::EDAnalyzer<> {
public:
  explicit HGCalDetIdLogicalMappingTester(const edm::ParameterSet &iConfig)
      : geomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
        modFilename_(iConfig.getParameter<std::string>("modFile")),
        sipmFilename_(iConfig.getParameter<std::string>("sipmFile")),
        siFilename_(iConfig.getParameter<std::string>("siFile")),
        makeTextFiles_(iConfig.getParameter<bool>("makeTextFiles")),
        foundSipmCells_(0),
        foundSiCells_(0),
        missingSipmCells_(0),
        missingSiCells_(0),
        missingModules_(0) {
    sipmCellLocator_.buildLocatorFrom(sipmFilename_);
    siCellLocator_.buildLocatorFrom(siFilename_, false, true);
    modLocator_.buildLocatorFrom(modFilename_, true);

    if (makeTextFiles_) {
      missSiCellFile_.open("missingSiCells.txt");
      missSiCellFile_ << "layer u v type isHD\n";
      missModFile_.open("missingModules.txt");
      missModFile_ << "layer modU modV isSiPM\n";
      missSipmCellFile_.open("missingSipmCells.txt");
      missSipmCellFile_ << "layer iring modiphi\n";
      misMatchSipmCellFile_.open("mismatchedSipmCells.txt");
      misMatchSipmCellFile_ << "layer iring iphi originalType newType originalSipmType newSipmType originalID newId\n";
      misMatchSiCellFile_.open("mismatchedSiCells.txt");
      misMatchSiCellFile_ << "layer modu modv originalID newId\n";
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("label", {});
    desc.add<std::string>("modFile", "");
    desc.add<std::string>("sipmFile", "");
    desc.add<std::string>("siFile", "");
    desc.add<bool>("makeTextFiles", false);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  /**
     @short the main analysis method to get DetIDs from the geometry and convert them to electronic ids then back to DetIDs again.
   */
  void analyze(const edm::Event &, const edm::EventSetup &iSetup) override;

  void fillMapFromGeom(DetId::Detector detType);

  HGCalElectronicsId detIdToEleId(const DetId detId);

  DetId eIdToDetId(const HGCalElectronicsId eleId);

  //tokens and record watches
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;

  std::string modFilename_;
  std::string sipmFilename_;
  std::string siFilename_;
  bool makeTextFiles_;

  std::map<DetId, HGCalElectronicsId> mappedIds_;
  int foundSipmCells_, foundSiCells_, missingSipmCells_, missingSiCells_, missingModules_;
  std::ofstream missSiCellFile_, missModFile_, missSipmCellFile_, misMatchSiCellFile_, misMatchSipmCellFile_;

  HGCalSiPMCellLocator sipmCellLocator_;
  HGCalSiCellLocator siCellLocator_;
  HGCalModuleLocator modLocator_;
  const CaloGeometry *caloGeom_;
  const HGCalGeometry *siGeom_;
  const HGCalGeometry *eeGeom_;
};

//
void HGCalDetIdLogicalMappingTester::analyze(const edm::Event &, const edm::EventSetup &iSetup) {
  caloGeom_ = &iSetup.getData(geomToken_);
  siGeom_ = static_cast<const HGCalGeometry *>(
      caloGeom_->getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty));
  eeGeom_ = static_cast<const HGCalGeometry *>(
      caloGeom_->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty));

  fillMapFromGeom(DetId::HGCalEE);
  fillMapFromGeom(DetId::HGCalHSi);
  edm::LogInfo("HGCalMapping") << "Found " << foundSiCells_ << " si cells in geom and siCellLocator, "
                               << missingSiCells_ << " cells in geometry but not cell map." << std::endl;

  fillMapFromGeom(DetId::HGCalHSc);
  edm::LogInfo("HGCalMapping") << "Found " << foundSipmCells_ << " sipm cells in geom and sipmCellLocator, "
                               << missingSipmCells_ << " cells in geometry but not cell map." << std::endl;
  edm::LogInfo("HGCalMapping") << "Number of modules in geometry but not found in module map: " << missingModules_
                               << std::endl;

  for (const auto &eleIdIt : mappedIds_) {
    DetId newDetId = eIdToDetId(eleIdIt.second);

    if (newDetId != eleIdIt.first) {
      if (newDetId.det() == DetId::HGCalHSc) {
        HGCScintillatorDetId originalId(eleIdIt.first.rawId());
        HGCScintillatorDetId sipmDetId(newDetId.rawId());
        if (makeTextFiles_) {
          misMatchSipmCellFile_ << originalId.layer() << " " << originalId.ring() << " " << originalId.iphi() << " "
                                << originalId.type() << " " << sipmDetId.type() << " " << originalId.sipm() << " "
                                << sipmDetId.sipm() << " " << std::hex << originalId.rawId() << std::dec << " "
                                << std::hex << sipmDetId.rawId() << std::dec << std::endl;
        } else {
          edm::LogError("HGCalMapping") << "original raw ID: " << std::hex << originalId.rawId()
                                        << " does not match new raw ID: " << sipmDetId.rawId() << std::dec
                                        << " from electronics ID." << std::endl;
          edm::Exception e(
              edm::errors::NotFound,
              "HGCalDetIdLogicalMappingTester::Found a scintillator logical ID that does not map back to its detID.");
          throw e;
        }
      } else {
        HGCSiliconDetId originalId(eleIdIt.first.rawId());
        if (makeTextFiles_) {
          misMatchSiCellFile_ << originalId.layer() << " " << originalId.waferU() << " " << originalId.waferV() << " "
                              << std::hex << originalId.rawId() << std::dec << " " << std::hex << newDetId.rawId()
                              << std::dec << std::endl;
        } else {
          edm::LogError("HGCalMapping") << "original raw ID: " << std::hex << originalId.rawId()
                                        << " does not match raw ID: " << newDetId.rawId() << std::dec
                                        << " from electronics id." << std::endl;
          edm::Exception e(
              edm::errors::NotFound,
              "HGCalDetIdLogicalMappingTester::Found a silicon logical ID that does not map back to its detID.");
          throw e;
        }
      }
    }
  }
}

void HGCalDetIdLogicalMappingTester::fillMapFromGeom(DetId::Detector detType) {
  const HGCalGeometry *geom =
      static_cast<const HGCalGeometry *>(caloGeom_->getSubdetectorGeometry(detType, ForwardSubdetector::ForwardEmpty));
  for (const auto &baseCellId : geom->getValidDetIds()) {
    HGCalElectronicsId eleId;
    try {
      eleId = detIdToEleId(baseCellId);
    } catch (...) {
      continue;
    }
    if (mappedIds_.count(baseCellId) == 0) {
      mappedIds_[baseCellId] = eleId;
    } else {
      edm::Exception e(edm::errors::NotFound, "HGCalDetIdLogicalMappingTester::Found a duplicate detID in HGCAL geom.");
      throw e;
    }
  }
}

HGCalElectronicsId HGCalDetIdLogicalMappingTester::detIdToEleId(DetId detId) {
  if (detId.det() == DetId::HGCalHSi || detId.det() == DetId::HGCalEE) {
    HGCSiliconDetId siDetId(detId.rawId());
    bool isSiPM = false;
    int modU = siDetId.waferU();
    int modV = siDetId.waferV();
    int layer = siDetId.layer();
    HGCalModuleInfo modInfo;
    try {
      modInfo = modLocator_.getModuleFromGeom(layer, modU, modV, siDetId.zside(), isSiPM);
    } catch (...) {
      missingModules_++;
      if (makeTextFiles_) {
        if (siDetId.zside() < 0) {
          missModFile_ << layer << " " << modU << " " << modV << " " << isSiPM << std::endl;
        }
      }
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Failed to find module in locator");
      throw e;
    }
    HGCalParameters::waferInfo wafInfo;
    if (detId.det() == DetId::HGCalHSi) {
      wafInfo = siGeom_->topology().dddConstants().waferInfo(layer, modU, modV);
    } else {
      wafInfo = eeGeom_->topology().dddConstants().waferInfo(layer, modU, modV);
    }
    HGCalSiCellChannelInfo cellInfo;
    try {
      cellInfo = siCellLocator_.locateCellByGeom(siDetId.cellU(), siDetId.cellV(), wafInfo.part, modInfo.isHD);
    } catch (...) {
      missingSiCells_++;
      if (makeTextFiles_) {
        missSiCellFile_ << siDetId.layer() << " " << siDetId.cellU() << " " << siDetId.cellV() << " " << wafInfo.part
                        << " " << modInfo.isHD << std::endl;
      }
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Failed to find si cell in map");
      throw e;
    }
    foundSiCells_++;
    int econdErx = cellInfo.chip * 2 + cellInfo.half;

    // The electronics ID only uses 4 bits for the captureblock, so we need to skip values over 15.
    if (modInfo.captureblock > 15) {
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Module with capture block>15");
      throw e;
    }
    HGCalElectronicsId eleId(
        modInfo.zside, modInfo.fedid, modInfo.captureblock, modInfo.econdidx, econdErx, cellInfo.seq);
    return eleId;
  } else if (detId.det() == DetId::HGCalHSc) {
    bool isSiPM = true;
    HGCScintillatorDetId sipmDetId(detId.rawId());
    // Cell locator returns <layer, iring, iphi>
    std::tuple<int, int, int> modLoc = sipmCellLocator_.getModuleLocation(sipmDetId);
    int layer = std::get<0>(modLoc);
    int modiring = std::get<1>(modLoc);
    int modiphi = std::get<2>(modLoc);

    // For module locator, iu->iring, iv->iphi
    HGCalModuleInfo modInfo;
    try {
      modInfo = modLocator_.getModuleFromGeom(layer, modiring, modiphi, sipmDetId.zside(), isSiPM);
    } catch (...) {
      missingModules_++;
      if (makeTextFiles_) {
        if (sipmDetId.zside() < 0) {
          missModFile_ << layer << " " << modiring << " " << modiphi << " " << isSiPM << std::endl;
        }
      }
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Failed to find module in locator");
      throw e;
    }
    HGCalSiPMTileInfo cellInfo;
    try {
      cellInfo = sipmCellLocator_.getCellByGeom(layer, sipmDetId.ring(), (sipmDetId.iphi() - 1) % 8);
    } catch (...) {
      missingSipmCells_++;
      if (makeTextFiles_) {
        missSipmCellFile_ << layer << " " << sipmDetId.ring() << " " << (sipmDetId.iphi() - 1) % 8 << std::endl;
      }
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Failed to find sipm cell in locator");
      throw e;
    }
    foundSipmCells_++;

    // Math here requires seq<36 and halfRoc 0 or 1. econdErx = ROC*2 + halfRoc
    int econdErx = 2 * int(cellInfo.sipmcell / 72) + int((cellInfo.sipmcell % 72) / 36);
    int seqHalfRocCh = cellInfo.sipmcell % 36;

    // The electronics ID only uses 4 bits for the captureblock, so we need to skip values over 15.
    if (modInfo.captureblock > 15) {
      edm::Exception e(edm::errors::NotFound,
                       "HGCalDetIdLogicalMappingTester::detIdToEleId: Module with capture block >15");
      throw e;
    }
    const HGCalElectronicsId sipmEleId(
        modInfo.zside, modInfo.fedid, modInfo.captureblock, modInfo.econdidx, econdErx, seqHalfRocCh);
    return sipmEleId;
  } else {
    edm::Exception e(edm::errors::NotFound, "HGCalDetIdLogicalMappingTester::detIdToEleId: Invalid detector type");
    throw e;
  }
}

DetId HGCalDetIdLogicalMappingTester::eIdToDetId(const HGCalElectronicsId eleId) {
  HGCalModuleInfo modInfo = modLocator_.getModule(eleId.econdIdx(), eleId.captureBlock(), eleId.fedId());
  if (modInfo.isSiPM) {
    //Geom uses integer zside, modInfo uses boolean
    int intzside = (modInfo.zside) ? -1 : 1;

    //Can attempt to use tileInfo to get tile type and sipm type
    //HGCalParameters::tileInfo tileInfo = sipmGeom->topology().dddConstants().tileInfo(intzside,modInfo.plane,modInfo.u);
    const HGCScintillatorDetId sipmDetId =
        sipmCellLocator_.getDetId(eleId, eleId.halfrocChannel(), intzside, modInfo.plane, modInfo.u, modInfo.v);
    return sipmDetId;
  } else {
    uint8_t roc = uint8_t(eleId.econdeRx() / 2.);
    uint8_t rocHalf = uint8_t(eleId.econdeRx() % 2);

    DetId::Detector detType = DetId::HGCalEE;
    if (modInfo.plane > 26) {
      detType = DetId::HGCalHSi;
    }

    HGCalParameters::waferInfo wafInfo;
    if (detType == DetId::HGCalEE) {
      wafInfo = eeGeom_->topology().dddConstants().waferInfo(modInfo.plane, modInfo.u, modInfo.v);
    } else {
      wafInfo = siGeom_->topology().dddConstants().waferInfo(modInfo.plane, modInfo.u, modInfo.v);
    }
    HGCalSiCellChannelInfo cellInfo;
    try {
      cellInfo = siCellLocator_.locateCellByChannel(roc, rocHalf, eleId.halfrocChannel(), wafInfo.part, modInfo.isHD);
    } catch (...) {
      detType = DetId::HGCalHSi;
      wafInfo = siGeom_->topology().dddConstants().waferInfo(modInfo.plane, modInfo.u, modInfo.v);
      cellInfo = siCellLocator_.locateCellByChannel(roc, rocHalf, eleId.halfrocChannel(), wafInfo.part, modInfo.isHD);
    }

    int intzside = (modInfo.zside) ? -1 : 1;
    const HGCSiliconDetId siDetId(
        detType, intzside, wafInfo.type, modInfo.plane, modInfo.u, modInfo.v, cellInfo.iu, cellInfo.iv);
    return siDetId;
  }
}

DEFINE_FWK_MODULE(HGCalDetIdLogicalMappingTester);
