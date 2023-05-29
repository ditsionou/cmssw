#include <algorithm>

#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableSiCellChannelInfo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

//
HGCalCondSerializableSiCellChannelInfo& HGCalCondSerializableSiCellChannelInfo::addParameter(HGCalSiCellChannelInfo& info) {
  params_.push_back(info);
  return *this;
}
