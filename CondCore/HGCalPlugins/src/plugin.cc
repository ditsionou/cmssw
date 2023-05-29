#include "CondCore/ESSources/interface/registration_macros.h"
#include "CondFormats/DataRecord/interface/HGCalLabTestConditionsRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalLabTestConditions.h"
#include "CondFormats/DataRecord/interface/HGCalCondSerializableSiCellChannelInfoRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableSiCellChannelInfo.h"

REGISTER_PLUGIN(HGCalLabTestConditionsRcd, HGCalLabTestConditions);
REGISTER_PLUGIN(HGCalCondSerializableSiCellChannelInfoRcd, HGCalCondSerializableSiCellChannelInfo);
