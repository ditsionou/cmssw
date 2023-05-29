import FWCore.ParameterSet.Config as cms

hgcconfigESSourceFromYAML = cms.ESProducer("HGCalConfigESSourceFromYAML",
                                           filename = cms.untracked.string('/eos/cms/store/group/dpg_hgcal/tb_hgcal/2023/labtest/module822/pedestal_run0.yaml')
)
