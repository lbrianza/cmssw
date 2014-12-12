import FWCore.ParameterSet.Config as cms
from CondCore.DBCommon.CondDBSetup_cfi import *

#### Please fill with comments
# Basic tag combination for 2012 Moriond conditions 
# Laser tag is fixed: no time to improve them
# Alpha tag is fixed: no improvement given by AlphaStudies
# Start with this for IC validation
# tag = cms.string('EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'),
# this tag is the one present in Prompt2012C   
# this is the reference of any IC test for 2012C

# final tag for PPD sign-off (conditions for A+B+C rereco)
# ADCtoGeV derived from this rereco validation (many IOVs to correct the scale drift)
RerecoGlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'),
                               globaltag = cms.string('GR_P_V42B::All'),
                               toGet = cms.VPSet(
    cms.PSet(record = cms.string('EcalLaserAPDPNRatiosRcd'),
             tag = cms.string('EcalLaserAPDPNRatios_20121020_447_p1_v2'),
             connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_42X_ECAL_LAS')
             ),
    cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
             tag = cms.string('EcalLaserAlphas_EB_sic1_btcp152_EE_sic1_btcp116'),
             connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_ECAL')
             ),
    cms.PSet(record = cms.string('EcalIntercalibConstantsRcd'),
             tag = cms.string('EcalIntercalibConstants_V20121215_NLT_MeanPizPhiABC_EleABC_HR9EtaScABC_PhisymCor2_EoPCorAvg1_mixed'),
             connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_ECAL')
             ),
    cms.PSet(record = cms.string('EcalClusterLocalContCorrParametersRcd'),
             tag = cms.string('EcalClusterLocalContCorrParameters_jun2012_offline'),
             connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_ECAL')
             ),
    cms.PSet(record = cms.string('EcalADCToGeVConstantRcd'),
             tag = cms.string('EcalADCToGeVConstant_Bon_RUN2012ABC_offline_V1'),
             connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_ECAL')
             )
#     cms.PSet(record = cms.string('EBAlignmentRcd'),
#              tag = cms.string('EBAlignment_measured_v07_offline'),
#              connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_34X_ECAL')
#              ),
#     cms.PSet(record = cms.string('EEAlignmentRcd'),
#              tag = cms.string('EEAlignment_measured_v07_offline'),
#              connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_34X_ECAL')
#              ),
#     cms.PSet(record = cms.string('ESAlignmentRcd'),
#              tag = cms.string('ESAlignment_measured_v07_offline'),
#              connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_PRESHOWER')
#              )
    ),
                               BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
                               )


