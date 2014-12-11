import FWCore.ParameterSet.Config as cms

        
from RecoLocalCalo.Configuration.RecoLocalCalo_cff import *
#ecalRecHit.EBuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB","ALCASKIM")
#ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEE","ALCASKIM")
ecalLocalRecoSeq = cms.Sequence( ecalRecHit * ecalCompactTrigPrim * ecalTPSkim + ecalPreshowerRecHit)

from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
#pfClusteringECAL = cms.Sequence(particleFlowRecHitECAL*particleFlowClusterECAL)
rerecoPFClusteringSeq = cms.Sequence(pfClusteringPS + pfClusteringECAL)

from  RecoEcal.Configuration.RecoEcal_cff import *

from Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi import *
#electronRecalibSCAssociator.scIslandCollection = cms.string('endcapRecalibSC')
#electronRecalibSCAssociator.scIslandProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower')
#electronRecalibSCAssociator.scProducer = cms.string('correctedHybridSuperClusters')
#electronRecalibSCAssociator.scCollection = cms.string('recalibSC')
#electronRecalibSCAssociator.electronProducer = 'gedGsfElectrons'
electronClusteringSeq = cms.Sequence(ecalClusters * electronRecalibSCAssociator)


sandboxRerecoSeq = cms.Sequence(ecalLocalRecoSeq * electronClusteringSeq)
sandboxPFRerecoSeq = cms.Sequence(ecalLocalRecoSeq * rerecoPFClusteringSeq * electronClusteringSeq)
