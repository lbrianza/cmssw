import FWCore.ParameterSet.Config as cms

StudyHLT = cms.EDAnalyzer("StudyHLT",
                           Verbosity       = cms.untracked.int32( 0 ),
                           Triggers        = cms.untracked.vstring([]),
                           TrackQuality    = cms.untracked.string("highPurity"),
                           MinTrackPt      = cms.untracked.double(10.0),
                           MaxDxyPV        = cms.untracked.double(0.02),
                           MaxDzPV         = cms.untracked.double(0.02),
                           MaxChi2         = cms.untracked.double(5.0),
                           MaxDpOverP      = cms.untracked.double(0.1),
                           MinOuterHit     = cms.untracked.int32(4),
                           MinLayerCrossed =cms.untracked.int32(8),
                           MaxInMiss       = cms.untracked.int32(0),
                           MaxOutMiss      = cms.untracked.int32(0),
                           minTrackP       = cms.untracked.double( 1.0 ),
                           maxTrackEta     = cms.untracked.double( 2.6 ),
                           TimeMinCutECAL  = cms.untracked.double(-500.0),
                           TimeMaxCutECAL  = cms.untracked.double(500.0),
                           TimeMinCutHCAL  = cms.untracked.double(-500.0),
                           TimeMaxCutHCAL  = cms.untracked.double(500.0),
                           IsItAOD         = cms.untracked.bool(False),
)
