import FWCore.ParameterSet.Config as cms

import L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClusters_cfi

pfClustersFromHGC3DClustersEM = L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClusters_cfi.pfClustersFromHGC3DClusters.clone(
    emOnly = cms.bool(True),
    etMin = cms.double(0.0), 
    corrector  = cms.string("L1Trigger/Phase2L1ParticleFlow/data/emcorr_hgc.root"),
    preEmId  = cms.string(""),
    resol = cms.PSet(
            etaBins = cms.vdouble( 1.900,  2.200,  2.500,  2.800,  2.950),
            offset  = cms.vdouble( 2.059,  1.310,  0.786,  0.490,  0.246),
            scale   = cms.vdouble( 0.034,  0.021,  0.020,  0.021,  0.046),
            kind    = cms.string('calo'),
    )
)
