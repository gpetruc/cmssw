import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi
import RecoTracker.IterativeTracking.InitialStep_cff
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cff
highPtCkfTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cff.ckfBaseTrajectoryFilter.clone(
    ComponentName = 'highPtCkfTrajectoryFilter',
)
highPtCkfTrajectoryFilter.filterPset.minPt = 0.9
highPtCkfTrajectoryFilter.filterPset.minimumNumberOfHits = 3
lowPtCkfTrajectoryFilter = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryFilter.clone(
    ComponentName = 'lowPtCkfTrajectoryFilter',
)
lowPtCkfTrajectoryFilter.filterPset.minPt = 0.3
lowPtCkfTrajectoryFilter.filterPset.minimumNumberOfHits = 3


htCkfTrajectoryBuilder = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryBuilder.clone(
    ComponentName = cms.string('htCkfTrajectoryBuilder'),
    minNrOfHitsForRebuild = cms.int32(0),
    requireSeedHitsInRebuild = cms.bool(False),
    trajectoryFilterName = 'highPtCkfTrajectoryFilter',
)
lowPtHtCkfTrajectoryBuilder = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryBuilder.clone(
    ComponentName = cms.string('lowPtHtCkfTrajectoryBuilder'),
    minNrOfHitsForRebuild = cms.int32(0),
    requireSeedHitsInRebuild = cms.bool(False),
    trajectoryFilterName = 'lowPtCkfTrajectoryFilter',
)
htChi2Est = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    MaxChi2 = cms.double(30.0),
    nSigma = cms.double(3.0),
    ComponentName = cms.string('htChi2Est')
)

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()
pixelVerticesFake = cms.EDProducer("HTVertexCorrectionStudy",
    tracks = cms.InputTag("generalTracks"),
    trackSelection = cms.string("pt > 0.5 && numberOfValidHits > 6"),
    dz = cms.double(0.3),
    dxy = cms.double(0.0),
    phi = cms.double(0.0),
)
import RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi
ttrhbwrBootstrap = RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi.ttrhbwr.clone(
    ComputeCoarseLocalPositionFromDisk = True,
    ComponentName = 'WithTrackAngleAndBootstrap'
)    


highPtHTSeeds = cms.EDProducer("TestHT",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #vertices = cms.InputTag("pixelVerticesZero"),
    #vertices = cms.InputTag("pixelVerticesFake"),
    vertices = cms.InputTag("pixelVertices"),
    #vertices = cms.InputTag("offlinePrimaryVertices"),
    vertexSelection = cms.string("ndof > 4"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    measurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    # seeding 
    seed2d = cms.bool(False),
    seed3d = cms.bool(False),
    seedMixed = cms.bool(True),
    # binning
    etabins2d = cms.uint32(64),  #32),
    phibins2d = cms.uint32(256), #64),
    etabins3d = cms.uint32(256), #128),
    phibins3d = cms.uint32(256), #64),
    # clustering
    layerSeedCut2d = cms.uint32(4),
    layerSeedCut3d = cms.uint32(3),
    layerCut2d     = cms.uint32(5),
    layerCut3d     = cms.uint32(5),
    layerMoreCut   = cms.uint32(6),
    # pt clustering
    #ptSteps = cms.vdouble( 0.0,       2.5,      1.0,      0.7,      0.5 ), 
    #ptEdges = cms.vdouble( 3.0,-3.0,  1.5,4.0,  0.7,2.0,  0.5,1.0,  0.4,0.8 ), 
    #ptSteps = cms.vdouble(0.0, 0.8),
    ptSteps = cms.vdouble( 0.0,       2.0),
    ptEdges = cms.vdouble( 1.5,-1.5,  1.0,3.0 ), 
    autoPtStep = cms.bool(True),
    autoPtStepR   = cms.double(40.), # transverse length of trajectory which should fit in one phi bin
    autoPtStepMin = cms.double(1.),  # minimum pt to consider (included)
    # seed builder
    seedBuilderConfig = cms.PSet(
        propagator = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        #TTRHBuilder = cms.string('WithTrackAngle'), 
        TTRHBuilder = cms.string('WithTrackAngleAndBootstrap'), 
        chi2MeasurementEstimator = cms.string("htChi2Est"),
        NavigationSchool  = cms.string('SimpleNavigationSchool'),
        TrajectoryBuilder = cms.string('htCkfTrajectoryBuilder'),
        TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
        cleanTrajectoryAfterInOut = cms.bool(True),
        useHitsSplitting = cms.bool(True),
        splitSeedHits = cms.bool(True),
        # cuts for pair seeding
        detaCutPair = cms.double(0.002),
        dphiCutPair = cms.double(0.01),
        pairsOnSeedCellOnly = cms.bool(True),
        # cuts for additional hits (before and after the refit)
        dphiCutHits = cms.vdouble(0.03, 0.01),
        detaCutHits = cms.vdouble(0.04, 0.02),
        minHits = cms.uint32(4),
        maxSeedHits = cms.uint32(9),
        popOneLayer = cms.bool(True),
        keepPixelTriplets = cms.bool(True),
        # minimum number of hits in a trajectory for which ckf failed
        minHitsIfNoCkf = cms.uint32(6), 
        # max consec failing microclusters
        maxFailMicroClusters = cms.uint32(5),
        startingCovariance = cms.vdouble(
            4., # EXCEPTION: this we interpret as relative uncertainty^2 on 1/min(pt, 2 GeV),
            0.01, # dxdz
            0.01, # dydz
            625., # local x: 25cm
            625., # local y: 25cm (they get shrinked by the update)
        ),
        transientInitialStateEstimatorParameters = cms.PSet(
            propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
            numberMeasurementsForFit = cms.int32(4),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
        ),
    ),
    # for debugging
    debugger = cms.untracked.bool(False),
    tracks = cms.InputTag("trueTracks"),
)

lowPtHTSeeds = highPtHTSeeds.clone(
    # binning
    etabins2d = cms.uint32(64),
    phibins2d = cms.uint32(128),
    etabins3d = cms.uint32(256),
    phibins3d = cms.uint32(128),
    # clustering
    layerSeedCut2d = cms.uint32(4),
    layerSeedCut3d = cms.uint32(3),
    layerCut2d     = cms.uint32(5),
    layerCut3d     = cms.uint32(5),
    layerMoreCut   = cms.uint32(6),
    # 
    ptSteps = cms.vdouble(1.1, 0.90, 0.70, 0.60, 0.50),
    ptEdges = cms.vdouble(),
    autoPtStep = cms.bool(False),
    autoPtStepR   = cms.double(30.), # transverse length of trajectory which should fit in one phi bin
    autoPtStepMin = cms.double(01.),  # minimum pt to consider (included)
    # seed builder
    seedBuilderConfig = highPtHTSeeds.seedBuilderConfig.clone(
        TrajectoryBuilder = "lowPtHtCkfTrajectoryBuilder",
        # cuts for pair seeding
        detaCutPair = cms.double(0.005),
        dphiCutPair = cms.double(0.02),
        pairsOnSeedCellOnly = cms.bool(True),
        # cuts for additional hits (before and after the refit)
        dphiCutHits = cms.vdouble(0.03, 0.01),
        detaCutHits = cms.vdouble(0.04, 0.02),
 
    ),
)

import RecoTracker.TrackProducer.TrackProducer_cfi
highPtHTTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = cms.InputTag("highPtHTSeeds"),
)
lowPtHTTracks = highPtHTTracks.clone(src = "lowPtHTSeeds")

import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
highPtHTSelector = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
    src='highPtHTTracks',
    useAnyMVA = cms.bool(False),
    GBRForestLabel = cms.string('MVASelectorIter1'),
    trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
            name = 'highPtHTLoose',
            ), #end of pset
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
            name = 'highPtHTTight',
            preFilterName = 'highPtHTLoose',
            ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
            name = 'highPtHT',
            preFilterName = 'highPtHTTight',
            ),
        ) #end of vpset
    ) #end of clone
lowPtHTSelector = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
    src='lowPtHTTracks',
    useAnyMVA = cms.bool(False),
    GBRForestLabel = cms.string('MVASelectorIter1'),
    trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
            name = 'lowPtHTLoose',
            ), #end of pset
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
            name = 'lowPtHTTight',
            preFilterName = 'lowPtHTLoose',
            ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
            name = 'lowPtHT',
            preFilterName = 'lowPtHTTight',
            ),
        ) #end of vpset
    ) #end of clone

highPtHTClusterMasks = cms.EDProducer("TrackClusterRemover",
    clusterLessSolution= cms.bool(True),
    trajectories = cms.InputTag("highPtHTTracks"),
    overrideTrkQuals = cms.InputTag('highPtHTSelector','highPtHT'),
    TrackQuality = cms.string('highPurity'),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    )
)
lowPtHTSeeds.clustersToSkip = cms.InputTag("highPtHTClusterMasks")

finalHTTracks = cms.EDProducer("TrackListMerger",
    ShareFrac = cms.double(0.19),
    writeOnlyTrkQuals = cms.bool(False),
    MinPT = cms.double(0.05),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    Epsilon = cms.double(-0.001),
    #selectedTrackQuals = cms.VInputTag(cms.InputTag("highPtHTSelector","highPtHT")),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("highPtHTSelector","highPtHT"),cms.InputTag("lowPtHTSelector","lowPtHT")),
    indivShareFrac = cms.vdouble(1.0, 0.19),
    makeReKeyedSeeds = cms.untracked.bool(False),
    MaxNormalizedChisq = cms.double(1000.0),
    FoundHitBonus = cms.double(5.0),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(True),
        tLists = cms.vint32(0,1)
        #tLists = cms.vint32(0)
    )),
    MinFound = cms.int32(3),
    hasSelector = cms.vint32(1,1),
    #hasSelector = cms.vint32(1),
    TrackProducers = cms.VInputTag(cms.InputTag("highPtHTTracks"), cms.InputTag("lowPtHTTracks")),
    #TrackProducers = cms.VInputTag(cms.InputTag("highPtHTTracks")),
    LostHitPenalty = cms.double(20.0),
    newQuality = cms.string('confirmed')
)

testHT = cms.Sequence(
    #pixelVerticesZero + 
    #pixelVerticesFake +
    highPtHTSeeds * highPtHTTracks + highPtHTSelector + highPtHTClusterMasks 
    + lowPtHTSeeds * lowPtHTTracks + lowPtHTSelector 
    + finalHTTracks 
)

from PhysicsTools.RecoAlgos.trackingParticleSelector_cfi import trackingParticleSelector
tracksSim = trackingParticleSelector.clone()

finalStudy = cms.EDProducer("HTFinalOutcomeStudy",
    tracksCkf = cms.InputTag("TrackMCQuality"),
    tracksHT  = cms.InputTag("TrackMCQualityHT"),
    #tracksHT  = cms.InputTag("TrackMCQuality"),
    trackSelectionCkf = cms.string("pt > 0.5 && numberOfValidHits > 3 && quality('highPurity')"),
    trackSelectionHT  = cms.string("pt > 0.5"),
    doMC = cms.bool(False),
    tracksSim = cms.InputTag("tracksSim"),
    simAssociator = cms.string("TrackAssociatorByHits"),
)


from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.TrackMCQuality_cfi import *
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import *
TrackMCQualityHT = TrackMCQuality.clone(label_tr = "finalHTTracks")
trueTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("TrackMCQuality"),
    cut = cms.string("quality('qualitySize')"),
)
testHTMC = cms.Sequence(
    TrackMCQuality + trueTracks +
    testHT +
    TrackMCQualityHT + tracksSim + simHitTPAssocProducer
    + finalStudy
)
binningStudy = cms.EDProducer("HTBinningStudies",
    tracks = cms.InputTag("TrackMCQuality"),
    trackSelection = cms.string("pt > 0.4 && abs(eta) < 2.4 && quality('highPurity') && quality('qualitySize')"),    
)
binEfficiencyStudy = cms.EDProducer("HTBinningEfficiencyStudy",
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    tracks = cms.InputTag("TrackMCQuality"),
    trackSelection = cms.string("pt > 0.5 && abs(eta) < 2.4 && quality('highPurity') && quality('qualitySize')"),    
    useVertices = cms.bool(True),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("pixelVertices"),
    # binning
    etabins2d = cms.uint32(32),
    phibins2d = cms.uint32(64),
    etabins3d = cms.uint32(128),
    phibins3d = cms.uint32(64),
    ptSteps = cms.vdouble(0.0, 2., 1.),
    ptEdges = cms.vdouble(), 
)

binEfficiencyStudy_e128_p128 = binEfficiencyStudy.clone( etabins3d = 128, phibins3d = 128 )
binEfficiencyStudy_e128_p128_e32_p128 = binEfficiencyStudy.clone( etabins3d = 128, phibins3d = 128, etabins2d = 32, phibins2d = 128 )
binEfficiencyStudy_e128_p256 = binEfficiencyStudy.clone( etabins3d = 128, phibins3d = 256 )
binEfficiencyStudy_e128_p256_e32_p128 = binEfficiencyStudy.clone( etabins3d = 128, phibins3d = 256, etabins2d = 32, phibins2d = 128 )
binEfficiencyStudy_e256_p256 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256 )
binEfficiencyStudy_e256_p128 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 128 )
binEfficiencyStudy_e256_p256_e32_p128 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256, etabins2d = 32, phibins2d = 128 )
binEfficiencyStudy_e256_p256_e32_p256 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256, etabins2d = 32, phibins2d = 256 )
binEfficiencyStudy_e256_p256_e64_p128 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256, etabins2d = 64, phibins2d = 128 )
binEfficiencyStudy_e256_p256_e64_p256 = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256, etabins2d = 64, phibins2d = 256 )
binEfficiencyStudy_e256_p256_e64_p256_5f = binEfficiencyStudy.clone( etabins3d = 256, phibins3d = 256, etabins2d = 64, phibins2d = 256, ptSteps = [ 0, 10., 5., ] )
binStudy = cms.Sequence(
    TrackMCQuality + trueTracks +
    #binningStudy +
    binEfficiencyStudy
    #+ binEfficiencyStudy1
    #+ binEfficiencyStudy2
    #+ binEfficiencyStudy21
    #+ binEfficiencyStudy_e128_p128
    + binEfficiencyStudy_e128_p128
    + binEfficiencyStudy_e128_p128_e32_p128
    + binEfficiencyStudy_e128_p256
    + binEfficiencyStudy_e128_p256_e32_p128
    + binEfficiencyStudy_e256_p256
    + binEfficiencyStudy_e256_p256_e32_p128
    + binEfficiencyStudy_e256_p256_e64_p128
    + binEfficiencyStudy_e256_p256_e32_p256
    + binEfficiencyStudy_e256_p256_e64_p256
    #+ binEfficiencyStudy_e256_p128
)

occupancyStudy = cms.EDProducer("HTOccupancyStudy",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("pixelVertices"),
    vertexSelection = cms.string("ndof > 4"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    # seeding 
    seed2d = cms.bool(False),
    seed3d = cms.bool(False),
    seedMixed = cms.bool(True),
    # binning
    etabins2d = cms.uint32(32),
    phibins2d = cms.uint32(64),
    etabins3d = cms.uint32(128),
    phibins3d = cms.uint32(64),
    ptSteps = cms.vdouble( 0.0,  2.0, 1.5, 1.0, 0.75),
    tracks = cms.InputTag("generalTracks"),
)

occupancyStudy_e128_p128 = occupancyStudy.clone( etabins3d = 128, phibins3d = 128 )
occupancyStudy_e128_p256 = occupancyStudy.clone( etabins3d = 128, phibins3d = 256 )
occupancyStudy_e256_p256 = occupancyStudy.clone( etabins3d = 256, phibins3d = 256 )
occupancyStudy_e256_p128 = occupancyStudy.clone( etabins3d = 256, phibins3d = 128 )
occStudy = cms.Sequence(
    occupancyStudy
    + occupancyStudy_e128_p128
    + occupancyStudy_e128_p256
    + occupancyStudy_e256_p256
    + occupancyStudy_e256_p128
)


