import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_forPhase2 as _miniAOD_forPhase2
modifyPhysicsToolsPatAlgosSlimmingMiniAODToolsForPhase2_ = eras.phase2_common.makeProcessModifier( _miniAOD_forPhase2 )
