//
// Creation date: 2017/08/09
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisFilterTrees.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TRandom.h>
#include <TF1.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisFilterTrees);


//___________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fMixingHandlerMult(),
  fMultBinsMixing(),
  fNMultBinsMixing(0),
  fEventCuts(),
  fTrackCuts(),
  fWriteFilteredTracks(kTRUE),
  fWriteFilteredTracksCandidatesOnly(kFALSE),
  fOptionRunMixing(kTRUE),
  fOptionRunMixingMult(kFALSE),
  fComputeMult(kTRUE),
  fReweightParticleComposition(kFALSE),
  fPairCuts(),
  fWriteFilteredPairs(kTRUE),
  fRejectEmptyEvents(kFALSE),
  fBuildCandidatePairs(kFALSE),
  fBuildCandidateLikePairs(kFALSE),
  fMCTruthJpsi2eeOnly(kFALSE),
  fRegionsToMCTruth(kFALSE),
  fDefaultRandomPhi(kTRUE),
  fMinPtLeading(0.0),
  fReweightCut(-1),
  fCandidateType(AliReducedPairInfo::kJpsiToEE),
  fLeg1Cuts(),
  fLeg2Cuts(),
  fCandidatePairCuts(),
  fRunCandidatePrefilter(kFALSE),
  fRunCandidatePrefilterOnSameCharge(kFALSE),
  fJpsiMassDist(),
  fMCParticleCompositionWeights(),
  fMeasuredMultTrackCuts(),
  fTrueMultTrackCuts(),
  fLeg1PrefilterCuts(),
  fLeg2PrefilterCuts(),
  fLeg1PairPrefilterCuts(),
  fLeg2PairPrefilterCuts(),
  fLeg1Tracks(),
  fLeg2Tracks(),
  fLeg1PrefilteredTracks(),
  fLeg2PrefilteredTracks(),
  fOptionRunOverMC(kFALSE),
  fLegCandidatesMCcuts(),
  fJpsiMotherMCcuts(),
  fMCJpsiPtWeights(0x0),
  fSkipMCEvent(kFALSE),
  fJpsiElectronMCcuts()
{
  //
  // default constructor
  //
}

//___________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fMixingHandlerMult(),  
  fMultBinsMixing(),
  fNMultBinsMixing(0),  
  fEventCuts(),
  fTrackCuts(),
  fWriteFilteredTracks(kTRUE),
  fWriteFilteredTracksCandidatesOnly(kFALSE),
  fRejectEmptyEvents(kFALSE),
  fOptionRunMixing(kTRUE),
  fOptionRunMixingMult(kFALSE),
  fComputeMult(kTRUE),
  fReweightParticleComposition(kFALSE),
  fPairCuts(),
  fWriteFilteredPairs(kTRUE),
  fBuildCandidatePairs(kFALSE),
  fBuildCandidateLikePairs(kFALSE),
  fMCTruthJpsi2eeOnly(kFALSE),
  fRegionsToMCTruth(kFALSE),
  fDefaultRandomPhi(kTRUE),
  fMinPtLeading(0.0),
  fReweightCut(-1),
  fCandidateType(AliReducedPairInfo::kJpsiToEE),
  fLeg1Cuts(),
  fLeg2Cuts(),
  fCandidatePairCuts(),
  fRunCandidatePrefilter(kFALSE),
  fRunCandidatePrefilterOnSameCharge(kFALSE),
  fJpsiMassDist(),
  fMCParticleCompositionWeights(),
  fMeasuredMultTrackCuts(),
  fTrueMultTrackCuts(),
  fLeg1PrefilterCuts(),
  fLeg2PrefilterCuts(),
  fLeg1PairPrefilterCuts(),
  fLeg2PairPrefilterCuts(),
  fLeg1Tracks(),
  fLeg2Tracks(),
  fLeg1PrefilteredTracks(),
  fLeg2PrefilteredTracks(),
  fOptionRunOverMC(kFALSE),
  fLegCandidatesMCcuts(),
  fJpsiMotherMCcuts(),
  fMCJpsiPtWeights(0x0),
  fSkipMCEvent(kFALSE),
  fJpsiElectronMCcuts()  
{
  //
  // named constructor
  //
   fMixingHandlerMult.SetOwner(kTRUE);
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
   fLeg1Cuts.SetOwner(kTRUE);
   fLeg2Cuts.SetOwner(kTRUE);
   fCandidatePairCuts.SetOwner(kTRUE);
   fMeasuredMultTrackCuts.SetOwner(kTRUE);
   fTrueMultTrackCuts.SetOwner(kTRUE);
   fLeg1PrefilterCuts.SetOwner(kTRUE);
   fLeg2PrefilterCuts.SetOwner(kTRUE);
   fLeg1PairPrefilterCuts.SetOwner(kTRUE);
   fLeg2PairPrefilterCuts.SetOwner(kTRUE);
   fLeg1Tracks.SetOwner(kFALSE);
   fLeg2Tracks.SetOwner(kFALSE);
   fLeg1PrefilteredTracks.SetOwner(kFALSE);
   fLeg2PrefilteredTracks.SetOwner(kFALSE);
}

//___________________________________________________________________________
AliReducedAnalysisFilterTrees::~AliReducedAnalysisFilterTrees() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); fTrackCuts.Clear("C"); fPairCuts.Clear("C");
   fLeg1Cuts.Clear("C"); fLeg2Cuts.Clear("C"); fCandidatePairCuts.Clear("C");
   fMeasuredMultTrackCuts.Clear("C"); fTrueMultTrackCuts.Clear("C");
   fLeg1PrefilterCuts.Clear("C"); fLeg2PrefilterCuts.Clear("C");
   fLeg1PairPrefilterCuts.Clear("C"); fLeg2PairPrefilterCuts.Clear("C");
   fLeg1Tracks.Clear("C"); fLeg2Tracks.Clear("C");
   fLeg1PrefilteredTracks.Clear("C"); fLeg2PrefilteredTracks.Clear("C");
   if(fHistosManager) delete fHistosManager;
   if(fMixingHandler) delete fMixingHandler;
   if (fJpsiMassDist) delete fJpsiMassDist;
   fMixingHandlerMult.Clear("C");

}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
   fMixingHandler->SetHistogramManager(fHistosManager);
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Process() {
  //
  // process the current event
  //  

  if(!fEvent) return;
  if(!(fEvent->IsA()==AliReducedEventInfo::Class())) {
     cout << "ERROR: AliReducedAnalysisFilterTrees::Process() needs AliReducedEventInfo events" << endl;
     return;
  }
 
  if(fEventCounter%10000==0) 
     cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;
  AliReducedVarManager::SetEvent(fEvent);

  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: the run wise data will be updated automatically in the VarManager in case a run number change is detected
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i) fValues[i]=-9999.;
  
  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);


  if(fComputeMult) FillMultiplicity(kFALSE);
  int nGlobalEstimators = (fComputeMult ? GetNMeasMultCuts() : 0) ;


  bool isEventUnbiased = fEvent->TestEventTag(14) || GetRunOverMC(); 
  //Assuming that we have selected 2% unbiased events only for data and not for MC

  if(isEventUnbiased){
      fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
      for(UShort_t ibit=0; ibit<64; ++ibit) {
         AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
         fHistosManager->FillHistClass("EventTag_BeforeCuts", fValues);
      }
      for(UShort_t ibit=0; ibit<64; ++ibit) {
            AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
            fHistosManager->FillHistClass("EventTriggers_BeforeCuts", fValues);
      }
  }


  bool isEventSelected = IsEventSelected(fEvent, fValues);

   // Trigger efficiency
   if (isEventUnbiased && isEventSelected && GetRunOverMC() && (fValues[AliReducedVarManager::kMCNch] > 0)) {
      fHistosManager->FillHistClass("Event_INELGT0", fValues);
      for (int icut = 0; icut < nGlobalEstimators; icut++) {
         fHistosManager->FillHistClass(Form("EventMult_%s_INELGT0", GetMeasMultcutName(icut)), fValues);
      }
      if (fValues[AliReducedVarManager::kINT7Triggered]) {
         fHistosManager->FillHistClass("Event_INELGT0_MB", fValues);
         for (int icut = 0; icut < nGlobalEstimators; icut++) {
            fHistosManager->FillHistClass(Form("EventMult_%s_INELGT0_MB", GetMeasMultcutName(icut)), fValues);
         }
      }
   }
   if (GetRunOverMC() && !(fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered])) isEventSelected = kFALSE;

   //Vtx reconstruction
   if(isEventUnbiased && isEventSelected) {
      fHistosManager->FillHistClass("Event_NoVtxRec", fValues);
      for (int icut = 0; icut < nGlobalEstimators; icut++) {
         if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered]) fHistosManager->FillHistClass(Form("EventMult_%s_NoVtxRec", GetMeasMultcutName(icut)), fValues);
      }
   }
   if (fValues[AliReducedVarManager::kNVtxContributors] < 0.1) isEventSelected = kFALSE;

   //Vtx z cut
   if(isEventUnbiased && isEventSelected) {
      fHistosManager->FillHistClass("Event_NoVtxzCut", fValues);
      for (int icut = 0; icut < nGlobalEstimators; icut++) {
         if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered]) fHistosManager->FillHistClass(Form("EventMult_%s_NoVtxzCut", GetMeasMultcutName(icut)), fValues);
      }
   }
   if (abs(fValues[AliReducedVarManager::kVtxZ]) > 10) isEventSelected = kFALSE;


  if(!isEventSelected) {
      for(int v=AliReducedVarManager::kMultiplicity; v<AliReducedVarManager::kMultiplicity + AliReducedVarManager::kNMultiplicityEstimators;v++)
         fValues[v]=-999. ; //Reset the multiplicity estimators
      fValues[AliReducedVarManager::kMCNch09 + 2]=-999.;
   }


  if (isEventUnbiased) {
      for (int cutMode = 0; cutMode < 4*nGlobalEstimators; cutMode = cutMode+4) {
         // For multiplicity unfolding
         // For MC, only the smearing matrix is important (supposed independent of the trigger) 
         // => HM filled even if no HM triggered events

         // MB triggered
         if (fValues[AliReducedVarManager::kINT7Triggered] || GetRunOverMC()) {
            fHistosManager->FillHistClass(Form("pp_13TeV_Data_cutMode_%d", cutMode+100),fValues);
            fHistosManager->FillHistClass(Form("pp_13TeV_MC_cutMode_%d", cutMode+100),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_Data_cutMode_%d", cutMode+100),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d", cutMode+100),fValues);
         }

         // HM triggered
         if (fValues[AliReducedVarManager::kHighMultV0Triggered] || GetRunOverMC()) {
            fHistosManager->FillHistClass(Form("pp_13TeV_Data_cutMode_%d", cutMode+101),fValues);
            fHistosManager->FillHistClass(Form("pp_13TeV_MC_cutMode_%d", cutMode+101),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_Data_cutMode_%d", cutMode+101),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d", cutMode+101),fValues);
         }

         // Inclusive (MB or HM triggered)
         if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered]) {
            fHistosManager->FillHistClass(Form("pp_13TeV_Data_cutMode_%d", cutMode+102),fValues);
            fHistosManager->FillHistClass(Form("pp_13TeV_MC_cutMode_%d", cutMode+102),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_Data_cutMode_%d", cutMode+102),fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d", cutMode+102),fValues);
         }

         // Fill for each jpsi counts
         for (int j = 0; j < fValues[AliReducedVarManager::kMCNJpsi]; j++) {
            fHistosManager->FillHistClass(Form("pp_13TeV_MC_cutMode_%d", cutMode+103), fValues);
            fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d", cutMode+103), fValues);
         }
      }
  }

  // apply event selection
  if(!isEventSelected) return;

  // Temporary stuff for MC to check pion pT distribution - should be removed soon
   /*AliReducedBaseTrack* track = 0x0;
   TClonesArray* trackList = fEvent->GetTracks2();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedBaseTrack*)nextTrack();
      //if (track->IsMCTruth() && track->TestMCFlag(14) && abs(track->Rapidity(0.13957)) < 0.5) {
      if (track->IsMCTruth() && track->TestMCFlag(7) && abs(track->Rapidity(0.13957)) < 0.5) {
         AliReducedVarManager::FillTrackInfo(track, fValues);
         fHistosManager->FillHistClass("PionPrimary", fValues);
      }
   }
return;*/

  if(fOptionRunOverMC) {
   fSkipMCEvent = kFALSE;
   FillMCTruthHistograms();
   if(fSkipMCEvent) return;
  }

  // fill event info histograms after cuts
  if(isEventUnbiased){
   fHistosManager->FillHistClass("Event_AfterCuts", fValues);
   for (int icut = 0; icut<nGlobalEstimators; icut++) {
         if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered] /*|| fValues[AliReducedVarManager::kHighMultSPDTriggered]*/) { 
            fHistosManager->FillHistClass(Form("EventMult_%s_Inclusive",GetMeasMultcutName(icut)), fValues);
         }
         if (fValues[AliReducedVarManager::kINT7Triggered]) fHistosManager->FillHistClass(Form("EventMult_%s_MB",GetMeasMultcutName(icut)), fValues);
         if (fValues[AliReducedVarManager::kHighMultV0Triggered] ) fHistosManager->FillHistClass(Form("EventMult_%s_HM",GetMeasMultcutName(icut)), fValues);
         if (fValues[AliReducedVarManager::kHighMultSPDTriggered] ) fHistosManager->FillHistClass(Form("EventMult_%s_HighMultSPD",GetMeasMultcutName(icut)), fValues);
   }

   if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered] /*|| fValues[AliReducedVarManager::kHighMultSPDTriggered]*/) fHistosManager->FillHistClass("CorrelMult_Inclusive", fValues); //Correlations between different multiplicity estimators
   if (fValues[AliReducedVarManager::kINT7Triggered]) fHistosManager->FillHistClass("CorrelMult_MB", fValues); //Correlations between different multiplicity estimators
   if (fValues[AliReducedVarManager::kHighMultV0Triggered]) fHistosManager->FillHistClass("CorrelMult_HM", fValues); //Correlations between different multiplicity estimators

   for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_AfterCuts", fValues);
   }
   for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_AfterCuts", fValues);
   }
  }

  CreateFilteredEvent();
  if(fRejectEmptyEvents && (fFilteredEvent->NPairs()+fFilteredEvent->NTracks1()+fFilteredEvent->NTracks2())==0)
      return;
  fFilteredTree->Fill();

}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::CreateFilteredEvent() {
   //
   // create the filtered event
   //
   // NOTE: The folowing information is filtered from the input event 
   //     1) Event header (either base event or full event header)
   //     2) selected V0 candidates -> see WriteFilteredPairs()
   //     3) selected tracks and legs of selected V0 candidates -> see WriteFilteredTracks()
   //     4) create candidates of the types specified in AliReducedPairInfo -> see BuildCandidatePairs()
   if(fFilteredTreeWritingOption==kFullEventsWithBaseTracks || fFilteredTreeWritingOption==kFullEventsWithFullTracks)
      ((AliReducedEventInfo*)fFilteredEvent)->CopyEventHeader((AliReducedEventInfo*)fEvent);
   else
      fFilteredEvent->CopyEventHeader(fEvent);
   
   if(fWriteFilteredPairs) WriteFilteredPairs();
   if(fBuildCandidatePairs) BuildCandidatePairs();
   if(fWriteFilteredTracks) WriteFilteredTracks();
   if(fWriteFilteredTracks && !fWriteFilteredTracksCandidatesOnly) WriteFilteredTracks(2);

   if(fComputeMult) {
      FillMultiplicity(kTRUE);
      for (int icut = 0; icut<GetNMeasMultCuts(); icut++) {

         if(fEvent->TestEventTag(14) || GetRunOverMC()) { //event is unbiased
            // Both triggers + MB + HM  (regions to jpsi)
            if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered] || fValues[AliReducedVarManager::kHighMultSPDTriggered]) fHistosManager->FillHistClass(Form("EventMultRegions_%s_Inclusive",GetMeasMultcutName(icut)), fValues);
            if (fValues[AliReducedVarManager::kINT7Triggered]) fHistosManager->FillHistClass(Form("EventMultRegions_%s_MB",GetMeasMultcutName(icut)), fValues);
            if (fValues[AliReducedVarManager::kHighMultV0Triggered] ) fHistosManager->FillHistClass(Form("EventMultRegions_%s_HM",GetMeasMultcutName(icut)), fValues);
            if (fValues[AliReducedVarManager::kHighMultSPDTriggered] ) fHistosManager->FillHistClass(Form("EventMultRegions_%s_HighMultSPD",GetMeasMultcutName(icut)), fValues);
            
            if(fValues[AliReducedVarManager::kPtLeading + icut] > fMinPtLeading) {
               // Both triggers + MB + HM  (regions to leading)
               if (fValues[AliReducedVarManager::kINT7Triggered] || fValues[AliReducedVarManager::kHighMultV0Triggered] || fValues[AliReducedVarManager::kHighMultSPDTriggered]) fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s_Inclusive",GetMeasMultcutName(icut)), fValues);
               if (fValues[AliReducedVarManager::kINT7Triggered]) fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s_MB", GetMeasMultcutName(icut)), fValues);
               if (fValues[AliReducedVarManager::kHighMultV0Triggered] ) fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s_HM", GetMeasMultcutName(icut)), fValues);
               if (fValues[AliReducedVarManager::kHighMultSPDTriggered] ) fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s_HighMultSPD", GetMeasMultcutName(icut)), fValues);
            }
         }

         fFilteredEvent->SetNGlobalTracks(fValues[AliReducedVarManager::kNGlobalTracks + icut], icut);
         // Fill mult regions relative to jpsi
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward + icut], 0, true, icut);
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse + icut], 1, true, icut);
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway + icut], 2, true, icut);
         // Fill mult regions relative to leading pt
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut], 0, false, icut);
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut], 1, false, icut);
         fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut], 2, false, icut);

         fFilteredEvent->SetLeadingParticle(fValues[AliReducedVarManager::kPtLeading + icut], fValues[AliReducedVarManager::kPhiLeading + icut], fValues[AliReducedVarManager::kEtaLeading + icut], icut);
      }
   }

}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::WriteFilteredPairs() {
   //
   // select and add filtered pair candidates to the filtered event
   //
   // loop over the pair list in the unfiltered event and evaluate all the pair cuts
   AliReducedPairInfo* pair = 0x0;
   TClonesArray* pairList = fEvent->GetPairs();
   if(!pairList) return;
   TIter nextPair(pairList);
   for(Int_t ip=0; ip<fEvent->NPairs(); ++ip) {
      pair = (AliReducedPairInfo*)nextPair();
      AliReducedVarManager::FillPairInfo(pair, fValues);
      fHistosManager->FillHistClass("Pair_BeforeCuts", fValues);
      for(UShort_t iflag=0; iflag<32; ++iflag) {
         AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
         fHistosManager->FillHistClass("PairQualityFlags_BeforeCuts", fValues);
      }
      
      if(IsPairSelected(pair, fValues)) {
         for(Int_t icut=0; icut<fPairCuts.GetEntries(); ++icut) {
            if(pair->TestFlag(icut)) {
               fHistosManager->FillHistClass(Form("Pair_%s", fPairCuts.At(icut)->GetName()), fValues);
               for(UShort_t iflag=0; iflag<32; ++iflag) {
                  AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
                  fHistosManager->FillHistClass(Form("PairQualityFlags_%s", fPairCuts.At(icut)->GetName()), fValues);
               }
            }
         }
         TClonesArray& pairs = *(fFilteredEvent->fCandidates);
         AliReducedPairInfo* filteredPair=NULL;
         filteredPair=new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo(*pair);
         fFilteredEvent->fNV0candidates[1] += 1;
      }
   }  // end loop over pairs
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::WriteFilteredTracks(Int_t array /*=1*/) {
   //
   // select and add filtered tracks to the filtered event
   //
   // loop over the track list and evaluate all the track cuts
   AliReducedBaseTrack* track = 0x0;
   TClonesArray* trackList = (array==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if(!trackList) return;
   if (fWriteFilteredTracksCandidatesOnly && (fFilteredEvent->NV0Candidates() == 0)) return;
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedBaseTrack*)nextTrack();
      AliReducedVarManager::FillTrackInfo(track, fValues);
      AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
      fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
      
      Bool_t writeTrack = IsTrackSelected(track, fValues);
      writeTrack |= (fWriteFilteredPairs && TrackIsCandidateLeg(track));

      if (fWriteFilteredTracksCandidatesOnly) writeTrack = TrackIsCandidateLeg(track) && !(track->IsMCTruth());
      
      if(writeTrack) {
         for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
            if(track->TestFlag(icut)) 
               fHistosManager->FillHistClass(Form("Track_%s", fTrackCuts.At(icut)->GetName()), fValues);
         }
         TClonesArray& tracks = (array==1 ? *(fFilteredEvent->fTracks) : *(fFilteredEvent->fTracks2));
      
         
         if(track->IsA() == AliReducedBaseTrack::Class()) {
            AliReducedBaseTrack* filteredParticle=NULL;
            filteredParticle=new(tracks[tracks.GetEntries()]) AliReducedBaseTrack(*track);
         }
         if(track->IsA() == AliReducedTrackInfo::Class()) {
            AliReducedTrackInfo* filteredParticle=NULL;
            AliReducedTrackInfo* tempTrack = dynamic_cast<AliReducedTrackInfo*>(track);
            //filteredParticle = (AliReducedTrackInfo*) tracks.ConstructedAt(tracks.GetEntries());
            new(tracks[tracks.GetEntries()]) AliReducedTrackInfo(*tempTrack);
            //filteredParticle = new(tracks[tracks.GetEntries()]) AliReducedTrackInfo(*tempTrack);
         }
         fFilteredEvent->fNtracks[1] += 1;
      }
   }  // end loop over tracks
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::BuildCandidatePairs() {
   //
   // Build candidate pairs and add them to the filtered event
   //
   // clear the track arrays
   fLeg1Tracks.Clear("C"); fLeg2Tracks.Clear("C"); 
   fLeg1PrefilteredTracks.Clear("C"); fLeg2PrefilteredTracks.Clear("C");
   RunCandidateLegsSelection(1);
   RunCandidateLegsSelection(2);
   
   if(fRunCandidatePrefilter) {
      RunCandidateLegsPrefilter(1);
      RunCandidateLegsPrefilter(2);
   }
   
   // Feed the selected tracks to the event mixing handler 
   if(fOptionRunMixing) {
      fMixingHandler->FillEvent(&fLeg1Tracks, &fLeg2Tracks, fValues, AliReducedPairInfo::kJpsiToEE);
      if(fOptionRunMixingMult){ // Mixed event in multiplicity bins
        AliMixingHandler* handler;
        TIter nextHandler(&fMixingHandlerMult); 
        for (int i=0; i<fMixingHandlerMult.GetEntries(); i++) {
            handler = (AliMixingHandler*) nextHandler();
            if (fValues[AliReducedVarManager::kNGlobalTracks]>=fMultBinsMixing[i] && fValues[AliReducedVarManager::kNGlobalTracks]<fMultBinsMixing[i+1]) {
               handler->FillEvent(&fLeg1Tracks, &fLeg2Tracks, fValues, AliReducedPairInfo::kJpsiToEE);
            }
        }
     }
   }
   
   if(fLeg1Tracks.GetEntries() + fLeg2Tracks.GetEntries() > 1) 
      RunSameEventPairing();
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunCandidateLegsSelection(Int_t arrayOption /*=1*/) {
   //
   // select leg candidates and prefilter tracks
   // NOTE: In the case of symmetric decay channels, the candidates are separated by charge in fLeg1Tracks and fLeg2Tracks
   //            For asymmetric decays, the candidates for each of the LEG1 and LEG2 are not a priori separated by charge.
   //                 It is the responsability of the analyzer to setup the track cuts such that these are separated
   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   Bool_t mcDecision = kTRUE;
   
   // loop over the track list and evaluate all the track cuts
   AliReducedBaseTrack* track = 0x0;
   TClonesArray* trackList = (arrayOption==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if (!trackList) return;
   TIter nextTrack(trackList);

   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if (track->IsMCKineParticle()) continue;
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      AliReducedVarManager::FillTrackInfo(track, fValues);
      AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
      // NOTE: mcDecision works with just one MC truth cut. It is not implemented properly for asymmetric decay channels
      if(fOptionRunOverMC && fLegCandidatesMCcuts) mcDecision = CheckReconstructedLegMCTruth(track);
      if (arrayOption == 1) { // candidates tracks are only in array 1, but prefilter tracks may be in array 2
         if(isAsymmetricDecayChannel) {
            if(IsCandidateLegSelected(track, fValues, 1)) {
               fLeg1Tracks.Add(track);
               FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, fValues, 1, isAsymmetricDecayChannel);
            }
            if(IsCandidateLegSelected(track, fValues, 2)) {
               fLeg2Tracks.Add(track);
               FillCandidateLegHistograms("Track_LEG2_BeforePrefilter", track, fValues, 2, isAsymmetricDecayChannel);
            }
         }
         else {
            // mcDecision is by default true, it can be false only if running on MC and the track fails the test 
            if(IsCandidateLegSelected(track, fValues,1) && mcDecision) {
               if(track->Charge()>0) {
                  fLeg1Tracks.Add(track);
                  FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, fValues, 1, isAsymmetricDecayChannel);
               }
               if(track->Charge()<0) {
                  fLeg2Tracks.Add(track);
                  FillCandidateLegHistograms("Track_LEG2_BeforePrefilter", track, fValues, 1, isAsymmetricDecayChannel);
               }
            } 
         }
      }
      
      if(!fRunCandidatePrefilter) continue;
      
      if(isAsymmetricDecayChannel) {
         if(IsCandidateLegPrefilterSelected(track, fValues, 1)) {
            fLeg1PrefilteredTracks.Add(track);
            fHistosManager->FillHistClass("Track_LEG1_PrefilterTrack", fValues);
         }
         if(IsCandidateLegPrefilterSelected(track, fValues, 2)) {
            fLeg2PrefilteredTracks.Add(track);
            fHistosManager->FillHistClass("Track_LEG2_PrefilterTrack", fValues);
         }
      }
      else {
         if(IsCandidateLegPrefilterSelected(track, fValues)) {
            if(track->Charge()>0) {
               fLeg1PrefilteredTracks.Add(track);
               fHistosManager->FillHistClass("Track_LEG1_PrefilterTrack", fValues);
            }
            if(track->Charge()<0) {
               fLeg2PrefilteredTracks.Add(track);
               fHistosManager->FillHistClass("Track_LEG2_PrefilterTrack", fValues);
            }
         }
      }
   }   // end loop over tracks

}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunCandidateLegsPrefilter(Int_t leg) {
   //
   // Run the prefilter selection
   // At this point it is assumed that the track lists are filled
   //
   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   
   // initialize iterators
   TIter iterLeg((leg==1 ? &fLeg1Tracks : &fLeg2Tracks));
   TIter iterPrefLeg1(&fLeg1PrefilteredTracks);
   TIter iterPrefLeg2(&fLeg2PrefilteredTracks);
  
   // Pair the LEG candidates with the prefilter selected tracks
   AliReducedBaseTrack* track=0;
   AliReducedBaseTrack* prefTrack=0;
   for(Int_t it = 0; it<(leg==1 ? fLeg1Tracks.GetEntries() : fLeg2Tracks.GetEntries()); ++it) {
      track = (AliReducedBaseTrack*)iterLeg();
      
      if((leg==1 && (isAsymmetricDecayChannel || (!isAsymmetricDecayChannel && fRunCandidatePrefilterOnSameCharge))) ||
         (leg==2 && !isAsymmetricDecayChannel)) 
      {
         iterPrefLeg1.Reset();
         for(Int_t it2 = 0; it2<fLeg1PrefilteredTracks.GetEntries(); ++it2) {
            prefTrack = (AliReducedBaseTrack*)iterPrefLeg1();
            
            if(track->TrackId()==prefTrack->TrackId()) continue;      // avoid self-pairing
            AliReducedVarManager::FillPairInfo(track, prefTrack, fCandidateType, fValues);
            if(!IsCandidateLegPairPrefilterSelected(fValues, 1)) {
               track->ResetFlags();
               break;
            }
         }  // end loop over tracks
      }  // end if
      
      if((leg==1 && !isAsymmetricDecayChannel) ||
         (leg==2 && (isAsymmetricDecayChannel || (!isAsymmetricDecayChannel && fRunCandidatePrefilterOnSameCharge))))
      {
         iterPrefLeg2.Reset();
         for(Int_t it2 = 0; it2<fLeg2PrefilteredTracks.GetEntries(); ++it2) {
            prefTrack = (AliReducedBaseTrack*)iterPrefLeg2();
            
            if(track->TrackId()==prefTrack->TrackId()) continue;       // avoid self-pairing
            AliReducedVarManager::FillPairInfo(track, prefTrack, fCandidateType, fValues);
            if(!IsCandidateLegPairPrefilterSelected(fValues, 2)) {
               track->ResetFlags();
               break;
            }
         }  // end loop over tracks
      }  // end if
   }  // end loop over tracks
   
   // remove tracks
   iterLeg.Reset();
   for(Int_t it = (leg==1 ? fLeg1Tracks.GetEntries() : fLeg2Tracks.GetEntries()) - 1; it >= 0; --it) {
      track = (AliReducedBaseTrack*)iterLeg();
      if(!track->GetFlags()) {
         if(leg==1) fLeg1Tracks.Remove(track);
         if(leg==2) fLeg2Tracks.Remove(track);
      }
   }
   
   // fill histograms after the prefilter
   iterLeg.Reset();
   for(Int_t it = 0; it<(leg==1 ? fLeg1Tracks.GetEntries() : fLeg2Tracks.GetEntries()); ++it) {
      track = (AliReducedBaseTrack*)iterLeg();
      AliReducedVarManager::FillTrackInfo(track, fValues);
      AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
      FillCandidateLegHistograms(Form("Track_LEG%d_AfterPrefilter", leg), track, fValues, (leg==2 && isAsymmetricDecayChannel ? 2 : 1), isAsymmetricDecayChannel);
   }
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunSameEventPairing() {
   //
   // Run the same event pairing
   //   

   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   TIter iterLeg1(&fLeg1Tracks);
   TIter iterLeg2(&fLeg2Tracks);
   
   AliReducedTrackInfo* leg1Track=0;
   AliReducedTrackInfo* leg2Track=0;
   AliReducedTrackInfo* leg1Track_2=0;

   for(Int_t it1=0; it1<fLeg1Tracks.GetEntries(); ++it1) {
      leg1Track = (AliReducedTrackInfo*)iterLeg1();
      
      iterLeg2.Reset();
      for(Int_t it2=0; it2<fLeg2Tracks.GetEntries(); ++it2) {
         leg2Track = (AliReducedTrackInfo*)iterLeg2();
         
         // verify that the two current tracks have at least 1 common bit
         ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg2Track, isAsymmetricDecayChannel);
         if(!compatibilityMask) continue;
         AliReducedVarManager::FillPairInfo(leg1Track, leg2Track, fCandidateType, fValues);

         if (fOptionRunOverMC){
            bool isJpsi = abs(leg1Track->MCPdg(0))==11 && abs(leg2Track->MCPdg(0))==11 && leg1Track->MCLabel(1)==leg2Track->MCLabel(1) && leg1Track->MCPdg(1)==443;
            bool isJpsiFromB = isJpsi && ((abs(leg1Track->MCPdg(2))>500 && abs(leg1Track->MCPdg(2))<599) || (abs(leg1Track->MCPdg(2))>5000 && abs(leg1Track->MCPdg(2))<5999));
            fValues[AliReducedVarManager::kPairMCMap] = isJpsi+2*isJpsiFromB;
         }
         if(!IsCandidatePairSelected(fValues)) continue;
         if(fOptionRunOverMC && fLegCandidatesMCcuts && (!CheckReconstructedLegMCTruth(leg1Track, leg2Track))) continue; 

         TClonesArray& pairs = *(fFilteredEvent->fCandidates);
         AliReducedPairInfo* candidatePair = new (pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
         candidatePair->SetLegIds(leg1Track->TrackId(), leg2Track->TrackId());
         candidatePair->SetFlags(compatibilityMask);
         FillCandidatePairHistograms("Pair_Candidate12", candidatePair, fValues, isAsymmetricDecayChannel);
         SetupPair(candidatePair, fValues);
         fFilteredEvent->fNV0candidates[1] += 1;

         // Fill legs histogram for QA
         AliReducedVarManager::FillTrackInfo(leg1Track, fValues);
         FillCandidateLegHistograms("Track_LEG1_AfterPairCuts", leg1Track, fValues, 1, isAsymmetricDecayChannel);
         AliReducedVarManager::FillTrackInfo(leg2Track, fValues);
         FillCandidateLegHistograms("Track_LEG2_AfterPairCuts", leg2Track, fValues, 1, isAsymmetricDecayChannel);

      }  // end loop over LEG2 tracks
      
      if(fBuildCandidateLikePairs) {
         for(Int_t it1_2=it1+1; it1_2<fLeg1Tracks.GetEntries(); ++it1_2) {
            leg1Track_2 = (AliReducedTrackInfo*)fLeg1Tracks.At(it1_2);

            // verify that the two current tracks have at least 1 common bit
            ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg1Track_2, isAsymmetricDecayChannel);
            if(!compatibilityMask) continue;
            AliReducedVarManager::FillPairInfo(leg1Track, leg1Track_2, fCandidateType, fValues);
            fValues[AliReducedVarManager::kPairMCMap] = 0;

            if(!IsCandidatePairSelected(fValues)) continue;
            TClonesArray& pairs = *(fFilteredEvent->fCandidates);
            AliReducedPairInfo* candidatePair=new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
            candidatePair->SetLegIds(leg1Track->TrackId(), leg1Track_2->TrackId());
            candidatePair->SetFlags(compatibilityMask);
            FillCandidatePairHistograms("Pair_Candidate11", candidatePair, fValues, isAsymmetricDecayChannel);
            SetupPair(candidatePair, fValues);
            fFilteredEvent->fNV0candidates[1] += 1;            
         }  // end loop over LEG1 tracks
      }  // end if(fBuildCandidateLikePairs)
   }  // end loop over LEG1 tracks
   
   if(fBuildCandidateLikePairs) {
      AliReducedBaseTrack* leg2Track_2 = 0;
      iterLeg2.Reset();
      for(Int_t it2=0; it2<fLeg2Tracks.GetEntries(); ++it2) {
         leg2Track = (AliReducedTrackInfo*)iterLeg2();
         
         for(Int_t it2_2=it2+1; it2_2<fLeg2Tracks.GetEntries(); ++it2_2) {
            leg2Track_2 = (AliReducedTrackInfo*)fLeg2Tracks.At(it2_2);
            
            // verify that the two current tracks have at least 1 common bit
            ULong_t compatibilityMask = CheckTrackCompatibility(leg2Track, leg2Track_2, isAsymmetricDecayChannel);
            if(!compatibilityMask) continue;
            AliReducedVarManager::FillPairInfo(leg2Track, leg2Track_2, fCandidateType, fValues);
            fValues[AliReducedVarManager::kPairMCMap] = 0;
            if(!IsCandidatePairSelected(fValues)) continue;
            TClonesArray& pairs = *(fFilteredEvent->fCandidates);
            AliReducedPairInfo* candidatePair=new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
            candidatePair->SetLegIds(leg2Track->TrackId(), leg2Track_2->TrackId());
            candidatePair->SetFlags(compatibilityMask);
            FillCandidatePairHistograms("Pair_Candidate22", candidatePair, fValues, isAsymmetricDecayChannel);
            SetupPair(candidatePair, fValues);
            fFilteredEvent->fNV0candidates[1] += 1;           
         }  // end loop over negative tracks
      }  // end loop over negative tracks
   }
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillMultiplicity(Bool_t regions /*= kFALSE*/, bool wrtPhiRef /*= kFALSE*/, float phiRef /*= -999*/){
   //
   // Fill multiplicity values (global if regions = false; in the regions if regions = true)
   // You can also give first an angle wrt to the regions should be calculated
   //

   // Fill global tracks (both signal and MC and Jpsi)
   if (!regions) {
      for (int icut = 0; icut<GetNMeasMultCuts(); icut++)  fValues[AliReducedVarManager::kNGlobalTracks + icut] = 0.;
      fValues[AliReducedVarManager::kMCNch09] = 0.;
      fValues[AliReducedVarManager::kMCNJpsi] = 0.;
      fValues[AliReducedVarManager::kMCNch] = 0.;
   }
   float phi = phiRef;
   if (regions) {
      if (!wrtPhiRef) {
         //Defining the phi reference 
         if (fDefaultRandomPhi) phi = TMath::TwoPi() * gRandom->Rndm();
         else phi = fValues[AliReducedVarManager::kPhiLeading];
         
         if (!fRegionsToMCTruth) {
            TClonesArray* pairs = fFilteredEvent->fCandidates;
            TIter nextPair(pairs);
            AliReducedPairInfo* jpsiPair;
            float maxProbJpsi = 0.;
            
            //Looping on jpsi candidates
            for (int i=0; i<pairs->GetEntries(); i++) {
               jpsiPair = (AliReducedPairInfo*) nextPair();
               if(jpsiPair->PairType()!=1) continue; // reject like-sign pairs
               if(fJpsiMassDist) { // Take the jpsi with the maximal probability to be a real Jpsi (given its mass, if possible)
                  float probJpsi = fJpsiMassDist->Eval(jpsiPair->Mass(0)); 
                  if(probJpsi >= maxProbJpsi) {
                     maxProbJpsi = probJpsi;
                     phi = jpsiPair->Phi();
                  }
               }
               else phi = jpsiPair->Phi();
            }
            if ((fValues[AliReducedVarManager::kMCNJpsi] == 1) && (pairs->GetEntries() > 0.)) {
               fValues[AliReducedVarManager::kPhiJpsiMCTruth] = abs(fValues[AliReducedVarManager::kPhiJpsiMCTruth] - phi);
               fHistosManager->FillHistClass("DeltaPhi_JpsiTruth_JpsiCandidate", fValues);
            }
         }
         if (fRegionsToMCTruth && (fValues[AliReducedVarManager::kMCNJpsi] == 1)) {
            phi = fValues[AliReducedVarManager::kPhiJpsiMCTruth];
         }
      }

      for (int icut = 0; icut < GetNMeasMultCuts(); icut++) {
         fValues[AliReducedVarManager::kNGlobalTracksToward + icut] = 0.;
         fValues[AliReducedVarManager::kNGlobalTracksTransverse + icut] = 0.;
         fValues[AliReducedVarManager::kNGlobalTracksAway + icut] = 0.;
         fValues[AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] = 0.;
         fValues[AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] = 0.;
         fValues[AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] = 0.;
      }
      fValues[AliReducedVarManager::kMCNch09Toward] = 0.;
      fValues[AliReducedVarManager::kMCNch09Transverse] = 0.;
      fValues[AliReducedVarManager::kMCNch09Away] = 0.;
      fValues[AliReducedVarManager::kMCNch09Toward + 1] = 0.;
      fValues[AliReducedVarManager::kMCNch09Transverse + 1] = 0.;
      fValues[AliReducedVarManager::kMCNch09Away + 1] = 0.;
   }


   std::vector<int> weightsTrue;
   std::vector<int> weightsMeas;
   // Run on both arrays 
   for(int iArray = 1; iArray <= 2; iArray++) {
      AliReducedTrackInfo* track;
      TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
      TIter nextTrack(tracklist);

      for(Int_t it = 0; it < tracklist->GetEntries(); ++it) {
         track = (AliReducedTrackInfo*)nextTrack();  
         for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
         AliReducedVarManager::FillTrackInfo(track, fValues);

         if (!regions && track->IsMCTruth() && track->Charge() && abs(track->Eta()) < 1.) fValues[AliReducedVarManager::kMCNch] ++;

         if(!regions && (iArray == 1) && track->IsMCTruth()) {
            if ((track->MCPdg(0) == 443) && CheckMotherMCTruth(track)) { //also fill MC truth jpsi 
               // find the jpsi daughters (to check if dielectron)
               if (fMCTruthJpsi2eeOnly) {
                  Int_t daughter1Label = 0; Int_t daughter2Label = 0;
                  FindJpsiTruthLegs(track, daughter1Label, daughter2Label);
                  if ((daughter1Label == 0) || (daughter2Label == 0)) continue; // Not dielectron decay   
               }   
               fValues[AliReducedVarManager::kMCNJpsi] += 1;
               fValues[AliReducedVarManager::kPhiJpsiMCTruth] = track->Phi();
            }
         }

         if (!track->IsMCTruth() && !regions) IsTrackMeasuredMultSelected(track, fValues); // set mult flags, if regions they have already been set
         if ( track->IsMCTruth() && !regions) IsTrackTrueMultSelected(track, fValues);

         if (!track->GetMultFlags()) continue;

         if(!(track->IsMCTruth())) {
            int particleWeight = 1;
            if (fOptionRunOverMC && fReweightParticleComposition) particleWeight = GetParticleWeight(track);

           /* if (fOptionRunOverMC && fSharePCCWeightsBetweenEvents && !regions) { // The additional weights are given to other events
               if (particleWeight != 1 && track->TestFlag(5)) { //TODO: The cut number is hardcoded, it should not be
                  weightsMeas.push_back(particleWeight - 1);
               }
               particleWeight = 1;
            }*/

            // Counting measured tracks
            for (int icut = 0; icut<GetNMeasMultCuts(); icut++) {
               if (track->TestMultFlag(icut)) {
                  if (!regions && (track->Pt() > fValues[AliReducedVarManager::kPtLeading+icut]) && !fRegionsToMCTruth) { 
                     //Look for leading particle
                     fValues[AliReducedVarManager::kPtLeading+icut] = track->Pt();
                     fValues[AliReducedVarManager::kPhiLeading+icut] = track->Phi();
                     fValues[AliReducedVarManager::kEtaLeading+icut] = track->Eta();
                  }
                  if(!regions) 
                     fValues[AliReducedVarManager::kNGlobalTracks + icut] += particleWeight;
                  else { 
                     // Region to Jpsi
                     float phiIcut = fValues[AliReducedVarManager::kPhiLeading+icut];
                     float delta_phi = abs(track->Phi() - phi);
                     if (phi == fValues[AliReducedVarManager::kPhiLeading]) delta_phi = abs(track->Phi() - phiIcut);
                     if (delta_phi < M_PI/3. || delta_phi > 5*M_PI/3.) 
                        fValues[AliReducedVarManager::kNGlobalTracksToward+icut] += particleWeight;
                     if ((delta_phi > M_PI/3. && delta_phi < 2*M_PI/3.) || (delta_phi > 4*M_PI/3. && delta_phi < 5*M_PI/3.)) 
                        fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut] += particleWeight;
                     if (delta_phi > 2*M_PI/3. && delta_phi < 4*M_PI/3.) 
                        fValues[AliReducedVarManager::kNGlobalTracksAway+icut] += particleWeight;

                     // Region to leading
                     float delta_phi_leading = abs(track->Phi() - phiIcut);
                     if (delta_phi_leading < M_PI/3. || delta_phi_leading > 5*M_PI/3.) 
                        fValues[AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] += particleWeight;
                     if ((delta_phi_leading > M_PI/3. && delta_phi_leading < 2*M_PI/3.) || (delta_phi_leading > 4*M_PI/3. && delta_phi_leading < 5*M_PI/3.)) 
                        fValues[AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] += particleWeight;
                     if (delta_phi_leading > 2*M_PI/3. && delta_phi_leading < 4*M_PI/3.)
                        fValues[AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut] += particleWeight;
                  }
               }
            } 
         }

         if(track->IsMCTruth()) { 
            int particleWeight = 1;
            if (fOptionRunOverMC && fReweightParticleComposition) particleWeight = GetParticleWeight(track);

            /*if (fOptionRunOverMC && fSharePCCWeightsBetweenEvents && !regions) { // The additional weights are given to other events
               if (particleWeight != 1) {
                  weightsTrue.push_back(particleWeight - 1);
               }
               particleWeight = 1;
            }*/

            // Counting true tracks
            if (!regions && (track->Pt() > fValues[AliReducedVarManager::kPtLeading]) && fRegionsToMCTruth) { 
               // Look for leading particle
               for (int icut = 0; icut<GetNMeasMultCuts(); icut++) {
                  fValues[AliReducedVarManager::kPtLeading + icut] = track->Pt();
                  fValues[AliReducedVarManager::kPhiLeading + icut] = track->Phi();
                  fValues[AliReducedVarManager::kEtaLeading + icut] = track->Eta();
               }
            }
            if(!regions) {
               fValues[AliReducedVarManager::kMCNch09] += particleWeight;
            }
            else {
               //region to Jpsi
               float delta_phi = abs(track->Phi() - phi);
               if (delta_phi < M_PI/3. || delta_phi > 5*M_PI/3.) 
                  fValues[AliReducedVarManager::kMCNch09Toward] += particleWeight;
               if ((delta_phi > M_PI/3. && delta_phi < 2*M_PI/3.) || (delta_phi > 4*M_PI/3. && delta_phi < 5*M_PI/3.)) 
                  fValues[AliReducedVarManager::kMCNch09Transverse] += particleWeight;
               if (delta_phi > 2*M_PI/3. && delta_phi < 4*M_PI/3.) 
                  fValues[AliReducedVarManager::kMCNch09Away] += particleWeight;
               //region to leading
               float delta_phi_leading = abs(track->Phi() - fValues[AliReducedVarManager::kPhiLeading]);
               if (delta_phi_leading < M_PI/3. || delta_phi_leading > 5*M_PI/3.) 
                  fValues[AliReducedVarManager::kMCNch09Toward + 1] += particleWeight;
               if ((delta_phi_leading > M_PI/3. && delta_phi_leading < 2*M_PI/3.) || (delta_phi_leading > 4*M_PI/3. && delta_phi_leading < 5*M_PI/3.)) 
                  fValues[AliReducedVarManager::kMCNch09Transverse + 1] += particleWeight;
               if (delta_phi_leading > 2*M_PI/3. && delta_phi_leading < 4*M_PI/3.) 
                  fValues[AliReducedVarManager::kMCNch09Away + 1] += particleWeight;
            }
         }

      }  
   }

/*
   if (fOptionRunOverMC && fSharePCCWeightsBetweenEvents && !regions) {
      float nTrueNoWeights = fValues[AliReducedVarManager::kMCNch09]; // Take previous weights from events with similar properties
      std::vector<int> prevWeightsTrue = fWeightsTrue[nTrueNoWeights];
      std::vector<int> prevWeightsMeas = fWeightsMeas[nTrueNoWeights];

      for (auto iterTrue = prevWeightsTrue.begin(); iterTrue != prevWeightsTrue.end(); iterTrue++) {
         int weight = *iterTrue;
         if (weight < 0) fValues[AliReducedVarManager::kMCNch09] -= 1;
         if (weight > 0) {
            fValues[AliReducedVarManager::kMCNch09] += 1;
            if (weight > 1) weightsTrue.push_back(weight - 1);
         }
      }

      for (auto iterMeas = prevWeightsMeas.begin(); iterMeas != prevWeightsMeas.end(); iterMeas++) {
         int weight = *iterMeas;
         if (weight < 0) fValues[AliReducedVarManager::kNGlobalTracks] -= 1; // For now only one cut
         if (weight > 0) {
            fValues[AliReducedVarManager::kNGlobalTracks] += 1;
            if (weight > 1) weightsMeas.push_back(weight - 1);
         }
      }

      fWeightsTrue[nTrueNoWeights] = weightsTrue;
      fWeightsMeas[nTrueNoWeights] = weightsMeas;
   }*/




   //some quantities necessary to compute efficiency/contamination
   fValues[AliReducedVarManager::kMCNch09+2] = fValues[AliReducedVarManager::kMCNch09];
   fValues[AliReducedVarManager::kMCNch09+1] = fValues[AliReducedVarManager::kMCNch09];

   float VtxZMC = fValues[AliReducedVarManager::kVtxZMC];
   bool isTriggered = kTRUE;
   bool isMCAccepted = (fValues[AliReducedVarManager::kMCNch] > 0) && (abs(VtxZMC) <= 10.);

   if(!isTriggered){
      fValues[AliReducedVarManager::kMCNch09+1] = -999.;
      fValues[AliReducedVarManager::kMCNch09+2] = -999.;
   }

   if(!isMCAccepted){
      fValues[AliReducedVarManager::kMCNch09] = -999.;
      fValues[AliReducedVarManager::kMCNch09+1] = -999.;
      fValues[AliReducedVarManager::kMCNch09+2] = -999.;
      if (regions) {
         fValues[AliReducedVarManager::kMCNch09Toward] = -999.;
         fValues[AliReducedVarManager::kMCNch09Away] = -999.;
         fValues[AliReducedVarManager::kMCNch09Transverse] = -999.;
      }
   }
   
}

//____________________________________________________________________________
Int_t AliReducedAnalysisFilterTrees::GetParticleWeight(AliReducedTrackInfo* track) {
   //
   // Get particle weight for multiplicity calculation in MC
   // The weight depends on pt, species, (multiplicity class)
   // It is used while reweighting particle composition in MC in order to match it to data
   // Particle type: pion = 0, proton = 1, kaon = 2, sigma- = 3, sigma+ = 4, rest = 5, weight_lambda = weight_xi = weight_sigma+
   // see $ALICE_PHYSICS/PWG/Tools/AliMCSpectraWeights.h and https://alice-notes.web.cern.ch/node/1311
   //

   float ptMC = track->Pt();
   int type = 5;
   // For now flags are hard-coded
   if (track->TestMCFlag(7)) type = 0; // pion primary
   if (track->TestMCFlag(8)) type = 2; // kaon primary
   if (track->TestMCFlag(9)) type = 1; // proton primary
   if (track->TestMCFlag(10)) type = 4; // sigma+ primary
   if (track->TestMCFlag(11)) type = 3; // sigma- primary

   /*if (track->TestMCFlag(14)) type = 0; // pion primary
   if (track->TestMCFlag(15)) type = 2; // kaon primary
   if (track->TestMCFlag(16)) type = 1; // proton primary
   if (track->TestMCFlag(17)) type = 4; // sigma+ primary
   if (track->TestMCFlag(18)) type = 3; // sigma- primary*/

   if (type == 5 && !track->IsMCTruth() && !track->TestMCFlag(6/*13*/)) {
      // Check if particle is a detected secondary from strange decay
      // In this case, the particle was put in array 1
      int pdgmother = track->MCPdg(1);
      if (abs(pdgmother) == 3222) type = 4; 
      if (abs(pdgmother) == 3112) type = 3;
      if (abs(pdgmother) == 310) type = 2;
      if (abs(pdgmother) == 3122) type = 4;

      if (type != 5) {
         // The mother is strange primary - we look for its pt in array 1
         int labelmother = track->MCLabel(1);
         // loop over the first track array
         TClonesArray* trackList = fEvent->GetTracks();
         TIter nextTrack(trackList);
         AliReducedTrackInfo* mother;
         for(Int_t it = 0; it < trackList->GetEntries(); ++it) {
            mother = (AliReducedTrackInfo*)nextTrack();
            if (mother->IsMCTruth() && (mother->MCLabel(0) == labelmother)){
               // Get the pt of the mother
               if (mother->TestMCFlag(6/*13*/)) // mother is primary
                  ptMC = mother->PtMC(); 
               else type = 5;
               break;
            }
            if (it == trackList->GetEntries() - 1) type = 5; // We did not find the mother: do nothing
         }
      }
   }
   if (type == 5) return 1;
      
   // For multiplicity class (only rough value needed, because not so big differences in pp)
   float mult = fValues[AliReducedVarManager::kMultEstimatorPercentileV0M];
   if (mult < 0. || mult > 100.) mult = 50.; // multiplicity estimator percentile was probably not kept

   int iptMC = fMCParticleCompositionWeights->GetXaxis()->FindBin(ptMC);
   int imult = fMCParticleCompositionWeights->GetYaxis()->FindBin(mult);
   int itype = fMCParticleCompositionWeights->GetZaxis()->FindBin(type);
   
   // In case out of the range
   if (ptMC < fMCParticleCompositionWeights->GetXaxis()->GetBinLowEdge(1)) 
      iptMC = 1;
   if (ptMC > fMCParticleCompositionWeights->GetXaxis()->GetBinUpEdge(fMCParticleCompositionWeights->GetNbinsX())) 
      iptMC = fMCParticleCompositionWeights->GetNbinsX();
   if (mult > fMCParticleCompositionWeights->GetYaxis()->GetBinUpEdge(fMCParticleCompositionWeights->GetNbinsY())) 
      imult = fMCParticleCompositionWeights->GetNbinsY();

   float weight = fMCParticleCompositionWeights->GetBinContent(iptMC, imult, itype);
   int intWeight = (int) weight;
   float rest = weight - intWeight;
   if (gRandom->Rndm() <= rest) intWeight++;

   return intWeight;
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
   //
   // apply event cuts
   //
   if(fEventCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical and between all cuts in the list
   for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
      if(values) { if(!cut->IsSelected(event, values)) return kFALSE; }
      else { if(!cut->IsSelected(event)) return kFALSE; }
   }
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
   //
   // apply track cuts
   //
   if(fTrackCuts.GetEntries()==0) return kTRUE;
   track->ResetFlags();
   
   // loop over all the cuts and toggle a filter bit if the track passes the corresponding cut
   for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
      if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
      else { if(cut->IsSelected(track)) track->SetFlag(i); }
   }
   return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::TrackIsCandidateLeg(AliReducedBaseTrack* track) {
   //
   // if the track is a leg belonging to at least one of the candidate pairs written, then return kTRUE
   //
   AliReducedPairInfo* pair = 0x0;
   TClonesArray* pairList = fFilteredEvent->GetPairs();
   TIter nextPair(pairList);
   for(Int_t ip=0; ip<fFilteredEvent->NPairs(); ++ip) {
      pair = (AliReducedPairInfo*)nextPair();
      if(track->TrackId()==pair->LegId(0)) return kTRUE;
      if(track->TrackId()==pair->LegId(1)) return kTRUE;
   }
   return kFALSE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsPairSelected(AliReducedPairInfo* pair, Float_t* values /*=0x0*/) {
   //
   // apply pair cuts
   // NOTE: Multiple cut sets are supported. The decisions are encoded in the fFlags inherited from AliReducedBaseTrack
   if(fPairCuts.GetEntries()==0) return kTRUE;
   pair->ResetFlags();
   
   // loop over all the cuts and toggle a filter bit if the pair passes the corresponding cut
   for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
      if(values) { if(cut->IsSelected(pair, values)) pair->SetFlag(i); }
      else { if(cut->IsSelected(pair)) pair->SetFlag(i); }
   }
   return (pair->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackMeasuredMultSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/) {
   //
   // apply cuts for determining measured multiplicity
   //

   if(fMeasuredMultTrackCuts.GetEntries()==0) return kTRUE;
   if (track->IsMCTruth()) return kFALSE;

   track->ResetMultFlags();  

   for(Int_t i=0; i<fMeasuredMultTrackCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*) fMeasuredMultTrackCuts.At(i) ;
      if(values) { 
         if(cut->IsSelected(track, values))
            track->SetMultFlag(i);
      }
      else { 
         if(cut->IsSelected(track))
            track->SetMultFlag(i);
      }
   }
   return (track->GetMultFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackTrueMultSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/) {
   //
   // apply cuts for determining MC true multiplicity
   //
   if(fTrueMultTrackCuts.GetEntries()==0) return kTRUE;
   if (!track->IsMCTruth()) return kFALSE;

   track->ResetMultFlags();  

   for(Int_t i=0; i<fTrueMultTrackCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*) fTrueMultTrackCuts.At(i) ;
      if(values) { 
         if(cut->IsSelected(track, values)) 
            track->SetMultFlag(i); 
      }
      else { 
         if(cut->IsSelected(track)) 
            track->SetMultFlag(i); 
      }
   }
   return (track->GetMultFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/, Int_t whichLeg /*=1*/) {
   //
   // apply cuts to the candidate leg
   //
   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   
   if(fLeg1Cuts.GetEntries()==0) return kTRUE;
   
   // reset the flags for the track only if these are the cuts on LEG1
   // For LEG2, we use the same bit map as for LEG1, which it is assumed was already evaluated so the ResetFlags() should not be called
   // IMPORTANT: In the case of asymmetric decay channels, the LEG1 cuts have to be always evaluated before LEG2 cuts
   if(whichLeg==1) track->ResetFlags();
   
   for(Int_t i=0; i<fLeg1Cuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel ? fLeg2Cuts.At(i) : fLeg1Cuts.At(i));
      if(values) { 
         if(cut->IsSelected(track, values)) 
            track->SetFlag(whichLeg==2 && isAsymmetricDecayChannel ? 32+i : i); 
      }
      else { 
         if(cut->IsSelected(track)) 
            track->SetFlag(whichLeg==2 && isAsymmetricDecayChannel ? 32+i : i); 
      }
   }
   return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPrefilterSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/, Int_t whichLeg /*=1*/) {
   //
   // apply prefilter cuts on the candidate legs
   //   
   if(fLeg1PrefilterCuts.GetEntries()==0) return kTRUE;
   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   
   for(Int_t i=0; i<fLeg1PrefilterCuts.GetEntries(); ++i) {
      // if there are more cuts specified, we apply an AND on all of them
      // NOTE: The analysis task could also be configured such that there is a prefilter track cut corresponding to each track cut
      //              in which case one needs to make sure the number of prefilter cuts is the same as the number of track cuts
      AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel ? fLeg2PrefilterCuts.At(i) : fLeg1PrefilterCuts.At(i));
      if(values) { if(!cut->IsSelected(track, values)) return kFALSE; }
      else { if(!cut->IsSelected(track)) return kFALSE; }
   }
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidatePairSelected(Float_t* values) {
   //
   // apply pair cuts to the built candidates
   //
   if(fCandidatePairCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical and between all cuts in the list
   for(Int_t i=0; i<fCandidatePairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fCandidatePairCuts.At(i);
      if(!cut->IsSelected(values)) return kFALSE;
   }
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPairPrefilterSelected(Float_t* values, Int_t whichLeg /*=1*/) {
   //
   // apply the prefilter pair cuts
   //   
   if(fLeg1PairPrefilterCuts.GetEntries()==0) return kTRUE;
   Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
   
   // loop over all the cuts and make a logical OR between all cuts in the list
   for(Int_t i=0; i<fLeg1PairPrefilterCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel ? fLeg2PairPrefilterCuts.At(i) : fLeg1PairPrefilterCuts.At(i));
      if(cut->IsSelected(values)) return kTRUE;
   }
   return kFALSE;
}

//___________________________________________________________________________
ULong_t AliReducedAnalysisFilterTrees::CheckTrackCompatibility(AliReducedBaseTrack* leg1, AliReducedBaseTrack* leg2, Bool_t isAsymmetricDecayChannel) {
   //
   // check whether the 2 tracks fulfill at least one set of common cuts
   // NOTE: In the case of asymmetric decay channels, a track can fulfill simultaneously both the cuts of LEG1 and LEG2
   //              and since we work with just one instance of the track object in memory, the cut flags map is shared
   //              between the two LEG cut sets: [0-31) cuts for LEG1, [32-63) cuts for LEG2.
   //              We should pair just tracks fulfilling cuts from the same "doublet": e.g. 0 - 32, 1-33, 2-34, ...
   //             In the case of symmetric decays, there is only one list of cuts for which the tracks are evaluated
   //               so here there is no problem. The flags of the tracks forming a pair are evaluated by a simple bitwise AND
   if(!isAsymmetricDecayChannel) return (leg1->GetFlags() & leg2->GetFlags());
   
   ULong_t mask = 0;
   for(Int_t i=0;i<32;++i)
      if(leg1->TestFlag(i) && leg2->TestFlag(32+i)) 
         mask |= (ULong_t(1)<<i);
   return mask;
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::SetupPair(AliReducedPairInfo* pair, Float_t* values) {
   //
   // setup pair information
   //   
   
   pair->Pt(values[AliReducedVarManager::kPt]);
   pair->Phi(values[AliReducedVarManager::kPhi]);
   pair->Eta(values[AliReducedVarManager::kEta]);
   pair->Charge(0);
   pair->CandidateId(fCandidateType);
   pair->PairType(values[AliReducedVarManager::kPairType]);
   pair->PairTypeSPD(values[AliReducedVarManager::kPairTypeSPD]);
   pair->SetMass(values[AliReducedVarManager::kMass]);
   pair->SetLxy(values[AliReducedVarManager::kPairLxy]);
   pair->SetPseudoProper(values[AliReducedVarManager::kPseudoProperDecayTime]);
   pair->SetPointingAngle(values[AliReducedVarManager::kPairPointingAngle]);
   pair->SetChisquare(values[AliReducedVarManager::kPairChi2prNDOF]);
   pair->SetMCMap(values[AliReducedVarManager::kPairMCMap]);

   pair->SetPairTopology(values[AliReducedVarManager::kPairLxy], 0);
   pair->SetPairTopology(values[AliReducedVarManager::kPairLxyz], 1);
   //pair->SetPairTopology(values[AliReducedVarManager::kPairLxy], 2);
   //pair->SetPairTopology(values[AliReducedVarManager::kPairLxy], 3);
   pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTime], 4);
   pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTimeXYZ], 5);
   if (values[AliReducedVarManager::kPseudoProperTimeError] > 0.) {
      pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTime] / values[AliReducedVarManager::kPseudoProperTimeError], 6);
   }
   if (values[AliReducedVarManager::kPseudoProperTimeXYZError] > 0.) {
      pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTimeXYZ] / values[AliReducedVarManager::kPseudoProperTimeXYZError], 7);
   }
   pair->SetPairTopology(values[AliReducedVarManager::kPairCosPointingAngle], 8);
   pair->SetPairTopology(values[AliReducedVarManager::kPairCosPointingAngleXY], 9);
   pair->SetPairTopology(values[AliReducedVarManager::kPairDCAXY], 10);
   pair->SetPairTopology(values[AliReducedVarManager::kPairDCAZ], 11);

   if (fComputeMult) {
      FillMultiplicity(kTRUE, kTRUE, values[AliReducedVarManager::kPhi]);
      for (int icut = 0; icut<GetNMeasMultCuts(); icut++) {
         // Fill mult regions relative to jpsi
         pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward + icut], 0, icut);
         pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse + icut], 1, icut);
         pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway + icut], 2, icut);
      }
   }
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsAsymmetricDecayChannel() {
   //
   // Check if the studied decay channel is assymetric
   // NOTE: If the decay channel studied is asymmetric, it is assumed that the cuts on the two legs are different
   //            and that the 2 lists of leg cuts hold the same number of cuts
   if(fCandidateType==AliReducedPairInfo::kLambda0ToPPi) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kALambda0ToPPi) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDplusToK0sPiplus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDplusToK0sKplus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDplusToPhiPiplus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDminusToK0sPiminus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDminusToK0sKminus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDminusToPhiPiminus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDzeroToKminusPiplus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kADzeroToKplusPiminus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDsplusToK0sKplus) return kTRUE;
   if(fCandidateType==AliReducedPairInfo::kDsminusToK0sKminus) return kTRUE;
   return kFALSE;
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillCandidateLegHistograms(TString histClass, AliReducedBaseTrack* track, Float_t* values, Int_t leg, Bool_t isAsymmetricDecayChannel) {
   //
   // fill track histogram lists according to the track flags 
   //
   for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
      if(track->TestFlag(leg==2 && isAsymmetricDecayChannel ? icut+32 : icut)) {
         fHistosManager->FillHistClass(Form("%s_%s", histClass.Data(), (leg==2 && isAsymmetricDecayChannel ? fLeg2Cuts.At(icut)->GetName() : fLeg1Cuts.At(icut)->GetName())), fValues);
      }
   }
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillCandidatePairHistograms(TString histClass, AliReducedPairInfo* pair, Float_t* values, Bool_t isAsymmetricDecayChannel) {
   //
   // fill track histogram lists according to the track flags 
   //
   for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
      if(pair->TestFlag(icut)) {
         fHistosManager->FillHistClass(Form("%s_%s%s", histClass.Data(), fLeg1Cuts.At(icut)->GetName(), (isAsymmetricDecayChannel ? Form("_%s", fLeg2Cuts.At(icut)->GetName()) : "")), fValues);
      }
   }
}

//___________________________________________________________________________
const Char_t* AliReducedAnalysisFilterTrees::GetCandidateLegCutName(Int_t i, Int_t leg) {
   //
   // get candidate leg cut name
   //
   if(leg==2 && IsAsymmetricDecayChannel())
      return (i<fLeg2Cuts.GetEntries() ? fLeg2Cuts.At(i)->GetName() : "");
   return (i<fLeg1Cuts.GetEntries() ? fLeg1Cuts.At(i)->GetName() : "");
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Finish() {
  //
  // run stuff after the event loop
  //

   if(fOptionRunMixing && !fOptionRunOverMC) {
     fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
     if(fOptionRunMixingMult) { // Mixed event in multiplicity bins
        AliMixingHandler* handler;
        TIter nextHandler(&fMixingHandlerMult); 
        for (int i=0; i<fMixingHandlerMult.GetEntries(); i++) {
            handler = (AliMixingHandler*) nextHandler();
            handler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
        }
     }
   }

}

//___________________________________________________________________________
// UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track) {
Bool_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track) {
   //
   // Check a reconstructed track against all the specified MC truth cuts
   //
   // TODO: In the fLegCandidatesMCcuts one can also add AliSignalMC objects which can then be tested
   //             using the AliReducedTrackInfo::fMCPdg[]
   //
   if(!fLegCandidatesMCcuts) return 0;

   Bool_t decision = fLegCandidatesMCcuts->IsSelected(track);
   //for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
      //AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
     // if(->IsSelected(track))
      //   decisionMap |= (UInt_t(1)<<i);
   //}

   return decision;
}

//___________________________________________________________________________
//UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack) {
Bool_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack) {
   //
   // check the pair of tracks to see if they match the defined MC cuts and in addition
   // that they have the same mother
   // NOTE: The condition for the 2 tracks to have the same mother requires information on the MC label,
   //             which is available just in the full track information (AliReducedTrackInfo::fMCLabels[]).
   //           The consequence is that for the jpsi2ee analysis, the reconstructed tracks need to be always written as full tracks 
   //

   // check that both tracks are full tracks
   if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
   if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;

   // check that the tracks have the same mother
   if(TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) != TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1))) return 0;

   // check the MC requirements on each of the leg and their logical intersection
   if(!fLegCandidatesMCcuts) return 0;
   UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
   if(!pTrackDecisions) return 0;
   UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);
   return (pTrackDecisions & nTrackDecisions);
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillMCTruthHistograms() {
  //
  // fill histograms with pure signal
  //   
   // loop over the first track array
  LoopOverMCTracks(1);
  // and over the second
  // NOTE: In the current model, handling the MC truth info requires the labels, which are properties of the full track,
  //         so there is no point in looping over the second track array which, if it exists, contains just base tracks
  //LoopOverMCTracks(2);
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::LoopOverMCTracks(Int_t trackArray /*=1*/) {
   //
   // loop over the track array and check the pure MC tracks against the defined MC selections
   //   
   AliReducedTrackInfo* mother=0x0;
   AliReducedTrackInfo* daughter1 = 0x0;
   AliReducedTrackInfo* daughter2 = 0x0;

   TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if(!trackList) return;
   TIter nextTrack(trackList); 
    // if the pt dependent weights were set, check the weight and reject randomly the event
   if(fMCJpsiPtWeights) {
      for(Int_t it=0; it<trackList->GetEntries(); ++it) {
         mother = (AliReducedTrackInfo*)nextTrack();
         if(!mother->IsMCKineParticle()) continue;
         if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue;

         // apply selections on the jpsi mother
         UInt_t motherDecisions = CheckMotherMCTruth(mother);
         if(!motherDecisions) continue;
         //Apply only selections to one of the MCJpsi cuts
         if((fReweightCut >= 0) && (fReweightCut < GetNJpsiMotherMCCuts()) && !(motherDecisions & UInt_t(1)<<fReweightCut) ) continue;

         Double_t pt = mother->Pt();
         if(pt>fMCJpsiPtWeights->GetXaxis()->GetXmax())
            pt = fMCJpsiPtWeights->GetXaxis()->GetXmax();
         Double_t weight = fMCJpsiPtWeights->GetBinContent(fMCJpsiPtWeights->FindBin(pt));
         if(weight>1.0) weight = 1.0;
         Double_t rnd = gRandom->Rndm();
         if(weight<rnd) {
            fSkipMCEvent = kTRUE; 
            return;
         }
      }
   }

   nextTrack.Reset();
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
      if(!mother->IsMCKineParticle()) continue;
      if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue; // This reduces a lot running time when all primaries are stored

      // apply selections on the jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother);
      if(!motherDecisions) continue;

      // find the jpsi daughters (needed to compute 2-track properties like the polarization, etc.)
      Int_t daughter1Label = 0; Int_t daughter2Label = 0;
      FindJpsiTruthLegs(mother, daughter1Label, daughter2Label); // if mother does not decay to dielectron, daughter1Label or daughter2Label is 0
      daughter1 = FindTrackByLabel(daughter1Label, kTRUE);
      daughter2 = FindTrackByLabel(daughter2Label, kTRUE);


      // reset track variables and fill info
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
      // loop over jpsi mother selections and fill histograms before the kine cuts on electrons
      for(Int_t iCut = 0; iCut < fJpsiMotherMCcuts.GetEntries(); ++iCut) {
         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("PureMCTRUTH_BeforeSelection_%s", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
      }
      if (fMCTruthJpsi2eeOnly && ((daughter1Label == 0) || (daughter2Label == 0))) continue;      
      if (!daughter1) continue;
      if (!daughter2) continue;

      // apply selections on pure MC daughter electrons (kine cuts)
      UInt_t daughter1Decisions = CheckDaughterMCTruth(daughter1);
      if(!daughter1Decisions) continue;
      UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
      if(!daughtersDecisions) continue;
      for(Int_t iCut = 0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
         AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);

         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         if(!(daughtersDecisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("PureMCTRUTH_AfterSelection_%s", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);

         for (int jcut = 0; jcut < GetNMeasMultCuts(); jcut++)
            fHistosManager->FillHistClass(Form("JpsiPtMultCorrel_%s_%s", GetMeasMultcutName(jcut), fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);

         if ((daughter1Label == 0) || (daughter2Label == 0)) continue; // We only want Jpsi2ee and not all Jpsi

         AliReducedTrackInfo* daughter1Det = FindTrackByLabel(daughter1Label, kFALSE); //detected track
         if(!daughter1Det) continue;
         if(!IsCandidateLegSelected(daughter1Det)) continue;
         AliReducedTrackInfo* daughter2Det = FindTrackByLabel(daughter2Label, kFALSE);
         if(!daughter2Det) continue;
         if(!IsCandidateLegSelected(daughter2Det)) continue;

         AliReducedVarManager::FillPairInfo(daughter1Det, daughter2Det, fCandidateType, fValues);

         //if ((fValues[AliReducedVarManager::kMass] > 3.16) || (fValues[AliReducedVarManager::kMass] < 2.92)) continue; 
         for (int icutleg = 0; icutleg < GetNCandidateLegCuts(); icutleg++) {
            if (!daughter1Det->TestFlag(icutleg) || !daughter2Det->TestFlag(icutleg)) continue;
            fHistosManager->FillHistClass(Form("PureMCTRUTH_DetectedDaughters_%s_%s",fJpsiMotherMCcuts.At(iCut)->GetName(), GetCandidateLegCutName(icutleg, 1)), fValues);

            for (int jcut = 0; jcut < GetNMeasMultCuts(); jcut++)
               fHistosManager->FillHistClass(Form("JpsiPtMultCorrelMeas_%s_%s_%s", GetMeasMultcutName(jcut), fJpsiMotherMCcuts.At(iCut)->GetName(), GetCandidateLegCutName(icutleg, 1)), fValues);
         }
      }

      
   }  // end loop over tracks
   return;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckMotherMCTruth(AliReducedTrackInfo* mother) {
   //
   // Check the mother pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiMotherMCcuts.GetEntries()==0) return 0;

   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);

      if(cut->IsSelected(mother))
         decisionMap |= (UInt_t(1)<<i);
   }

   return decisionMap;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckDaughterMCTruth(AliReducedTrackInfo* daughter) {
   //
   // Check the daughter pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiElectronMCcuts.GetEntries()==0) return 0;

   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiElectronMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiElectronMCcuts.At(i);
      if(cut->IsSelected(daughter))
         decisionMap |= (UInt_t(1)<<i);
   }

   return decisionMap;
}

//___________________________________________________________________________
AliReducedTrackInfo* AliReducedAnalysisFilterTrees::FindTrackByLabel(Int_t label, bool isTruth) {
   //
   // search the track list for pure MC track with label and return the track pointer
   //
   AliReducedTrackInfo* track=0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(isTruth && (!track->IsMCKineParticle())) continue;
      if(!isTruth && track->IsMCKineParticle()) continue;
      if(track->MCLabel(0)==label) return track;
   }
   return 0x0;
}

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label) {
   //
   // find the jpsi legs in the list of pure MC truth particles
   //
   Int_t mLabel = mother->MCLabel(0);
   AliReducedTrackInfo* track=0x0;

   Int_t legsFound = 0;

   // loop over the first track array
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      if(legsFound==2) return;
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCKineParticle()) continue;
      if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
         legsFound += 1;
         if(legsFound==1) leg1Label = track->MCLabel(0);
         if(legsFound==2) leg2Label = track->MCLabel(0);
      }
   }
   return;
}