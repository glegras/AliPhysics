//
// Creation date: 2017/08/09
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISFILTERTREES_H
#define ALIREDUCEDANALYSISFILTERTREES_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisFilterTrees : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisFilterTrees();
  AliReducedAnalysisFilterTrees(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisFilterTrees();
  
  virtual void Init();
  virtual void Process();
  virtual void Finish();
  
  // setters


  void AddTrackCut(AliReducedInfoCut* cut) {
      // Add track cut and setup mixing handler      
      fTrackCuts.Add(cut);
      fMixingHandler->SetNParallelCuts(fMixingHandler->GetNParallelCuts()+1);
      TString histClassNames = fMixingHandler->GetHistClassNames();
      fMixingHandlerTRD->SetNParallelCuts(fMixingHandlerTRD->GetNParallelCuts()+1);
      TString histClassNamesTRD = fMixingHandlerTRD->GetHistClassNames();
      if (fPairCuts.GetEntries()>1) {
         histClassNames = "";
         histClassNamesTRD = "";
         for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
            for (Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); iTrackCut++) {
               histClassNames += Form("PairMEPP_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNames += Form("PairMEPM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNames += Form("PairMEMM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesTRD += Form("PairMEPP_%s_%s_TRD;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesTRD += Form("PairMEPM_%s_%s_TRD;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesTRD += Form("PairMEMM_%s_%s_TRD;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            }
         }
      } else {
         histClassNames += Form("PairMEPP_%s;", cut->GetName());
         histClassNames += Form("PairMEPM_%s;", cut->GetName());
         histClassNames += Form("PairMEMM_%s;", cut->GetName());
         histClassNamesTRD += Form("PairMEPP_%s_TRD;", cut->GetName());
         histClassNamesTRD += Form("PairMEPM_%s_TRD;", cut->GetName());
         histClassNamesTRD += Form("PairMEMM_%s_TRD;", cut->GetName());
      }
      fMixingHandler->SetHistClassNames(histClassNames.Data());
      fMixingHandlerTRD->SetHistClassNames(histClassNamesTRD.Data());
   }


  void AddPairCut(AliReducedInfoCut* cut) {
      // Add pair cut and setup mixing handler
      fPairCuts.Add(cut);  
      fMixingHandler->SetNParallelPairCuts(fMixingHandler->GetNParallelPairCuts()+1);
      fMixingHandlerTRD->SetNParallelPairCuts(fMixingHandlerTRD->GetNParallelPairCuts()+1);
      if (fPairCuts.GetEntries()>1) {
         TString histClassNamesNew = "";
         TString histClassNamesNewTRD = "";
         for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
            for (Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); iTrackCut++) {
               histClassNamesNew += Form("PairMEPP_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesNew += Form("PairMEPM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesNew += Form("PairMEMM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesNewTRD += Form("PairMEPP_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesNewTRD += Form("PairMEPM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
               histClassNamesNewTRD += Form("PairMEMM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            }
         }
         fMixingHandler->SetHistClassNames(histClassNamesNew.Data());
         fMixingHandlerTRD->SetHistClassNames(histClassNamesNewTRD.Data());
      }
   }


  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void SetWriteFilteredTracks(Bool_t option=kTRUE) {fWriteFilteredTracks=option;}
  void SetWriteFilteredPairs(Bool_t option=kTRUE) {fWriteFilteredPairs=option;}
  void SetWriteFilteredTracksCandidatesOnly(Bool_t option=kTRUE) {fWriteFilteredTracksCandidatesOnly=option;}
  void SetRejectEmptyEvents(Bool_t option=kTRUE) {fRejectEmptyEvents=option;}
  void SetMCJpsiPtWeights(TH1F* weights) {fMCJpsiPtWeights = weights;}
  void SetReweightCut(Int_t ncut){fReweightCut = ncut;}
  void SetMCParticleCompositionWeights(TH3F* weights) {fMCParticleCompositionWeights = weights;}
  
  void SetBuildCandidatePairs(AliReducedPairInfo::CandidateType type) {fBuildCandidatePairs=kTRUE; fCandidateType=type;}
  void SetBuildCandidateLikePairs(Bool_t option=kTRUE) {fBuildCandidateLikePairs=option;}
  void AddCandidateLeg1Cut(AliReducedInfoCut* cut) {fLeg1Cuts.Add(cut);}
  void AddCandidateLeg2Cut(AliReducedInfoCut* cut) {fLeg2Cuts.Add(cut);}
  void AddCandidatePairCut(AliReducedInfoCut* cut) {fCandidatePairCuts.Add(cut);}
  void SetRunCandidatePrefilter(Bool_t option=kTRUE) {fRunCandidatePrefilter=option;}
  void SetRunCandidatePrefilterOnSameCharge(Bool_t option=kTRUE) {fRunCandidatePrefilterOnSameCharge=option;}
  void SetRunEventMixing(Bool_t option) {fOptionRunMixing = option;}
  void SetRunEventMixingMult(Bool_t option) {fOptionRunMixingMult = option;}
  void SetComputeMult(Bool_t option) {fComputeMult = option;}
  void SetReweightParticleComposition(Bool_t option) {fReweightParticleComposition = option;}
  void SetSharePCCWeightsBetweenEvents(Bool_t option) {fSharePCCWeightsBetweenEvents = option;}

  void AddMeasuredMultTrackCut(AliReducedInfoCut* cut, TH2F* hWeights = nullptr) {
   fMeasuredMultTrackCuts.Add(cut);
   if (hWeights) fWeightsTrackCuts.Add(hWeights);
   else {
      TH2F* hWeightsTrackCuts = new TH2F(Form("hWeightsTrackCuts_%s", cut->GetName()), "weights", 1, 0, 3e5, 1, 0, 1e3);
      hWeightsTrackCuts->SetBinContent(1, 1, 1.);
      fWeightsTrackCuts.Add(hWeightsTrackCuts);
   }
  }
  void AddTrueMultTrackCut(AliReducedInfoCut* cut) {fTrueMultTrackCuts.Add(cut);}
  void AddCandidateLeg1PrefilterCut(AliReducedInfoCut* cut) {fLeg1PrefilterCuts.Add(cut);}
  void AddCandidateLeg2PrefilterCut(AliReducedInfoCut* cut) {fLeg2PrefilterCuts.Add(cut);}
  void AddCandidateLeg1PairPrefilterCut(AliReducedInfoCut* cut) {fLeg1PairPrefilterCuts.Add(cut);}
  void AddCandidateLeg2PairPrefilterCut(AliReducedInfoCut* cut) {fLeg2PairPrefilterCuts.Add(cut);}

  void AddMixingHandler(AliMixingHandler* handler) {fMixingHandlerMult.Add(handler);}
  void SetMultBinsMixing(int nBins, Float_t* bins) {fNMultBinsMixing = nBins; for(int i = 0;i<nBins+1;i++) fMultBinsMixing[i] = bins[i];}
  void SetJpsiMassDist(TF1* massDist) {fJpsiMassDist = massDist;}
  void SetMCTruthJpsi2eeOnly(Bool_t option) {fMCTruthJpsi2eeOnly = option;}
  void SetRegionsToMCTruth(Bool_t option) {fRegionsToMCTruth = option;}
  void SetDefaultRandomPhi(Bool_t option) {fDefaultRandomPhi = option;}
  void SetMinPtLeading(float minpt) {fMinPtLeading = minpt;}
  void SetVertexCorrection(Bool_t option) {fVtxCorrection = option;}
  
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;};
  // getters
  virtual AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  virtual AliMixingHandler* GetMixingHandler() const {return fMixingHandler;}
  virtual AliMixingHandler* GetMixingHandlerTRD() const {return fMixingHandlerTRD;}
  Bool_t GetWriteFilteredTracks() const {return fWriteFilteredTracks;}
  Int_t GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t* GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");} 
  Bool_t GetWriteFilteredPairs() const {return fWriteFilteredPairs;}
  Bool_t GetRejectEmptyEvents() const {return fRejectEmptyEvents;}
  Int_t GetNPairCuts() const {return fPairCuts.GetEntries();}
  const Char_t* GetPairCutName(Int_t i) const {return (i<fPairCuts.GetEntries() ? fPairCuts.At(i)->GetName() : "");} 
  Bool_t GetBuildCandidatePairs() const {return fBuildCandidatePairs;}
  Bool_t GetBuildCandidateLikePairs() const {return fBuildCandidateLikePairs;}
  Bool_t GetComputeMult() const {return fComputeMult;};
  Bool_t GetRunEventMixing() const {return fOptionRunMixing;}
  Bool_t GetRunEventMixingMult() const {return fOptionRunMixingMult;}
  Int_t GetNMultBinsMixing() const {return (fOptionRunMixingMult ? fNMultBinsMixing : 0);}
  Float_t GetMultBinsMixing(int i) const {return (fOptionRunMixingMult ? fMultBinsMixing[i] : 0);}
  Int_t GetCandidateType() const {return fCandidateType;}
  Bool_t GetRunCandidatePrefilter() const {return fRunCandidatePrefilter;}
  Int_t GetNCandidateLegCuts() const {return fLeg1Cuts.GetEntries();}
  const Char_t* GetCandidateLegCutName(Int_t i, Int_t leg);
  Bool_t IsAsymmetricDecayChannel();
  Bool_t GetRunOverMC() const {return fOptionRunOverMC;};
  //  Int_t GetNLegCandidateMCcuts() const {return fLegCandidatesMCcuts.GetEntries();}
  //const Char_t* GetLegCandidateMCcutName(Int_t i) const {return (i<fLegCandidatesMCcuts.GetEntries() ? fLegCandidatesMCcuts.At(i)->GetName() : "");}
  const Char_t* GetLegCandidateMCcutName() const {return (fLegCandidatesMCcuts ? fLegCandidatesMCcuts->GetName() : "");}
  Int_t GetNJpsiMotherMCCuts() const {return fJpsiMotherMCcuts.GetEntries();}
  const Char_t* GetJpsiMotherMCcutName(Int_t i) const {return (i<fJpsiMotherMCcuts.GetEntries() ? fJpsiMotherMCcuts.At(i)->GetName() : "");}
  Int_t GetNMeasMultCuts() const {return fMeasuredMultTrackCuts.GetEntries();}
  const Char_t* GetMeasMultcutName(Int_t i) const {return (i<fMeasuredMultTrackCuts.GetEntries() ? fMeasuredMultTrackCuts.At(i)->GetName() : "");}
  TF1* GetJpsiMassDist() const {return fJpsiMassDist;}
  Bool_t GetMCTruthJpsi2eeOnly() const {return fMCTruthJpsi2eeOnly;}
  Bool_t GetReweightParticleComposition() const {return fReweightParticleComposition;}
  Bool_t GetVertexCorrection () const {return fVtxCorrection;}
 
 /* void AddLegCandidateMCcut(AliReducedInfoCut* cut) {
     if(fLegCandidatesMCcuts.GetEntries()>=32) return;
     fLegCandidatesMCcuts.Add(cut);
  }*/
 // NOTE: The MC truth selection works with just one MC truth cut. It is not implemented properly for asymmetric decay channels,
 //         just one MC selection is applied to both legs.
 //       In the case that an MC truth selection is applied, then only built pairs which fulfill the MC truth will be written in the filtered trees.
 void SetLegCandidateMCcut(AliReducedInfoCut* cut) {
     //if(fLegCandidatesMCcuts.GetEntries()>=32) return;
     fLegCandidatesMCcuts  = cut;
  }

  void AddJpsiMotherMCCut(AliReducedInfoCut* cutMother, AliReducedInfoCut* cutElectron) {
     if(fJpsiMotherMCcuts.GetEntries()>=32) return;
     fJpsiMotherMCcuts.Add(cutMother);
     fJpsiElectronMCcuts.Add(cutElectron);
  }

  void FillMCTruthHistograms();
 
protected:
   std::map<int, std::vector<int>> fWeightsTrue; // Save some weights to be used in next events
   std::map<int, std::vector<int>> fWeightsMeas; // Save some weights to be used in next events

   AliHistogramManager* fHistosManager;   // Histogram manager
   AliMixingHandler*    fMixingHandler;       // mixing handler
   AliMixingHandler*    fMixingHandlerTRD;       // mixing handler for TRD (event mixing from MB and TRD gives different results)
   TList                fMixingHandlerMult;       // mixing handlers in multiplicity bins
   Float_t                fMultBinsMixing[1000];  // Multiplicity bins for mixing in multiplicity bins
   Int_t                fNMultBinsMixing;   // Number of multiplicity bins for mixing in multiplicity bins
   
   TList fEventCuts;               // array of event cuts used for filtering
   TList fTrackCuts;               // array of track cuts used for filtering
   Bool_t fWriteFilteredTracks;   // filter the track list
   Bool_t fWriteFilteredTracksCandidatesOnly; // Write only the tracks associated to a pair candidate
   TList fPairCuts;                  // array of pair cuts used for filtering
   Bool_t fWriteFilteredPairs;   // filter the pair list
   Bool_t fRejectEmptyEvents;     // if true, do not write events without tracks or pairs
   Bool_t fMCTruthJpsi2eeOnly;    //if true, store only MCtruth for Jpsi which decay into dielectron
   Bool_t fRegionsToMCTruth;    //if true, the regions are defined relative to true particles, else they are defined relative to measured tracks
   Bool_t fDefaultRandomPhi;    //if true, when regions are calculated relative to Jpsi and there is no Jpsi, we choose a random phi, else we take phi leading
   Float_t fMinPtLeading;
   Bool_t fVtxCorrection;

   TF1* fJpsiMassDist;            // Real jpsi mass distribution (Crystal-ball)
   Bool_t fSharePCCWeightsBetweenEvents;
   Bool_t fReweightParticleComposition;
   TH3F* fMCParticleCompositionWeights;

   Bool_t fComputeMult;            //if true, count the tracks to compute true and measured multiplicity
   Bool_t fBuildCandidatePairs;   // if true, build additional candidate pairs from selected tracks 
   Bool_t fBuildCandidateLikePairs;  // if true, build also like pairs (e.g. like-sign for symmetric decay channels)
   Int_t fCandidateType;             // candidate type, see AliReducedPairInfo::CandidateType 
   TList fLeg1Cuts;                      // list of track cuts for LEG1  (these cuts will be used also for LEG2 if the decay channel is symmetric)
   TList fLeg2Cuts;                      // list of tracks cuts for LEG2 (NOTE: fLeg1Cuts and fLeg2Cuts must contain the same number of cuts)
   TList fCandidatePairCuts;          // list of cuts for pair candidates
   Bool_t fRunCandidatePrefilter;   // if true, run a prefilter on the selected legs
   Bool_t fRunCandidatePrefilterOnSameCharge;   // default FALSE (unlike charged pairs only);
                                                                // if true, run the prefilter on same charge pairs also;
   
   TList fMeasuredMultTrackCuts;            // cuts for tracks for determination of measured multiplicity
   TList fWeightsTrackCuts;             // histograms giving weights applied to tracks, to smear the efficiency, vs pt vs run number 
   TList fTrueMultTrackCuts;            // cuts for MC tracks for determination of true multiplicity
   TList fLeg1PrefilterCuts;            // cuts for tracks used in the prefilter for LEG1  
   TList fLeg2PrefilterCuts;            // cuts for tracks used in the prefilter for LEG2
   TList fLeg1PairPrefilterCuts;      // cuts on the prefilter pairs for LEG1
   TList fLeg2PairPrefilterCuts;      // cuts on the prefilter pairs for LEG2
   TList fLeg1Tracks;                    // list of selected LEG1 tracks in the current event
   TList fLeg2Tracks;                    // list of selected LEG2 tracks in the current event
   TList fLeg1PrefilteredTracks;    // list of prefilter selected LEG1 tracks in the current event
   TList fLeg2PrefilteredTracks;    // list of prefilter selected LEG2 tracks in the current event

   Bool_t fOptionRunMixingMult;
   Bool_t fOptionRunMixing;    // true: run event mixing, false: no event mixing
   Bool_t fOptionRunOverMC;  // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
   // selection based on the MC truth information of the reconstructed leg candidates
   // NOTE:    The list is a list of AliReducedInfoCut objects which can be used to 
   //              apply cuts on the MC flags of the tracks.
   // NOTE: The names of the cuts are used in the naming of the histogram classes
   AliReducedInfoCut *fLegCandidatesMCcuts; 
 
   // selection cuts for the pure MC truth (select the J/psi from stack)
   // the list should contains cuts which can be applied to a pure MC truth particle (no reconstructed information)
   //  e.g. cuts on the MC flags and on kinematics
   //  For each selection, a separate histogram directory will be created
   TList fJpsiMotherMCcuts;
   TH1F*  fMCJpsiPtWeights;            //! weights vs pt to reject events depending on the jpsi true pt (needed to re-weights jpsi Pt distribution)
   Int_t fReweightCut;                //The number of the MC cut on which the weights should be applied (ex: only on prompt Jpsi)
   Bool_t fSkipMCEvent; // if true MC event is skipped
   // Selection on the MC truth of the electrons from the jpsi decay
   //  Tipically, here one can specify the kinematic selection on the electrons from jpsi decay
   //       so dividing the jpsi yield at this step by the yield of jpsi selected by the fJpsiMotherMCcuts, one can obtain the
   //       acceptance efficiency.
   //  NOTE: The number of selections on the jpsi electron needs to be the same and in sync with the number of fJpsiMotherMCcuts cuts
   TList fJpsiElectronMCcuts;
   
   Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
   Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
   Bool_t IsPairSelected(AliReducedPairInfo* pair, Float_t* values=0x0);
   void CreateFilteredEvent();
   Bool_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack);
 Bool_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* track);
  void    FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label);
  AliReducedTrackInfo* FindTrackByLabel(Int_t label, bool isTruth = true);
  void    LoopOverMCTracks(Int_t trackArray =1);
  UInt_t CheckMotherMCTruth(AliReducedTrackInfo* mother);
  UInt_t CheckDaughterMCTruth(AliReducedTrackInfo* daughter); 

   void FillMultiplicity(Bool_t regions = kFALSE, bool wrtPhiRef = kFALSE, float phiRef = -999.);
   Float_t GetParticleWeight(AliReducedTrackInfo* track);
   Bool_t TrackIsCandidateLeg(AliReducedBaseTrack* track);
   void WriteFilteredPairs();
   void WriteFilteredTracks(Int_t array=1);
   Bool_t IsTrackMeasuredMultSelected(AliReducedBaseTrack* track, Float_t* values=0x0); 
   Bool_t IsTrackTrueMultSelected(AliReducedBaseTrack* track, Float_t* values=0x0); 
   Bool_t IsCandidateLegSelected(AliReducedBaseTrack* track, Float_t* values=0x0, Int_t whichLeg=1); 
   Bool_t IsCandidatePairSelected(Float_t* values);
   Bool_t IsCandidateLegPrefilterSelected(AliReducedBaseTrack* track, Float_t* values=0x0, Int_t whichLeg=1);
   Bool_t IsCandidateLegPairPrefilterSelected(Float_t* values, Int_t whichLeg=1);
   void BuildCandidatePairs();
   void RunCandidateLegsSelection(Int_t arrayOption /*=1*/);
   void RunCandidateLegsPrefilter(Int_t leg);
   void RunSameEventPairing();
   void SetupPair(AliReducedPairInfo* pair, Float_t* values);
   ULong_t CheckTrackCompatibility(AliReducedBaseTrack* leg1, AliReducedBaseTrack* leg2, Bool_t isAsymmetricDecayChannel);
   void FillCandidateLegHistograms(TString histClass, AliReducedBaseTrack* track, Float_t* values, Int_t leg, Bool_t isAsymmetricDecayChannel);
   void FillCandidatePairHistograms(TString histClass, AliReducedPairInfo* pair, Float_t* values, Bool_t isAsymmetricDecayChannel);

  ClassDef(AliReducedAnalysisFilterTrees,2);
};

#endif
