#include <iostream>
#include <fstream>
#include <string>

#include "TF1.h"


void Setup(AliReducedAnalysisFilterTrees* processor, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/, Bool_t isInjected);
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/,  TString prod /*="LHC10h"*/);
void SetupMixingHandler(AliReducedAnalysisFilterTrees* task);
void DefineHistograms(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/);


//__________________________________________________________________________________________
AliAnalysisTask* AddTask_glegras_FilterTrees(Bool_t isAliRoot=kTRUE, Int_t runMode=1, Bool_t isMC = kTRUE, Bool_t isInjected = kFALSE, Bool_t ispp = kTRUE, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_glegras_FilterTrees(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisFilterTrees* filterTask = new AliReducedAnalysisFilterTrees("FilterTrees","filter DST trees");
  filterTask->Init();
  filterTask->SetFilteredTreeWritingOption(AliReducedAnalysisTaskSE::kBaseEventsWithFullTracks);
  filterTask->SetWriteFilteredTracks(kTRUE);
  filterTask->SetWriteFilteredTracksCandidatesOnly(kTRUE);
  filterTask->SetWriteFilteredPairs(kFALSE);
  filterTask->SetBuildCandidatePairs(AliReducedPairInfo::kJpsiToEE);
  filterTask->SetBuildCandidateLikePairs(kTRUE);
  filterTask->SetRunCandidatePrefilter(kTRUE);
  filterTask->SetRunCandidatePrefilterOnSameCharge(kFALSE);
  filterTask->SetRejectEmptyEvents(kTRUE);
  filterTask->SetRunEventMixing(kTRUE);
  filterTask->SetComputeMult(kTRUE);
  filterTask->SetRunOverMC(isMC);
  filterTask->SetMCTruthJpsi2eeOnly(kFALSE);
  filterTask->SetReweightParticleComposition(kTRUE);
  filterTask->SetRegionsToMCTruth(kFALSE);
  filterTask->SetDefaultRandomPhi(kTRUE);
  if (isMC && isInjected) {
    bool fDoImpParCorr = kTRUE;
    // Look where we are: maf or local
    TFile* file = new TFile("/home/glegras/alice/data/ImpParCorr/16k/ITSgraphs_Current.root");
    TGraph* gTest; file->GetObject("D0RPResE", gTest);
    if (gTest) {// We are local
      AliReducedVarManager::SetDoImpParCorr(fDoImpParCorr, "/home/glegras/alice/data/ImpParCorr/");
    }
    else { // We are probably in maf
      AliReducedVarManager::SetDoImpParCorr(fDoImpParCorr, "/gluster1/glegras/ImpParCorr/");
    }
    file->Close();
  }

  Setup(filterTask, ispp, prod, isInjected);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisFilterTrees such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode, kTRUE);
  task->AddTask(filterTask);
  
  if(isAliRoot){
     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) {
        Error("AddTask_iarsene_dst", "No analysis manager found.");
        return 0;
     }
     
     AliAnalysisDataContainer* cReducedEvent = NULL;
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
       printf("INFO on AddTask_iarsene_FilterTrees(): use on the fly events ");
       cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
       if(!cReducedEvent) {
         printf("ERROR: In AddTask_iarsene_FilterTrees(), couldn't find exchange container with ReducedEvent! ");
         printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
         return 0x0;
       }
     }
            
     mgr->AddTask(task);
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
        mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
       mgr->ConnectInput(task, 0, cReducedEvent);
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("filterQAhistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
     mgr->ConnectOutput(task, 1, cOutputHist );
     
     AliAnalysisDataContainer *cOutputTree = mgr->CreateContainer("filteredDstTree", TTree::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstTree.root");
     mgr->ConnectOutput(task, 2, cOutputTree );
  }
  else {
    // nothing at the moment  

    
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

      AliAnalysisDataContainer* cReducedEvent = NULL;

      cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");

      mgr->AddTask(task);

      mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());

      AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("Histos", THashList::Class(),
                     AliAnalysisManager::kOutputContainer, "AnalysisHistograms_Jpsi2ee.root");
      

      AliAnalysisDataContainer *cOutputTree = mgr->CreateContainer("filteredDstTree", TTree::Class(),
                             AliAnalysisManager::kOutputContainer, "JpsiCandidates.root");
     
      mgr->ConnectOutput(task, 1, cOutputHist ); 
      mgr->ConnectOutput(task, 2, cOutputTree ); 


  }
  
  return task;
}

//______________________________________________________________
Double_t V0cut(Double_t *x, Double_t *par){
    double runNo = x[0];

    if (runNo<=254604) return 560.;
    if (runNo<=255079) return 530.;
    if (runNo<=256504) return 500.;
    if (runNo<=257850) return 450.;
    if (runNo<=264076) return 400.;
    if (runNo<=272151) return 570.;
    if (runNo<=273824) return 480.;
    if (runNo<=275515) return 540.;
    if (runNo<=279879) return 520.;
    if (runNo<=285980) return 500.;
    if (runNo<=290323) return 580.;
    if (runNo<=292012) return 570.;
    if (runNo<=292192) return 560.;
    return 540.;

}


//_________________________________________________________________
void Setup(AliReducedAnalysisFilterTrees* processor, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/, Bool_t isInjected) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  // 
  
  bool isMC = processor->GetRunOverMC();

  bool correctTPCNsig = kTRUE;
  if (isMC) correctTPCNsig = kFALSE;

  TF1* fV0cut = new TF1("fV0cut", V0cut, 2e5, 3e5);

  TF1* NtracksCutMC = new TF1("NtracksCutMC", "0.558*x - 0.43 + 5*(0.564 * sqrt(x) - 0.14)", 0, 200); // (mean + 5 sigma) cut - for new mult estimator vs vtx contributors
  TF1* NtracksCut = new TF1("NtracksCut", "0.541*x - 0.18 + 5*(0.614 * sqrt(x) - 0.50)", 0, 200);

  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  //For now, vtxz and vertex requirement cuts are dealt with in the analysis task
  //evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  if (isMC && !isInjected) evCut1->AddCut(AliReducedVarManager::kVtxZDeltaTrue, -0.5, 0.5);
  if (isMC && !isInjected) evCut1->AddCut(AliReducedVarManager::kNGlobalTracks+1, 0., NtracksCutMC, kFALSE, AliReducedVarManager::kNVtxContributors, 0., 1e3);
  if (!isMC) evCut1->AddCut(AliReducedVarManager::kNGlobalTracks+1, 0., NtracksCut, kFALSE, AliReducedVarManager::kNVtxContributors, 0., 1e3);
  //evCut1->AddCut(AliReducedVarManager::kNVtxContributors, -0.5, 0.5, kTRUE);
  //evCut1->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  //evCut1->AddEventTagFilterBit(14);
  if (!isMC) evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, 0., fV0cut, kTRUE, AliReducedVarManager::kRunNo, 0., 0., kTRUE, AliReducedVarManager::kHighMultV0Triggered, 0.5, 1.5);
  processor->AddEventCut(evCut1);

  processor->SetMinPtLeading(5.0);

  if (isMC && isInjected) {
    TF1* fPtJpsi = new TF1("fPtJpsi","[0]*x/pow(1+(x/[1])*(x/[1]),[2])", 0., 30.);
    //fPtJpsi->SetParameters(1, 4.09, 3.04); // arxiv:2108.01906 - standard -inclusive spectrum
    fPtJpsi->SetParameters(1, 4.68, 3.56); // ML - prompt spectrum
    TFile* filePtWeights = new TFile("/gluster1/glegras/InjectedJpsiPtWeights.root","read");
    TH1F* hpt;
    //const char* histName = "Pt_inclusive"; // standard
    const char* histName = "Pt_prompt"; // Only if ML
    filePtWeights->GetObject(histName, hpt);
    if(!hpt) {
       filePtWeights = new TFile("~/alice/AliPhysics/PWGDQ/reducedTree/macros/InjectedJpsiPtWeights.root","read");
       filePtWeights->GetObject(histName, hpt);
    }
    TH1F* hPtWeights = new TH1F("hPtWeights", "Weights for rescaling injected Jpsi", hpt->GetNbinsX(),hpt->GetBinLowEdge(1),hpt->GetBinCenter(hpt->GetNbinsX())+hpt->GetBinCenter(hpt->GetNbinsX())/2-hpt->GetBinLowEdge(hpt->GetNbinsX())/2);
    for (int n = 1; n < hpt->GetNbinsX()+1; n++) hPtWeights->SetBinContent(n, fPtJpsi->Eval(hPtWeights->GetBinCenter(n))/hpt->GetBinContent(n));
    hPtWeights->Scale(1./hPtWeights->GetMaximum());

    processor->SetMCJpsiPtWeights(hPtWeights); //only if MC and injected Jpsi
  }

  // --------------------------------- Set track cuts for multiplicity info ------------------------------------------

  // Measured track cuts (max. 8 cuts) - names should not be included one in the other
  // Standard cut
  AliReducedTrackCut* measMultCut = new AliReducedTrackCut("standardCut","Cut for measured mult");
  measMultCut->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCut->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  /*if(isMC && !isInjected) { // old cuts from pass1
      measMultCut->SetTrackFilterBit(8); // Hadron standard
      measMultCut->SetTrackFilterBit(16); // TPC chi2
  }
  else*/
  measMultCut->SetTrackFilterBit(3); // old cuts
  //processor->AddMeasuredMultTrackCut(measMultCut);

  // New Standard cut
  AliReducedTrackCut* newMultCut = new AliReducedTrackCut("newMultCut","Cut for measured mult");
  newMultCut->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  newMultCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  newMultCut->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  //if(isMC && !isInjected) { // old cuts from pass1
  //  newMultCut->SetTrackFilterBit(19);
  //}
  //else
  newMultCut->SetTrackFilterBit(4);
  //processor->AddMeasuredMultTrackCut(newMultCut);

  /*AliReducedTrackCut* looseDCAXY = new AliReducedTrackCut("looseDCAXY","Cut for measured mult");
  looseDCAXY->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseDCAXY->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseDCAXY->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseDCAXY->SetTrackFilterBit(24);
  processor->AddMeasuredMultTrackCut(looseDCAXY);

  AliReducedTrackCut* tightDCAXY = new AliReducedTrackCut("tightDCAXY","Cut for measured mult");
  tightDCAXY->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightDCAXY->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightDCAXY->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightDCAXY->SetTrackFilterBit(19);
  processor->AddMeasuredMultTrackCut(tightDCAXY);

  AliReducedTrackCut* looseTPCConstrainedGlobal = new AliReducedTrackCut("looseTPCConstrainedGlobal","Cut for measured mult");
  looseTPCConstrainedGlobal->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseTPCConstrainedGlobal->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseTPCConstrainedGlobal->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseTPCConstrainedGlobal->SetTrackFilterBit(31);
  processor->AddMeasuredMultTrackCut(looseTPCConstrainedGlobal);

  AliReducedTrackCut* tightTPCConstrainedGlobal = new AliReducedTrackCut("tightTPCConstrainedGlobal","Cut for measured mult");
  tightTPCConstrainedGlobal->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightTPCConstrainedGlobal->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightTPCConstrainedGlobal->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightTPCConstrainedGlobal->SetTrackFilterBit(30);
  tightTPCConstrainedGlobal->SetTrackFilterBit(4);
  processor->AddMeasuredMultTrackCut(tightTPCConstrainedGlobal);

  AliReducedTrackCut* noSPDany = new AliReducedTrackCut("noSPDany","Cut for measured mult");
  noSPDany->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  noSPDany->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  noSPDany->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  noSPDany->SetTrackFilterBit(25);
  noSPDany->SetTrackFilterBit(4);
  noSPDany->SetUseANDonTrackFilterMap(kFALSE);
  processor->AddMeasuredMultTrackCut(noSPDany);*/


  AliReducedTrackCut* looseDCAZ = new AliReducedTrackCut("looseDCAZ","Cut for measured mult");
  looseDCAZ->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseDCAZ->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseDCAZ->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseDCAZ->SetTrackFilterBit(16);
  processor->AddMeasuredMultTrackCut(looseDCAZ);

  AliReducedTrackCut* tightDCAZ = new AliReducedTrackCut("tightDCAZ","Cut for measured mult");
  tightDCAZ->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightDCAZ->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightDCAZ->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightDCAZ->SetTrackFilterBit(15);
  processor->AddMeasuredMultTrackCut(tightDCAZ);

  AliReducedTrackCut* looseGeoLength = new AliReducedTrackCut("looseGeoLength","Cut for measured mult");
  looseGeoLength->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseGeoLength->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseGeoLength->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseGeoLength->SetTrackFilterBit(26);
  processor->AddMeasuredMultTrackCut(looseGeoLength);

  AliReducedTrackCut* tightGeoLength = new AliReducedTrackCut("tightGeoLength","Cut for measured mult");
  tightGeoLength->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightGeoLength->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightGeoLength->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightGeoLength->SetTrackFilterBit(27);
  processor->AddMeasuredMultTrackCut(tightGeoLength);

  AliReducedTrackCut* looseDeadZone = new AliReducedTrackCut("looseDeadZone","Cut for measured mult");
  looseDeadZone->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseDeadZone->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseDeadZone->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseDeadZone->SetTrackFilterBit(29);
  processor->AddMeasuredMultTrackCut(looseDeadZone);

  AliReducedTrackCut* tightDeadZone = new AliReducedTrackCut("tightDeadZone","Cut for measured mult");
  tightDeadZone->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightDeadZone->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightDeadZone->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightDeadZone->SetTrackFilterBit(28);
  processor->AddMeasuredMultTrackCut(tightDeadZone);

  AliReducedTrackCut* Pt01 = new AliReducedTrackCut("Pt01","Cut for measured mult");
  Pt01->AddCut(AliReducedVarManager::kPt, 0.1, 1e5);
  Pt01->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  Pt01->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  Pt01->SetTrackFilterBit(4);
  processor->AddMeasuredMultTrackCut(Pt01);

  AliReducedTrackCut* Pt02 = new AliReducedTrackCut("Pt02","Cut for measured mult");
  Pt02->AddCut(AliReducedVarManager::kPt, 0.2, 1e5);
  Pt02->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  Pt02->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  Pt02->SetTrackFilterBit(4);
  processor->AddMeasuredMultTrackCut(Pt02);


  /*AliReducedTrackCut* looseChi2ITS = new AliReducedTrackCut("looseChi2ITS","Cut for measured mult");
  looseChi2ITS->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseChi2ITS->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseChi2ITS->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseChi2ITS->SetTrackFilterBit(14);
  processor->AddMeasuredMultTrackCut(looseChi2ITS);

  AliReducedTrackCut* tightChi2ITS = new AliReducedTrackCut("tightChi2ITS","Cut for measured mult");
  tightChi2ITS->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightChi2ITS->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightChi2ITS->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightChi2ITS->SetTrackFilterBit(13);
  processor->AddMeasuredMultTrackCut(tightChi2ITS);

  AliReducedTrackCut* looseSharedClusters = new AliReducedTrackCut("looseSharedClusters","Cut for measured mult");
  looseSharedClusters->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseSharedClusters->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseSharedClusters->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseSharedClusters->SetTrackFilterBit(6);
  processor->AddMeasuredMultTrackCut(looseSharedClusters);

  AliReducedTrackCut* tightSharedClusters = new AliReducedTrackCut("tightSharedClusters","Cut for measured mult");
  tightSharedClusters->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightSharedClusters->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightSharedClusters->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightSharedClusters->SetTrackFilterBit(5);
  processor->AddMeasuredMultTrackCut(tightSharedClusters);

  AliReducedTrackCut* looseCrossedRows = new AliReducedTrackCut("looseCrossedRows","Cut for measured mult");
  looseCrossedRows->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseCrossedRows->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseCrossedRows->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseCrossedRows->SetTrackFilterBit(7);
  processor->AddMeasuredMultTrackCut(looseCrossedRows);

  AliReducedTrackCut* tightCrossedRows = new AliReducedTrackCut("tightCrossedRows","Cut for measured mult");
  tightCrossedRows->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightCrossedRows->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightCrossedRows->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightCrossedRows->SetTrackFilterBit(10);
  processor->AddMeasuredMultTrackCut(tightCrossedRows);

  AliReducedTrackCut* looseChi2TPC = new AliReducedTrackCut("looseChi2TPC","Cut for measured mult");
  looseChi2TPC->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  looseChi2TPC->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseChi2TPC->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  looseChi2TPC->SetTrackFilterBit(12);
  processor->AddMeasuredMultTrackCut(looseChi2TPC);

  AliReducedTrackCut* tightChi2TPC = new AliReducedTrackCut("tightChi2TPC","Cut for measured mult");
  tightChi2TPC->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  tightChi2TPC->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  tightChi2TPC->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  tightChi2TPC->SetTrackFilterBit(11);
  processor->AddMeasuredMultTrackCut(tightChi2TPC);*/



  // old Track quality cut 1
  AliReducedTrackCut* measMultCutQuality1 = new AliReducedTrackCut("Quality1","Cut for measured mult");
  measMultCutQuality1->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCutQuality1->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutQuality1->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  //if(isMC) {
  //  measMultCutQuality1->SetTrackFilterBit(8); // Hadron standard
  //  measMultCutQuality1->SetTrackFilterBit(16); // TPC chi2
  //  measMultCutQuality1->SetTrackFilterBit(17); // Quality cut 1
  //} 
  //else {
    measMultCutQuality1->SetTrackFilterBit(7); // Hadron standard 
    measMultCutQuality1->SetTrackFilterBit(15); // TPC Nclusters 
    measMultCutQuality1->SetTrackFilterBit(16); // TPC chi2
    measMultCutQuality1->SetTrackFilterBit(17); // Quality cut 1
  //}
  //measMultCut->SetTrackFilterBit(6); //for gsi data
  //processor->AddMeasuredMultTrackCut(measMultCutQuality1);

  // old data standard cut
  AliReducedTrackCut* measMultCutOld = new AliReducedTrackCut("standardCut","Cut for measured mult");
  measMultCutOld->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCutOld->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutOld->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  measMultCutOld->SetTrackFilterBit(7); // old cuts
  measMultCutOld->SetTrackFilterBit(15); // old cuts
  measMultCutOld->SetTrackFilterBit(16); // old cuts
  //processor->AddMeasuredMultTrackCut(measMultCutOld);

  // old Track quality cut 2
  AliReducedTrackCut* measMultCutQuality2 = new AliReducedTrackCut("newMultCut","Cut for measured mult");
  measMultCutQuality2->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCutQuality2->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutQuality2->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  //if(isMC) {
  //  measMultCutQuality2->SetTrackFilterBit(8); // Hadron standard
  //  measMultCutQuality2->SetTrackFilterBit(16); // TPC chi2
  //  measMultCutQuality2->SetTrackFilterBit(20); // Quality cut 2
  //}
  //else {
    measMultCutQuality2->SetTrackFilterBit(7); // Hadron standard 
    //measMultCutQuality2->SetTrackFilterBit(15); // TPC Nclusters 
    //measMultCutQuality2->SetTrackFilterBit(16); // TPC chi2
    measMultCutQuality2->SetTrackFilterBit(19); // Quality cut 2
  //}
  //measMultCut->SetTrackFilterBit(6); //for gsi data
  //processor->AddMeasuredMultTrackCut(measMultCutQuality2);


  // MC cut on true multiplicity  
  AliReducedTrackCut* trueMultCut = new AliReducedTrackCut("trueMult","Cut for true mult");
  trueMultCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  //trueMultCut->AddCut(AliReducedVarManager::kPt, 0.15, 100);
  trueMultCut->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  trueMultCut->SetRejectPureMC(kFALSE);
  //trueMultCut->SetMCFilterBit(13); // physical primaries
  trueMultCut->SetMCFilterBit(6);
  processor->AddTrueMultTrackCut(trueMultCut);




  // ----------------------------------------- Set track cuts on electrons ---------------------------------------

  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kPt, 1., 100.);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    standardCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    standardCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, 3., 30000.0); //standard
    standardCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, 3., 30000.0); //standard
  }  
  else {
    standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3., 30000.0); //standard
    standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3., 30000.0); //standard
  }
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  standardCut->SetRequestSPDany();
  processor->AddTrackCut(standardCut);
  processor->AddCandidateLeg1Cut(standardCut);

  AliReducedTrackCut* tightCut = new AliReducedTrackCut("tight","");
  tightCut->AddCut(AliReducedVarManager::kPt, 1.,100.);
  tightCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  tightCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  tightCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  tightCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  tightCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    tightCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    tightCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, 3.5, 30000.0);
    tightCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, 3.5, 30000.0);
  }  
  else {
    tightCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    tightCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 30000.0);
    tightCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 30000.0);
  }  
  tightCut->SetRejectKinks();
  tightCut->SetRequestITSrefit();
  tightCut->SetRequestTPCrefit();
  tightCut->SetRequestSPDany();
  processor->AddTrackCut(tightCut); 
  processor->AddCandidateLeg1Cut(tightCut);

  AliReducedTrackCut* looseCut = new AliReducedTrackCut("loose","");
  looseCut->AddCut(AliReducedVarManager::kPt, 1.,100.);
  looseCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  looseCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  looseCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  looseCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  looseCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    looseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    looseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, 2.5, 30000.0); //standard
    looseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, 2.5, 30000.0); //standard
  }
  else{
    looseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    looseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 2.5, 30000.0); //standard
    looseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 2.5, 30000.0); //standard
  }
  looseCut->SetRejectKinks();
  looseCut->SetRequestITSrefit();
  looseCut->SetRequestTPCrefit();
  looseCut->SetRequestSPDany();
  processor->AddTrackCut(looseCut); 
  processor->AddCandidateLeg1Cut(looseCut);

  // we also try with very loose cut in the TPC if we have a tof or trd signal

  // default at 3, proton at 2 if tof, pion at 2 if trd
  TF1* fCutLowNsigProton2 = new TF1("fCutLowNsigProton2", "(x <= -50 || x >= 100) * 3. + (x > -50 && x < 100) * 2", -1e3, 1e3);
  TF1* fCutHighNsigProton = new TF1("fCutHighNsigProton", "1000", -1e3, 1e3);

  TF1* fCutLowNsigPion2 = new TF1("fCutLowNsigPion2", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. + ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 2", -1e3, 1e3);
  TF1* fCutHighNsigPion = new TF1("fCutHighNsigPion", "1000", -1e3, 1e3);


  AliReducedTrackCut* looseCuttoftrd = new AliReducedTrackCut("looseProtonPion","");
  looseCuttoftrd->AddCut(AliReducedVarManager::kPt, 1., 100.);
  looseCuttoftrd->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  looseCuttoftrd->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  looseCuttoftrd->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  looseCuttoftrd->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  looseCuttoftrd->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, fCutLowNsigProton2, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, fCutLowNsigPion2, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }  
  else {
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, fCutLowNsigProton2, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    looseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, fCutLowNsigPion2, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }
  looseCuttoftrd->SetRejectKinks();
  looseCuttoftrd->SetRequestITSrefit();
  looseCuttoftrd->SetRequestTPCrefit();
  looseCuttoftrd->SetRequestSPDany();
  processor->AddTrackCut(looseCuttoftrd);
  processor->AddCandidateLeg1Cut(looseCuttoftrd);

  /*AliReducedTrackCut* tightCutMomDep = new AliReducedTrackCut("tightMomDep","");
  tightCutMomDep->AddCut(AliReducedVarManager::kPt, 1.,100.);
  tightCutMomDep->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  tightCutMomDep->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  tightCutMomDep->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  tightCutMomDep->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  tightCutMomDep->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);

  TF1* cutPIDTightMin = new TF1("cutPIDTightMin", "3.5 * (x < 5) + 2.5 * (x >= 5)", 0., 1000.);
  TF1* cutPIDTightMax = new TF1("cutPIDTightMax", "30000", 0., 1000.);
  if (correctTPCNsig) {
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, cutPIDTightMin, cutPIDTightMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, cutPIDTightMin, cutPIDTightMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
  }
  else {
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, cutPIDTightMin, cutPIDTightMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
    tightCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, cutPIDTightMin, cutPIDTightMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
  }
  tightCutMomDep->SetRejectKinks();
  tightCutMomDep->SetRequestITSrefit();
  tightCutMomDep->SetRequestTPCrefit();
  tightCutMomDep->SetRequestSPDany();
  processor->AddTrackCut(tightCutMomDep); 
  processor->AddCandidateLeg1Cut(tightCutMomDep);

  AliReducedTrackCut* standardCutMomDep = new AliReducedTrackCut("standardMomDep","");
  standardCutMomDep->AddCut(AliReducedVarManager::kPt, 1.,100.);
  standardCutMomDep->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCutMomDep->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCutMomDep->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCutMomDep->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCutMomDep->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);

  TF1* cutPIDMin = new TF1("cutPIDMin", "3. * (x < 5) + 2 * (x >= 5)", 0., 1000.);
  TF1* cutPIDMax = new TF1("cutPIDMax", "30000", 0., 1000.);

  if (correctTPCNsig) {
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, cutPIDMin, cutPIDMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, cutPIDMin, cutPIDMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
  }
  else {
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, cutPIDMin, cutPIDMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
    standardCutMomDep->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, cutPIDMin, cutPIDMax, kFALSE, AliReducedVarManager::kPin, 0., 1000.); //standard
  }
  standardCutMomDep->SetRejectKinks();
  standardCutMomDep->SetRequestITSrefit();
  standardCutMomDep->SetRequestTPCrefit();
  standardCutMomDep->SetRequestSPDany();
  processor->AddTrackCut(standardCutMomDep); 
  processor->AddCandidateLeg1Cut(standardCutMomDep);*/

  /*AliReducedTrackCut* noRejCut = new AliReducedTrackCut("noRejCut","");
  noRejCut->AddCut(AliReducedVarManager::kPt, 1.,100.);
  noRejCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  noRejCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  noRejCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  noRejCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  noRejCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  noRejCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); //standard
  noRejCut->SetRejectKinks();
  noRejCut->SetRequestITSrefit();
  noRejCut->SetRequestTPCrefit();
  noRejCut->SetRequestSPDany();
  processor->AddTrackCut(noRejCut); 
  processor->AddCandidateLeg1Cut(noRejCut);*/


  AliReducedTrackCut* veryLooseCut = new AliReducedTrackCut("veryLoose","");
  veryLooseCut->AddCut(AliReducedVarManager::kPt, 1., 100.);
  veryLooseCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  veryLooseCut->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  veryLooseCut->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  veryLooseCut->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  veryLooseCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, 2., 30000.0);
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, 2., 30000.0);
  }  
  else {
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 2., 30000.0);
    veryLooseCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 2., 30000.0);
  }
  veryLooseCut->SetRejectKinks();
  veryLooseCut->SetRequestITSrefit();
  veryLooseCut->SetRequestTPCrefit();
  veryLooseCut->SetRequestSPDany();
  //processor->AddTrackCut(veryLooseCut);
  //processor->AddCandidateLeg1Cut(veryLooseCut);



  // we also try with very loose cut in the TPC if we have a tof or trd signal

  // Default at 3, pion cut at 2 if trd, no proton cut if tof
  TF1* fCutLowNsigProton = new TF1("fCutLowNsigProton", "(x <= -10 || x >= 10) * 3. - (x > -10 && x < 10) * 1000", -1e3, 1e3);
  TF1* fCutLowNsigPion = new TF1("fCutLowNsigPion", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. + ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 2", -1e3, 1e3);

  // Default at 3, pion cut at 2 if trd, no proton cut if tof, proton cut at 2 if trd
  TF1* fCutLowNsigProton3 = new TF1("fCutLowNsigProton3", "(x <= -10 || x >= 10) * 3. - (x > -10 && x < 10) * 1000", -1e3, 1e3);
  TF1* fCutLowNsigProton3_TRD = new TF1("fCutLowNsigProton3_TRD", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. + ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 2", -1e3, 1e3);
  TF1* fCutLowNsigPion3 = new TF1("fCutLowNsigPion3", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. + ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 2", -1e3, 1e3);
 
  // Default at 3, no pion cut if trd, no proton cut if tof, proton cut at 2 if trd
  TF1* fCutLowNsigProton4_TRD = new TF1("fCutLowNsigProton4_TRD", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. + ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 2", -1e3, 1e3);
  TF1* fCutLowNsigPion4 = new TF1("fCutLowNsigPion4", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. - ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 1000", -1e3, 1e3);

  // Default at 3, no pion cut if trd, no proton cut if tof
  TF1* fCutLowNsigProton5 = new TF1("fCutLowNsigProton5", "(x <= -10 || x >= 10) * 3. - (x > -10 && x < 10) * 1000", -1e3, 1e3);
  TF1* fCutLowNsigPion5 = new TF1("fCutLowNsigPion5", "(x < 0. || (x >= 0.199999 && x <= 0.2000001)) * 3. - ((x >= 0. && x <= 0.199999) || (x >= 0.2000001 && x <= 1.)) * 1000", -1e3, 1e3);


  AliReducedTrackCut* veryLooseCuttoftrd = new AliReducedTrackCut("veryLooseProtonPion","");
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kPt, 1., 100.);
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, fCutLowNsigProton, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, fCutLowNsigPion, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }  
  else {
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, fCutLowNsigProton, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    veryLooseCuttoftrd->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, fCutLowNsigPion, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }
  veryLooseCuttoftrd->SetRejectKinks();
  veryLooseCuttoftrd->SetRequestITSrefit();
  veryLooseCuttoftrd->SetRequestTPCrefit();
  veryLooseCuttoftrd->SetRequestSPDany();
  //processor->AddTrackCut(veryLooseCuttoftrd);
  //processor->AddCandidateLeg1Cut(veryLooseCuttoftrd);

  AliReducedTrackCut* veryLooseCuttoftrd3 = new AliReducedTrackCut("veryLooseProtonPion3","");
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kPt, 1., 100.);
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, fCutLowNsigProton3_TRD, fCutHighNsigProton, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -10, 10, kTRUE);
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, fCutLowNsigPion3, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }  
  else {
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, fCutLowNsigProton3_TRD, fCutHighNsigProton, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -10, 10, kTRUE);
    veryLooseCuttoftrd3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, fCutLowNsigPion3, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }
  veryLooseCuttoftrd3->SetRejectKinks();
  veryLooseCuttoftrd3->SetRequestITSrefit();
  veryLooseCuttoftrd3->SetRequestTPCrefit();
  veryLooseCuttoftrd3->SetRequestSPDany();
  //processor->AddTrackCut(veryLooseCuttoftrd3);
  //processor->AddCandidateLeg1Cut(veryLooseCuttoftrd3);

  AliReducedTrackCut* veryLooseCuttoftrd4 = new AliReducedTrackCut("veryLooseProtonPion4","");
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kPt, 1., 100.);
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, fCutLowNsigProton4_TRD, fCutHighNsigProton, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -10, 10, kTRUE);
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, fCutLowNsigPion4, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }  
  else {
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, fCutLowNsigProton4_TRD, fCutHighNsigProton, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -10, 10, kTRUE);
    veryLooseCuttoftrd4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, fCutLowNsigPion4, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }
  veryLooseCuttoftrd4->SetRejectKinks();
  veryLooseCuttoftrd4->SetRequestITSrefit();
  veryLooseCuttoftrd4->SetRequestTPCrefit();
  veryLooseCuttoftrd4->SetRequestSPDany();
  //processor->AddTrackCut(veryLooseCuttoftrd4);
  //processor->AddCandidateLeg1Cut(veryLooseCuttoftrd4);

    AliReducedTrackCut* veryLooseCuttoftrd5 = new AliReducedTrackCut("veryLooseProtonPion5","");
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kPt, 1., 100.);
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kDcaXY, -1.0, 1.0);
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kDcaZ, -3.0, 3.0);
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCncls, 70., 160.0);
  veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  if (correctTPCNsig) {
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kProton, fCutLowNsigProton5, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kPion, fCutLowNsigPion5, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }  
  else {
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, fCutLowNsigProton5, fCutHighNsigProton, false, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron, -1e4, 1e4);
    veryLooseCuttoftrd5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, fCutLowNsigPion5, fCutHighNsigPion, false, AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron, -1e4, 1e4);
  }
  veryLooseCuttoftrd5->SetRejectKinks();
  veryLooseCuttoftrd5->SetRequestITSrefit();
  veryLooseCuttoftrd5->SetRequestTPCrefit();
  veryLooseCuttoftrd5->SetRequestSPDany();
  //processor->AddTrackCut(veryLooseCuttoftrd5);
  //processor->AddCandidateLeg1Cut(veryLooseCuttoftrd5);



  
  
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt07","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kPt, 0.2,100.0);
  //prefTrackCut1->SetRequestTPCrefit();
  /*if(ispp && isMC)*/ prefTrackCut1->SetTrackFilterBit(2); //electron prefilter
  //if(ispp && !isMC) prefTrackCut1->SetTrackFilterBit(1); //electron prefilter
  //if (!isMC) {
  //  prefTrackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.5, 3.5);
  //}
  //prefTrackCut1->SetRequestITSrefit();
  //prefTrackCut1->SetRejectPureMC(kTRUE);
  processor->AddCandidateLeg1PrefilterCut(prefTrackCut1);  
  
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddCandidateLeg1PairPrefilterCut(prefPairCut);

  
  AliReducedTrackCut* pairCut = new AliReducedTrackCut("JpsiCut","Pt pair selection");
  pairCut->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  pairCut->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  pairCut->AddCut(AliReducedVarManager::kMass, 1.5, 4.5);
  processor->AddCandidatePairCut(pairCut);  
  
  AliReducedTrackCut* pairCutDalitzPosField = new AliReducedTrackCut("DalitzCut","Pt pair selection");
  pairCutDalitzPosField->AddCut(AliReducedVarManager::kMass, 0.0, 0.015, false, AliReducedVarManager::kPt, 0., 1.);
  pairCutDalitzPosField->AddCut(AliReducedVarManager::kMass, 0.0, 0.035, false, AliReducedVarManager::kPt, 0., 1., true);
  TF1* fcutHigh = new TF1("f1", "[0] - [0]/[1]*x", -1.5, 1.5);
  fcutHigh->SetParameters(0.6, 0.12);
  TF1* fcutLow = new TF1("f1", "-[0] + [0]/[1]*x", -1.5, 1.5);
  fcutLow->SetParameters(0.6, 0.12);
  pairCutDalitzPosField->AddCut(AliReducedVarManager::kPairPsi, fcutLow, fcutHigh, true, AliReducedVarManager::kPairDeltaPhi, 0., 0.12);
  //processor->AddCandidatePairCut(pairCutDalitzPosField);
  // There should not be prefilter with this cut - and activate DCA histograms with many bins 

  AliReducedTrackCut* pairCutDalitzNegField = new AliReducedTrackCut("DalitzCut","Pt pair selection");
  pairCutDalitzNegField->AddCut(AliReducedVarManager::kMass, 0.0, 0.015, false, AliReducedVarManager::kPt, 0., 1.);
  pairCutDalitzNegField->AddCut(AliReducedVarManager::kMass, 0.0, 0.035, false, AliReducedVarManager::kPt, 0., 1., true);
  TF1* fcutHigh2 = new TF1("f1", "[0] + [0]/[1]*x", -1.5, 1.5);
  fcutHigh2->SetParameters(0.6, 0.12);
  TF1* fcutLow2 = new TF1("f1", "-[0] - [0]/[1]*x", -1.5, 1.5);
  fcutLow2->SetParameters(0.6, 0.12);
  pairCutDalitzNegField->AddCut(AliReducedVarManager::kPairPsi, fcutLow2, fcutHigh2, true, AliReducedVarManager::kPairDeltaPhi, -0.12, 0.);
  //processor->AddCandidatePairCut(pairCutDalitzNegField);  

  AliReducedTrackCut* jpsiMCCut = new AliReducedTrackCut("jpsiCutMC","Rapidity JpsiSelection");
  jpsiMCCut->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  jpsiMCCut->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC, 442.5, 443.5); 
  jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC+1, 442.5, 443.5, kTRUE);
  jpsiMCCut->SetMCFilterBit(0);
  jpsiMCCut->SetRejectPureMC(kFALSE);

  AliReducedTrackCut* jpsiMCCutPrompt = new AliReducedTrackCut("jpsiCutMCPrompt","Rapidity JpsiSelection");
  jpsiMCCutPrompt->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  jpsiMCCutPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  jpsiMCCutPrompt->AddCut(AliReducedVarManager::kPdgMC, 442.5, 443.5); 
  jpsiMCCutPrompt->AddCut(AliReducedVarManager::kPdgMC+1, 442.5, 443.5, kTRUE);
  jpsiMCCutPrompt->SetMCFilterBit(2);
  jpsiMCCutPrompt->SetRejectPureMC(kFALSE);

  AliReducedTrackCut* jpsiMCCutNonPrompt = new AliReducedTrackCut("jpsiCutMCNonPrompt","Rapidity JpsiSelection");
  jpsiMCCutNonPrompt->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  jpsiMCCutNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  jpsiMCCutNonPrompt->AddCut(AliReducedVarManager::kPdgMC, 442.5, 443.5); 
  jpsiMCCutNonPrompt->AddCut(AliReducedVarManager::kPdgMC+1, 442.5, 443.5, kTRUE);
  jpsiMCCutNonPrompt->SetMCFilterBit(1);
  jpsiMCCutNonPrompt->SetRejectPureMC(kFALSE);


  AliReducedTrackCut* electronMCCut = new AliReducedTrackCut("standardMC","Pt electron Selection");
  processor->AddJpsiMotherMCCut(jpsiMCCut, electronMCCut);
  processor->AddJpsiMotherMCCut(jpsiMCCutPrompt, electronMCCut);
  processor->AddJpsiMotherMCCut(jpsiMCCutNonPrompt, electronMCCut);
  if (isInjected) processor->SetReweightCut(1); // ML - Only Reweighting on prompt Jpsi
  //if (isInjected) processor->SetReweightCut(0); // Standard - Reweighting on inclusive Jpsi

  AliReducedTrackCut* jpsiMCCutEta = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEta");
  AliReducedTrackCut* electronMCCutEta = new AliReducedTrackCut("standardMCEta","Eta electron Selection");
  electronMCCutEta->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEta, electronMCCutEta);


  AliReducedTrackCut* jpsiMCCutEtaPt = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEtaPt");
  AliReducedTrackCut* electronMCCutEtaPt = new AliReducedTrackCut("standardMCEtaPt","Eta and Pt electron Selection");
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kPt, 1., 1e5);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEtaPt, electronMCCutEtaPt);

  TFile* fileJpsiMass = new TFile("/gluster1/glegras/JpsiMassDist.root","read");
  TF1* fJpsiMass; fileJpsiMass->GetObject("signalFitCB", fJpsiMass);
  if(!fJpsiMass) {
    fileJpsiMass = new TFile("~/alice/AliPhysics/PWGDQ/reducedTree/macros/JpsiMassDist.root","read");
    fileJpsiMass->GetObject("signalFitCB", fJpsiMass);
  }
  processor->SetJpsiMassDist(fJpsiMass);
  if(!processor->GetJpsiMassDist()) cout<<"No jpsi mass distribution setup"<<endl;

  if (processor->GetReweightParticleComposition() && processor->GetRunOverMC()) {

    std::cout<<"Reweighting of particle composition is activated"<<std::endl;
    bool weightsPythia = kTRUE; // alternative: EPOS
    bool weightsWrtData = kTRUE; // alternative: weight with respect to EPOS
    bool corrPtDist = kFALSE;
    bool corrMult = kTRUE;

    TFile* fDataWeights = new TFile("/gluster1/glegras/AllPublishedFractions.root", "read");
    if (!fDataWeights->IsOpen() || fDataWeights->GetSize() < 1000.) fDataWeights = new TFile("~/alice/data/PCC/AllPublishedFractions.root", "read"); // EPOS
    
    // Choose either EPOS or Pythia
    TFile* fMCWeights;
    if (weightsPythia) {
      fMCWeights = new TFile("/gluster1/glegras/FirstTrain_pp_Pythia8SO.root", "read"); // Pythia
      if (!fMCWeights->IsOpen() || fMCWeights->GetSize() < 1000.) fMCWeights = new TFile("~/alice/data/PCC/FirstTrain_pp_Pythia8SO.root", "read"); // Pythia
    }
    else {
      fMCWeights = new TFile("/gluster1/glegras/FirstTrain_pp_LHC17d20ab12_EPOSLHC.root", "read"); // EPOS
      if (!fMCWeights->IsOpen() || fMCWeights->GetSize() < 1000.) fMCWeights = new TFile("~/alice/data/PCC/FirstTrain_pp_LHC17d20ab12_EPOSLHC.root", "read"); // EPOS
    }

    TDirectoryFile* dir1; TList* dir2;
    fMCWeights->GetObject("AliMCWeightsTask", dir1);
    dir1->GetObject("AliMCWeightsTask", dir2);
    TH3F* hMCWeights = (TH3F*) dir2->FindObject("fHistMCGenPrimTrackParticle");

    // For MC closure, if we want to unfold EPOS with detector response matrix from Pythia, the weights to apply to Pythia need to be EPOS/Pythia
    TFile* fileEPOS = new TFile("/gluster1/glegras/FirstTrain_pp_LHC17d20ab12_EPOSLHC.root", "read"); // EPOS
    if (!fileEPOS->IsOpen() || fileEPOS->GetSize() < 1000.) fileEPOS = new TFile("~/alice/data/PCC/FirstTrain_pp_LHC17d20ab12_EPOSLHC.root", "read"); // EPOS

    fileEPOS->GetObject("AliMCWeightsTask", dir1);
    dir1->GetObject("AliMCWeightsTask", dir2);
    TH3F* hMCWeightsEPOS = (TH3F*) dir2->FindObject("fHistMCGenPrimTrackParticle");


    const char* dataSyst = "Bylinkin"; // Function used for pt distribution extrapolation
    const char* partNames[5] = {"Pion", "Proton", "Kaon", "SigmaPlus", "SigmaMinus"};
    Double_t binsMult[11] = {0., 3.5, 5.4, 7.1, 8.99, 11.03, 12.91, 14.97, 18.09, 23.02, 49.};
    Double_t binsCent[11] = {0., 0.92, 4.6, 9.2, 13.8, 18.4, 27.6, 36.8, 46, 64.5, 100};
    TH3F* fWeightsPCC = new TH3F("fWeightsPCC", "Weights for Particle Composition Correction",              
              hMCWeights->GetNbinsX(),  hMCWeights->GetXaxis()->GetXbins()->GetArray(),
              10, binsCent,
              hMCWeights->GetNbinsZ(), hMCWeights->GetZaxis()->GetXbins()->GetArray());

    // Compute fractions MC identified from absolute abundancies (PCC = data/MC)
    for (int imult = 1; imult <= fWeightsPCC->GetNbinsY(); imult++) {
      for (int ipt = 1; ipt <= fWeightsPCC->GetNbinsX(); ipt++) {
        float totalIdentifiedMC = hMCWeights->Integral(ipt, ipt, 
                  hMCWeights->GetYaxis()->FindBin(binsMult[imult-1]) + 1, hMCWeights->GetYaxis()->FindBin(binsMult[imult]), 1, 5);
        for (int itype = 1; itype <= 5; itype++) {
          if (ipt > 3 && ipt < 51) {// Otherwise we might have 0, inf, nan
            float fractionMC = hMCWeights->Integral(ipt, ipt, hMCWeights->GetYaxis()->FindBin(binsMult[imult-1]) + 1, 
                                      hMCWeights->GetYaxis()->FindBin(binsMult[imult]), itype, itype) / totalIdentifiedMC;
            fWeightsPCC->SetBinContent(ipt, 11 - imult, itype, 1./fractionMC);
          }
        }
      }
    }

    // Fractions data identified
    if (weightsWrtData) {
      for (int imult = 1; imult <= fWeightsPCC->GetNbinsY(); imult++) {
        for (int itype = 1; itype <= 5; itype++) {
          const char* histName = Form("pp%s%d%s", partNames[itype-1], imult - 1, dataSyst);
          TH1D* hDataWeights; fDataWeights->GetObject(histName, hDataWeights);     
          for (int ipt = 1; ipt <= fWeightsPCC->GetNbinsX(); ipt++) {
            float fractionData = hDataWeights->GetBinContent(ipt);
            fWeightsPCC->SetBinContent(ipt, imult, itype, fractionData * fWeightsPCC->GetBinContent(ipt, imult, itype)); 
            // Protecting from 0, inf, nan
            if (ipt <= 3) fWeightsPCC->SetBinContent(ipt, imult, itype, hDataWeights->GetBinContent(4) * fWeightsPCC->GetBinContent(4, imult, itype));
            if (ipt >= 51) fWeightsPCC->SetBinContent(ipt, imult, itype, fWeightsPCC->GetBinContent(50, imult, itype)); 
          }
        }
      }
    }
    else { // Rather use this one as fractions from data if you are doing MC closure test with EPOS treated as data
      for (int imult = 1; imult <= fWeightsPCC->GetNbinsY(); imult++) {
        for (int ipt = 4; ipt <= fWeightsPCC->GetNbinsX(); ipt++) {
          float totalIdentifiedMC = hMCWeightsEPOS->Integral(ipt, ipt, 
                    hMCWeightsEPOS->GetYaxis()->FindBin(binsMult[imult-1]) + 1, hMCWeightsEPOS->GetYaxis()->FindBin(binsMult[imult]), 1, 5);
          for (int itype = 1; itype <= 5; itype++) {
            float fractionEPOS = hMCWeightsEPOS->Integral(ipt, ipt, hMCWeightsEPOS->GetYaxis()->FindBin(binsMult[imult-1]) + 1, 
                                        hMCWeightsEPOS->GetYaxis()->FindBin(binsMult[imult]), itype, itype) / totalIdentifiedMC;
            fWeightsPCC->SetBinContent(ipt, 11 - imult, itype, fractionEPOS * fWeightsPCC->GetBinContent(ipt, 11 - imult, itype)); 
            // Protecting from 0, inf, nan
            if (ipt >= 51) fWeightsPCC->SetBinContent(ipt, 11 - imult, itype, fWeightsPCC->GetBinContent(50, 11 - imult, itype)); 
          }
        }
        for (int ipt = 1; ipt <= 3; ipt++) {
          for (int itype = 1; itype <= 5; itype++) {
            fWeightsPCC->SetBinContent(ipt, 11 - imult, itype, fWeightsPCC->GetBinContent(4, 11 - imult, itype));
          }
        }
      }
    }


    // Correct also absolute pT distributions, we use pions for that
    if (corrPtDist) {
      TH2F* hMCPtDist; fMCWeights->GetObject("Mult_Pt", hMCPtDist);
      for (int imult = 1; imult <= 10; imult++) {
        TH1F* hWeightsVsPt = new TH1F(Form("hWeightsVsPt%d", imult-1), Form("hWeightsVsPt%d", imult-1), 
                                                  hMCWeights->GetNbinsX(),  hMCWeights->GetXaxis()->GetXbins()->GetArray());
        const char* histName = Form("PtDist_pion_%d", imult - 1);
        TH1F* hDataWeights; 
        if (weightsWrtData) fDataWeights->GetObject(histName, hDataWeights);
        else {// Rather use this one as pt dist from data if you are doing MC closure test with EPOS treated as data
          hDataWeights = (TH1F*) ((TH2F*) fileEPOS->FindObjectAny("Mult_Pt"))->ProjectionX(Form("proj%d", imult), imult, imult);
          for (int i = 1; i <= hDataWeights->GetNbinsX(); i++) 
            hDataWeights->SetBinContent(i, hDataWeights->GetBinContent(i) / (hDataWeights->GetXaxis()->GetBinUpEdge(i) - hDataWeights->GetXaxis()->GetBinLowEdge(i)));
        }

        for (int ipt = 1; ipt <= hMCWeights->GetNbinsX(); ipt++) {
          hWeightsVsPt->SetBinContent(ipt, hDataWeights->GetBinContent(hDataWeights->FindBin(hMCWeights->GetXaxis()->GetBinCenter(ipt)))
                                                                          / hMCPtDist->GetBinContent(ipt, imult) 
                                                                          * (hMCPtDist->GetXaxis()->GetBinUpEdge(ipt) - hMCPtDist->GetXaxis()->GetBinLowEdge(ipt)));
          if (ipt <= 3 && weightsWrtData) hWeightsVsPt->SetBinContent(ipt, hDataWeights->GetBinContent(hDataWeights->FindBin(hMCWeights->GetXaxis()->GetBinCenter(4)))
                                                                          / hMCPtDist->GetBinContent(4, imult)
                                                                          * (hMCPtDist->GetXaxis()->GetBinUpEdge(4) - hMCPtDist->GetXaxis()->GetBinLowEdge(4)));
          if (ipt >= 51) hWeightsVsPt->SetBinContent(ipt, hDataWeights->GetBinContent(hDataWeights->FindBin(hMCWeights->GetXaxis()->GetBinCenter(50)))
                                                                          / hMCPtDist->GetBinContent(50, imult)
                                                                          * (hMCPtDist->GetXaxis()->GetBinUpEdge(50) - hMCPtDist->GetXaxis()->GetBinLowEdge(50)));
        }
        // We want sum(weights*scale*counts) = sum(counts) -> scale = sum(counts)/sum(weights*counts)
        float sumCounts = hMCPtDist->Integral(weightsWrtData ? 4 : 1, 50, imult, imult);
        TH1F* hMCWeightedCounts = (TH1F*) hMCPtDist->ProjectionX(Form("weightedcounts%d", imult), imult, imult);
        hMCWeightedCounts->Multiply(hWeightsVsPt);
        float sumWeightedCounts = hMCWeightedCounts->Integral(weightsWrtData ? 4 : 1, 50);
        hWeightsVsPt->Scale(sumCounts / sumWeightedCounts);

        for (int itype = 1; itype <= 5; itype++) {
          for (int ipt = 1; ipt <= hMCWeights->GetNbinsX(); ipt++) {
            fWeightsPCC->SetBinContent(ipt, imult, itype, fWeightsPCC->GetBinContent(ipt, imult, itype) * hWeightsVsPt->GetBinContent(ipt));
          }
        }
      }
    }

    if (!corrMult) {// if we want to remove multiplicity dependence
      for (int imult = 1; imult <= fWeightsPCC->GetNbinsY(); imult++) {
        for (int itype = 1; itype <= 5; itype++) {
          for (int ipt = 1; ipt <= fWeightsPCC->GetNbinsX(); ipt++) {
            fWeightsPCC->SetBinContent(ipt, imult, itype, fWeightsPCC->GetBinContent(ipt, 4, itype));
          }
        }
      }
    }
    processor->SetMCParticleCompositionWeights(fWeightsPCC);
  }

  if (correctTPCNsig) {
    AliReducedVarManager::Variables vars[4] = {AliReducedVarManager::kPin, AliReducedVarManager::kEta, AliReducedVarManager::kRunNo, AliReducedVarManager::kNothing};
    AliReducedVarManager::SetTPCpidCalibDepVars(vars);

    THnF* centroid_el; THnF* width_el; THnI* status_el;
    THnF* centroid_pi; THnF* width_pi; THnI* status_pi;
    THnF* centroid_pr; THnF* width_pr; THnI* status_pr;

    TFile* file = new TFile("/gluster1/glegras/PostCalib.root");
    file->GetObject("centroid_el", centroid_el);
    if (!centroid_el) {
      file = new TFile("/home/glegras/alice/AliPhysics/PWGDQ/reducedTree/macros/PostCalib.root");
      file->GetObject("centroid_el", centroid_el);
    }
    file->GetObject("centroid_pi", centroid_pi);
    file->GetObject("centroid_pr", centroid_pr);
    file->GetObject("width_el", width_el);
    file->GetObject("width_pi", width_pi);
    file->GetObject("width_pr", width_pr);
    file->GetObject("status_el", status_el);
    file->GetObject("status_pi", status_pi);
    file->GetObject("status_pr", status_pr);

    AliReducedVarManager::SetTPCpidCalibMaps(0, centroid_el, width_el, status_el);
    AliReducedVarManager::SetTPCpidCalibMaps(1, centroid_pi, width_pi, status_pi);
    AliReducedVarManager::SetTPCpidCalibMaps(2, centroid_pr, width_pr, status_pr);
  }

  
  SetupMixingHandler(processor);
  SetupHistogramManager(processor,ispp, prod);
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, ispp, prod);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kRunID);
}

//_______________________________________________________________________
void SetupMixingHandler(AliReducedAnalysisFilterTrees* task) {
   //
   // setup the mixing handler
   AliMixingHandler* handler = task->GetMixingHandler();  
   handler->SetPoolDepth(100);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1.0);
   handler->SetDownscaleTracks(1);
   handler->SetMixLikeSign(kTRUE);
   const Int_t nZbinLimits = 2;
   Float_t zLims[nZbinLimits] = {-10., /*  -8.,  -6.,  -4.,  -2., 0.,  2.,  4.,  6.,  8., */10. };
   //r(Int_t i=0;i<=10;++i) zLims[nZbinLimits] = i
   const int nMultBins = 8; Float_t multBins[nMultBins+1] = {0,10,20,30,40,50,60,70,150};

  //handler->AddMixingVariable(AliReducedVarManager::kVtxZ, nZbinLimits, zLims);
  handler->AddMixingVariable(AliReducedVarManager::kNGlobalTracks, nMultBins+1, multBins); 
  TString histClassNames = handler->GetHistClassNames();
  TObjArray* histClassArr = histClassNames.Tokenize(";");

  //if (!task->GetRunOverMC()) task->SetRunEventMixingMult(kTRUE);

  
  task->SetMultBinsMixing(nMultBins, multBins);
  for (int i = 0; i < nMultBins; i++){
    AliMixingHandler* handlerMult = (AliMixingHandler*) handler->Clone();
    handlerMult->SetPoolDepth(200);
    if (i == 0) handlerMult->SetPoolDepth(5000);
    TString histNamesMult = "";
    for (int j = 0; j<histClassArr->GetEntries(); j++) {
        histNamesMult += histClassArr->At(j)->GetName();
        histNamesMult += Form("_%d_%d;", (int) multBins[i],(int) multBins[i+1]-1);
    }
    handlerMult->SetHistClassNames(histNamesMult);   
    handlerMult->SetHistogramManager(task->GetHistogramManager());
    task->AddMixingHandler(handlerMult);
  }
}


//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts
  //           are used in the histogram lists
  // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
  //  TODO: make needed changes such that this becomes less prone to mistakes
   
  TString runNumbers = ""; const char* fileRuns; const char* fileRuns2;
  int nRuns = 0;
  if(ispp) {
    fileRuns = "/home/glegras/alice/data/runLists/list_all.txt";
    fileRuns2 = "/gluster1/glegras/list_all.txt";
  }
  else {
    fileRuns = "/home/glegras/alice/data/runLists/runspPb.txt";
    fileRuns2 = "/gluster1/glegras/runspPb.txt";
  }
  fstream newfile;
  newfile.open(fileRuns,ios::in);
  if(!(newfile.is_open())) newfile.open(fileRuns2,ios::in);
  if(newfile.is_open()) {
    string temp;
    while(getline(newfile, temp)){
      runNumbers+=temp;
      runNumbers+=";";
      nRuns++;
    }
    newfile.close(); 
  }

  AliReducedVarManager::SetRunNumbers(runNumbers);

  AliHistogramManager* man = task->GetHistogramManager(); 

  bool isMC = task->GetRunOverMC();
  
  TString histClasses = "";
  if(isMC) histClasses += "Event_BeforeCuts;";
  //if(isMC) histClasses += "Event_INELGT0;";
  //if(isMC) histClasses += "Event_INELGT0_MB;";
  if(!isMC) histClasses += "Event_NoVtxRec;";
  if(!isMC) histClasses += "Event_NoVtxzCut;";
  histClasses += "Event_AfterCuts;";
  //histClasses += "PionPrimary;";



  int nEstimators = (task->GetComputeMult() ? task->GetNMeasMultCuts() : 0);

  for (int icut = 0; icut<nEstimators; icut++) {
    //if(isMC) histClasses += Form("EventMult_%s_INELGT0;", task->GetMeasMultcutName(icut));
    //if(isMC) histClasses += Form("EventMult_%s_INELGT0_MB;", task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMult_%s_NoVtxRec;", task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMult_%s_NoVtxzCut;", task->GetMeasMultcutName(icut));
    
    if(!isMC) histClasses += Form("EventMult_%s_Inclusive;", task->GetMeasMultcutName(icut));
    histClasses += Form("EventMult_%s_MB;", task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMult_%s_HM;", task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMult_%s_HighMultSPD;", task->GetMeasMultcutName(icut));

    if(!isMC) histClasses += Form("EventMultRegions_%s_Inclusive;",task->GetMeasMultcutName(icut));
    histClasses += Form("EventMultRegions_%s_MB;",task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMultRegions_%s_HM;",task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMultRegions_%s_HighMultSPD;",task->GetMeasMultcutName(icut));

    if(!isMC) histClasses += Form("EventMultRegions2Leading_%s_Inclusive;",task->GetMeasMultcutName(icut));
    histClasses += Form("EventMultRegions2Leading_%s_MB;",task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMultRegions2Leading_%s_HM;",task->GetMeasMultcutName(icut));
    if(!isMC) histClasses += Form("EventMultRegions2Leading_%s_HighMultSPD;",task->GetMeasMultcutName(icut));
  }

  /*if (!isMC) histClasses += "CorrelMult_Inclusive;";
  histClasses += "CorrelMult_MB;";
  if (!isMC) histClasses += "CorrelMult_HM;";
*/
  
  //histClasses += "EventTag_BeforeCuts;";   
  //histClasses += "EventTag_AfterCuts;";   
  // histClasses += "EventTriggers_BeforeCuts;";
  //histClasses += "EventTriggers_AfterCuts;";   
  
  if(task->GetWriteFilteredTracks()) {
    //histClasses += "Track_BeforeCuts;";
    for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
    }
  }
  
  if(task->GetWriteFilteredPairs()) {
    //histClasses += "Pair_BeforeCuts;";
    // histClasses += "PairQualityFlags_BeforeCuts;";
    for(Int_t i=0; i<task->GetNPairCuts(); ++i) {
      TString cutName = task->GetPairCutName(i);
      histClasses += Form("Pair_%s;", cutName.Data());
      //histClasses += Form("PairQualityFlags_%s;", cutName.Data());
    }
  }
  
  if(task->GetBuildCandidatePairs()) {
    for(Int_t i=0; i<task->GetNCandidateLegCuts(); ++i) {
      //histClasses += Form("Track_LEG1_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,1));
      //histClasses += Form("Track_LEG2_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,2));
      if(task->GetRunCandidatePrefilter()) {
        histClasses += Form("Track_LEG1_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,1));
        histClasses += Form("Track_LEG2_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,2));
      }
      //histClasses += Form("Track_LEG1_AfterPairCuts_%s;", task->GetCandidateLegCutName(i,1));
      //histClasses += Form("Track_LEG2_AfterPairCuts_%s;", task->GetCandidateLegCutName(i,2));
      histClasses += Form("Pair_Candidate12_%s%s;", 
			  task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      if(task->GetBuildCandidateLikePairs()) {
        	histClasses += Form("Pair_Candidate11_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
	        histClasses += Form("Pair_Candidate22_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      }
      if(task->GetRunEventMixing()) {
          histClasses += Form("PairMEPM_%s%s;PairMEPP_%s%s;PairMEMM_%s%s;",
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));

      }
      if(task->GetRunEventMixingMult()) {
        int nBins = task->GetNMultBinsMixing();
        for (int j = 0; j<nBins; j++) {
          int lowMult = (int) task->GetMultBinsMixing(j); int highMult = (int) task->GetMultBinsMixing(j+1)-1;
          histClasses += Form("PairMEPM_%s%s_%d_%d;PairMEPP_%s%s_%d_%d;PairMEMM_%s%s_%d_%d;",
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult,highMult,
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult,highMult,
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult, highMult);
        }
      }
    }
    if(task->GetRunCandidatePrefilter()) {
        //histClasses += "Track_LEG1_PrefilterTrack;";
        //histClasses += "Track_LEG2_PrefilterTrack;";
	}
  }

  if(isMC) {
    for (int i=0; i<task->GetNJpsiMotherMCCuts();i++) {
      //histClasses += Form("PureMCTRUTH_BeforeSelection_%s;", task->GetJpsiMotherMCcutName(i));
      histClasses += Form("PureMCTRUTH_AfterSelection_%s;", task->GetJpsiMotherMCcutName(i));
      for (int jcut = 0; jcut<task->GetNMeasMultCuts(); jcut ++) {
        histClasses += Form("JpsiPtMultCorrel_%s_%s;", task->GetMeasMultcutName(jcut), task->GetJpsiMotherMCcutName(i));
        for (int icutleg = 0; icutleg < task->GetNCandidateLegCuts(); icutleg++) 
          histClasses += Form("JpsiPtMultCorrelMeas_%s_%s_%s;", task->GetMeasMultcutName(jcut), task->GetJpsiMotherMCcutName(i), task->GetCandidateLegCutName(icutleg, 1));
      }
      for (int icutleg = 0; icutleg < task->GetNCandidateLegCuts(); icutleg++) 
        histClasses += Form("PureMCTRUTH_DetectedDaughters_%s_%s;", task->GetJpsiMotherMCcutName(i), task->GetCandidateLegCutName(icutleg, 1));
    }
    histClasses += "DeltaPhi_JpsiTruth_JpsiCandidate;";
  } 
   
 
  for (int cutMode = 0; cutMode < 4*nEstimators; cutMode++) { 
    // 0 for MB multiplicity  unfolding
    // 1 for HM multiplicity unfolding
    // 2 for MB+HM multiplicity unfolding
    // 3 is for Jpsi vs multiplicity unfolding (HM+MB)
    //if (cutMode >= 8) continue;
    if(ispp){
      if(isMC) histClasses += Form("pp_13TeV_MC_cutMode_%d;",100+cutMode);
      else histClasses += Form("pp_13TeV_Data_cutMode_%d;",100+cutMode);
    }
    else {
      if(isMC) histClasses += Form("pPb_5TeV_MC_cutMode_%d;",100+cutMode);
      else histClasses += Form("pPb_5TeV_Data_cutMode_%d;",100+cutMode);
    }
  }
  
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};


  int nMultBins; float minMult; float maxMult;
  int nMultBinspp = 151; float minMultpp = -0.5; float maxMultpp = 150.5;
  int nMultBinspPb = 251; float minMultpPb = -0.5; float maxMultpPb = 250.5;

  if (ispp) {
    nMultBins = nMultBinspp; minMult = minMultpp; maxMult = maxMultpp;
  }
  else {
    nMultBins = nMultBinspPb; minMult = minMultpPb; maxMult = maxMultpPb;
  }
  
  // Pb-Pb from 2010 run range is default: LHC10h
  /* runNBins = 2500;
  runHistRange[0] = 137100.;
  runHistRange[1] = 139600.;*/
  if(ispp){
    runNBins = 42800;
    runHistRange[0] = 252200.;
    runHistRange[1] = 295000.;
  }
  else{
    runNBins = 300;
    runHistRange[0] = 265300.;
    runHistRange[1] = 265600.;
  }
  
  // Pb-Pb of 2011
  if(prod.Contains("LHC11h")) {
    runNBins = 2700;
    runHistRange[0] = 167900.;
    runHistRange[1] = 170600.;
  }
  
  // Pb-Pb of 2015
  if(prod.Contains("LHC15o")) {
    runNBins = 2100;
    runHistRange[0] = 244900.;
    runHistRange[1] = 247000.;
  }
  
  // p-Pb of 2013
  if(prod.Contains("LHC13b") || prod.Contains("LHC13c")) {
    runNBins = 400;
    runHistRange[0] = 195300.;
    runHistRange[1] = 195700.;
  }
  
  // pp at 13 TeV
  if(prod.Contains("LHC16l")) {
    runNBins = 1140;
    runHistRange[0] = 258880.;
    runHistRange[1] = 260020.;
  }
  
  // p-Pb at 8.16 TeV
  if(prod.Contains("LHC16r")) {
    runNBins = 1000;
    runHistRange[0] = 265400.;
    runHistRange[1] = 266400.;
  }
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();

    // Event wise histograms
    if(classStr.Contains("EventTag_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString tagNames = "";
      tagNames += "AliAnalysisUtils 2013 selection;";
      tagNames += "AliAnalysisUtils MV pileup;";
      tagNames += "AliAnalysisUtils MV pileup, no BC check;";
      tagNames += "AliAnalysisUtils MV pileup, min wght dist 10;";
      tagNames += "AliAnalysisUtils MV pileup, min wght dist 5;";
      tagNames += "IsPileupFromSPD(3,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(4,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(5,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(6,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(3,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(4,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(5,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(6,0.8,3.,2.,5.);";
      tagNames += "vtx distance selected;";
      man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
			20, -0.5, 19.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
      man->AddHistogram(classStr.Data(), "EventTags_CentVZERO", "Event tags vs VZERO centrality", kFALSE,
			20, -0.5, 19.5, AliReducedVarManager::kEventTag, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
      continue;
    }
    
    if(classStr.Contains("EventTriggers_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
       
      man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
			64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
			64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      continue;
    }
    
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      man->AddHistogram(classStr.Data(),"RunID","Run ID",kFALSE, nRuns, 0, nRuns, AliReducedVarManager::kRunID);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X", kFALSE, 1000, -0.5, 0.5, AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y", kFALSE, 1000, -0.5, 0.5, AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z", kFALSE, 1000, -15., 15., AliReducedVarManager::kVtxZ);
      if (isMC) {
        man->AddHistogram(classStr.Data(),"VtxXDeltaTrue","Vtx X reconstructed - true", kFALSE, 1000, -0.1, 0.1, AliReducedVarManager::kVtxXDeltaTrue);
        man->AddHistogram(classStr.Data(),"VtxYDeltaTrue","Vtx Y reconstructed - true", kFALSE, 1000, -0.1, 0.1, AliReducedVarManager::kVtxYDeltaTrue);
        man->AddHistogram(classStr.Data(),"VtxZDeltaTrue","Vtx Z reconstructed - true", kFALSE, 1000, -5., 5., AliReducedVarManager::kVtxZDeltaTrue);
      }
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,30000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxContributors","", kFALSE, 200, 0., 200., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"NVtxContributors_SPDNtracklets", "", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"NVtxContributors_V0M","", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, 200, 0., 200., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TPCpileupContributorsAC", "SPD #tracklets in |#eta|<1.0 vs TPCpileupContributorsAC", kFALSE, 800, 0, 800, AliReducedVarManager::kTPCpileupContributorsAC, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TPCpileupZAC", "SPD #tracklets in |#eta|<1.0 vs TPCpileupZAC", kFALSE, 400, -200, 200, AliReducedVarManager::kTPCpileupZAC, 200, -0., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_V0M", "SPD #tracklets in |#eta|<1.0 vs V0M", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_V0A", "SPD #tracklets in |#eta|<1.0 vs V0A", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROATotalMult, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_V0C", "SPD #tracklets in |#eta|<1.0 vs V0C", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROCTotalMult, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      continue;
    }  // end if className contains "Event" 
            
    for (int icut = 0; icut<task->GetNMeasMultCuts(); icut++) {
      const char* cutName = task->GetMeasMultcutName(icut);
      if(classStr.Contains(Form("EventMult_%s",cutName))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(),"NGlobalTracks","Number of global tracks per event", kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NVtxContributors_NGlobalTracks", "", kFALSE, 200, 0., 200., AliReducedVarManager::kNGlobalTracks + icut, 200, 0., 200., AliReducedVarManager::kNVtxContributors);
        man->AddHistogram(classStr.Data(),"Mean_NGlobalTracks_V0Mult","Mean Number of global tracks per event vs V0 Mult", kTRUE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"Mean_V0Mult_NGlobalTracks","Mean V0 Mult vs Number of global tracks per event", kTRUE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_V0M","Number of global tracks per event vs V0M", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_V0A","Number of global tracks per event vs V0A", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROATotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_V0C","Number of global tracks per event vs V0C", kFALSE, 800, 0, 800, AliReducedVarManager::kVZEROCTotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        if (isMC) {
          man->AddHistogram(classStr.Data(),"Multiplicity","Multiplicity in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kTRUE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
          man->AddHistogram(classStr.Data(),"NGlobalTracks_Multiplicity_TProfile","Mean global tracks in |#eta|<0.9 vs multiplicity (in |#eta|<0.9)",kTRUE, nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks","Multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09);
        }
        man->AddHistogram(classStr.Data(),"SPDNtracklets_NGlobalTracks","SPDNtracklets vs Global tracks (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kSPDntracklets);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_vtxz_TProfile","Mean number of global tracks (in |#eta|<0.9) vs vtxz",kTRUE,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_RunNumber_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run number",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_RunID_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run ID",kTRUE, nRuns, 0, nRuns, AliReducedVarManager::kRunID,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        continue;
      }
    }

    if(classStr.Contains("CorrelMult")) { //Correlations between different multiplicity estimators
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"NGlobalTracks2_NGlobalTracks","Number of global tracks with 2nd mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 1);
      man->AddHistogram(classStr.Data(),"NGlobalTracks3_NGlobalTracks","Number of global tracks with 3rd mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 2);
      man->AddHistogram(classStr.Data(),"NGlobalTracks4_NGlobalTracks","Number of global tracks with 4th mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 3);
      man->AddHistogram(classStr.Data(),"NGlobalTracks5_NGlobalTracks","Number of global tracks with 5th mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 4);
      man->AddHistogram(classStr.Data(),"NGlobalTracks6_NGlobalTracks","Number of global tracks with 6th mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 5);
      man->AddHistogram(classStr.Data(),"NGlobalTracks7_NGlobalTracks","Number of global tracks with 7th mult cut vs 1st mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 6);
      man->AddHistogram(classStr.Data(),"NGlobalTracks4_NGlobalTracks2","Number of global tracks with 4th mult cut vs 2nd mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 1, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 3);
      man->AddHistogram(classStr.Data(),"NGlobalTracks6_NGlobalTracks2","Number of global tracks with 6th mult cut vs 2nd mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 1, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 5);
      man->AddHistogram(classStr.Data(),"NGlobalTracks5_NGlobalTracks3","Number of global tracks with 5th mult cut vs 3rd mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 2, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 4);
      man->AddHistogram(classStr.Data(),"NGlobalTracks7_NGlobalTracks3","Number of global tracks with 7th mult cut vs 3rd mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 2, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 6);
      man->AddHistogram(classStr.Data(),"NGlobalTracks3_NGlobalTracks2","Number of global tracks with 3rd mult cut vs 2nd mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 1, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 2);
      man->AddHistogram(classStr.Data(),"NGlobalTracks5_NGlobalTracks4","Number of global tracks with 5th mult cut vs 4th mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 3, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 4);
      man->AddHistogram(classStr.Data(),"NGlobalTracks7_NGlobalTracks6","Number of global tracks with 7th mult cut vs 6th mult cut",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 5, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + 6);

      continue;
    }

    for (int icut = 0; icut<task->GetNMeasMultCuts(); icut++) {
      const char* cutName = task->GetMeasMultcutName(icut);
      if(classStr.Contains(Form("EventMultRegions_%s",cutName))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(),"NGlobalTracksToward","Number of global tracks (Toward) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse","Number of global tracks (Transverse) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracksAway","Number of global tracks (Away) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksToward","Number of global tracks (all regions) vs Number of global tracks (Toward) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksTransverse","Number of global tracks (all regions) vs Number of global tracks (Transverse) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksAway","Number of global tracks (all regions) vs Number of global tracks (Away) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        if (isMC) {
          man->AddHistogram(classStr.Data(),"MultiplicityToward","Multiplicity (toward) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Toward_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (toward) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Toward","Multiplicity in |#eta|<0.9 vs global tracks (toward) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
          man->AddHistogram(classStr.Data(),"MultiplicityTransverse","Multiplicity (Transverse) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Transverse_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (Transverse) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Transverse","Multiplicity in |#eta|<0.9 vs global tracks (Transverse) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
          man->AddHistogram(classStr.Data(),"MultiplicityAway","Multiplicity (Away) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Away_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (Away) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_Away","Multiplicity in |#eta|<0.9 vs global tracks (Away) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
        }
        continue;
      }
    }

    for (int icut = 0; icut<task->GetNMeasMultCuts(); icut++) {
      const char* cutName = task->GetMeasMultcutName(icut);
      if(classStr.Contains(Form("EventMultRegions2Leading_%s",cutName))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(),"NGlobalTracksTowardLeading","Number of global tracks (toward relative to leading) per event",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracksTransverseLeading","Number of global tracks (Transverse relative to leading) per event",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracksAwayLeading","Number of global tracks (Away relative to leading) per event",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksTowardLeading","Number of global tracks (all regions) vs Number of global tracks (Toward relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksTransverseLeading","Number of global tracks (all regions) vs Number of global tracks (Transverse relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_NGlobalTracksAwayLeading","Number of global tracks (all regions) vs Number of global tracks (Away relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        if (isMC) {
          man->AddHistogram(classStr.Data(),"MultiplicityTowardLeading","Multiplicity (toward relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TowardLeading_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (toward relative to leading) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Toward + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TowardLeading","Multiplicity in |#eta|<0.9 vs global tracks (toward relative to leading) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksToward + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Toward + 1);
          man->AddHistogram(classStr.Data(),"MultiplicityTransverseLeading","Multiplicity (Transverse relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TransverseLeading_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (Transverse relative to leading) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Transverse + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TransverseLeading","Multiplicity in |#eta|<0.9 vs global tracks (Transverse relative to leading) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksTransverse + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Transverse + 1);
          man->AddHistogram(classStr.Data(),"MultiplicityAwayLeading","Multiplicity (Away relative to leading) in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_AwayLeading_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (Away relative to leading) (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Away + 1);
          man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_AwayLeading","Multiplicity in |#eta|<0.9 vs global tracks (Away relative to leading) (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracksAway + AliReducedVarManager::kNMaxCutsGlobalTracks + icut,nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09Away + 1);

        }
        continue;
      }
    }

    // Track histograms
    if(classStr.Contains("TrackITSclusterMap_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
			6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
			180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
			100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE,
      //                  1000, -0.5, 0.5, AliReducedVarManager::kDcaXY, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE,
      //                  1800, -1.0, 1.0, AliReducedVarManager::kDcaZ, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackTPCclusterMap_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
			8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
			180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
			100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
			100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      continue;
    }  // end of TPCclusterMap histogram definitions
    
    TString trkStatusNames = "";
    for(Int_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
      trkStatusNames += AliReducedVarManager::fgkTrackingStatusNames[iflag];
      trkStatusNames += ";";
    }
    if(classStr.Contains("TrackStatusFlags_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
			AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, 
			0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trkStatusNames.Data());
      /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TrackingFlags","Corrected TPC N_{#sigma} electron vs. inner param P vs tracking flags;;",kFALSE,
	50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data());
      man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_TrackingFlags_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
	43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data(),  AliReducedVarManager::kTPCncls);*/
      continue;
    }
    
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kP, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 1000, -0.5, 0.5, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 1000, -0.5, 0.5, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "TPCNSigEl", "NSigEl TPC", kFALSE, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSig);
      man->AddHistogram(classStr.Data(), "TPCNSigEl_Pin", "NSigEl TPC vs Pin", kFALSE, 90, 1., 10., AliReducedVarManager::kPin, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSig);
      man->AddHistogram(classStr.Data(), "TPCNSigEl_Eta", "NSigEl TPC vs Eta", kFALSE, 18, -0.9, 0.9, AliReducedVarManager::kEta, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSig);
      man->AddHistogram(classStr.Data(), "TPCNSigElCorr", "NSigEl TPC", kFALSE, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSigCorrected);
      man->AddHistogram(classStr.Data(), "TPCNSigElCorr_Pin", "NSigEl TPC vs Pin", kFALSE, 90, 1., 10., AliReducedVarManager::kPin, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSigCorrected);
      man->AddHistogram(classStr.Data(), "TPCNSigElCorr_Eta", "NSigEl TPC vs Eta", kFALSE, 18, -0.9, 0.9, AliReducedVarManager::kEta, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSigCorrected);
      man->AddHistogram(classStr.Data(), "TOFNSigEl", "NSigEl TOF", kFALSE, 2000, -50.0, 50.0, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TOFNSigEl_Pin", "NSigEl TOF vs Pin", kFALSE, 90, 1., 10., AliReducedVarManager::kPin, 400, -20.0, 20.0, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TOFNSigEl_Eta", "NSigEl TOF vs Eta", kFALSE, 18, -0.9, 0.9, AliReducedVarManager::kEta, 400, -20.0, 20.0, AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TRDElLQ2D", "TRD electron likelihood", kFALSE, 1000, 0., 1., AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TRDElLQ2D_Pin", "TRD electron likelihood vs Pin", kFALSE, 90, 1., 10., AliReducedVarManager::kPin, 200, 0., 1., AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TRDElLQ2D_Eta", "TRD electron likelihood vs Eta", kFALSE, 18, -0.9, 0.9, AliReducedVarManager::kEta, 200, 0, 1., AliReducedVarManager::kTRDpidProbabilitiesLQ2D+AliReducedVarManager::kElectron);

        
      if(classStr.Contains("MCTruth")) {
          man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "PtRec_PtMC", "p_{T} MC vs p_{T} reconstructed", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC, 150, 0., 15.0, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "PhiMC", "#varphi MC", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC);
          man->AddHistogram(classStr.Data(), "PhiRec_PhiMC", "#varphi MC vs #varphi reconstructed", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC, 180, 0., 6.3, AliReducedVarManager::kPhi);
          man->AddHistogram(classStr.Data(), "EtaMC", "#eta MC", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC);
          man->AddHistogram(classStr.Data(), "EtaRec_EtaMC", "#eta MC vs #eta reconstructed", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC, 100, -1.0, 1.0, AliReducedVarManager::kEta);          
          man->AddHistogram(classStr.Data(), "PDGcode0", "PDG code of the track", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC);
          man->AddHistogram(classStr.Data(), "PDGcode1", "PDG code of the track's mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+1);
          man->AddHistogram(classStr.Data(), "PDGcode2", "PDG code of the track's grand-mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+2);
          man->AddHistogram(classStr.Data(), "PDGcode3", "PDG code of the track's grand-grand mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+3);
      }
      continue;
    }  // end if "TrackQA"
    
    if(classStr.Contains("PairQualityFlags")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString pairQualityFlagNames = " ;K^{0}_{S}#rightarrow#pi^{+}#pi^{-};#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};#gamma#rightarrow e^{+}e^{-};";
       man->AddHistogram(classStr.Data(), "PairQualityFlags", "Pair quality flags;;", kFALSE,
                         32, -0.5, 31.5, AliReducedVarManager::kPairQualityFlag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, pairQualityFlagNames.Data());
       continue;
    }
    
    Double_t massBinWidth = 0.01;     // *GeV/c^2
    Double_t massRange[2] = {0.0,5.0};
    Int_t nMassBins = TMath::Nint((massRange[1]-massRange[0])/massBinWidth);
    
    TString candidateNames = "#gamma#rightarrow e^{+}e^{-};K^{0}_{S}#rightarrow#pi^{+}#pi^{-};";
    candidateNames += "#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};";
    
    // Histograms for pairs
    if(classStr.Contains("Pair_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "CandidateId", "Candidate id", kFALSE,
                         5, -0.5, 4.5, AliReducedVarManager::kCandidateId, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, candidateNames.Data());
       man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
       man->AddHistogram(classStr.Data(), "PairChi2", "Pair #chi^{2}", kFALSE, 200, 0.0, 50, AliReducedVarManager::kPairChisquare);
       man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMass);
       man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE, 2000, -0.2, 0.2, AliReducedVarManager::kPseudoProperDecayTime);
      /* man->AddHistogram(classStr.Data(), "Mass_V0K0s", "Invariant mass, K^{0}_{s} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0);
       man->AddHistogram(classStr.Data(), "Mass_V0Lambda", "Invariant mass, #Lambda^{0} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+1);
       man->AddHistogram(classStr.Data(), "Mass_V0ALambda", "Invariant mass, #bar{#Lambda^{0}} assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+2);
       man->AddHistogram(classStr.Data(), "Mass_V0Gamma", "Invariant mass, #gamma conversion assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+3);*/ 
       man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
       man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPx);
       man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPy);
       man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPz);
       man->AddHistogram(classStr.Data(), "P", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kP);
       man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
       man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEta);
       man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
       man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kTheta);
       man->AddHistogram(classStr.Data(), "LxyOrR", "L_{xy}/Decay Radius", kFALSE, 1000, -10.0, 10.0, AliReducedVarManager::kPairLxy);
       man->AddHistogram(classStr.Data(), "OpeningAngle", "Opening angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairOpeningAngle);
       man->AddHistogram(classStr.Data(), "PointingAngle", "Pointing angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairPointingAngle);
       man->AddHistogram(classStr.Data(), "PsiPair", "Psi pair", kFALSE, 200, -1.5, 1.5, AliReducedVarManager::kPairPsi);
       man->AddHistogram(classStr.Data(), "PsiPair_DeltaPhi", "Psi pair vs delta phi", kFALSE, 100, -0.5, 0.5, AliReducedVarManager::kPairDeltaPhi, 200, -1.5, 1.5, AliReducedVarManager::kPairPsi);      
       continue;
    }   // end if for Pair classes of histograms 
    
    if(classStr.Contains("PairME")){
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, 5000, 0.0, 5.0, AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Mass_Pt", "Pt vs Invariant mass", kFALSE, 5000, 0.0, 5.0, AliReducedVarManager::kMass, 150, 0.0, 15.0, AliReducedVarManager::kPt);
      // man->AddHistogram(classStr.Data(), "Multip", "Invariant mass vs multiplicity", kFALSE, 150, 0., 150., AliReducedVarManager::kNtracksSelected, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 150, 0.0, 15.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 300, -1.5, 1.5, AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Eta", "Eta", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 630, 0.0, 6.3, AliReducedVarManager::kPhi);
      continue;
    } 
    if(classStr.Contains("PionPrimary")) {
        man->AddHistClass(classStr.Data());
        Double_t ptbins[58] = {0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 
                              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5,
                                                                                                           7, 8, 10, 13, 20, 30, 50, 80, 100, 200};
        Double_t* binsPt = ptbins;
        Double_t multbins[11] = {0., 0.92, 4.6, 9.2, 13.8, 18.4, 27.6, 36.8, 46, 64.5, 100};
        Double_t* binsMult = multbins;

        man->AddHistogram(classStr.Data(), "Mult_Pt", "", kFALSE, 57, binsPt, AliReducedVarManager::kPt, 10, binsMult, AliReducedVarManager::AliReducedVarManager::kMultEstimatorPercentileV0M);
        continue;
    }

    if(classStr.Contains("MCTRUTH")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassMC);
      man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTimeMC);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1500, 0.0, 30.0, AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "Mass_Pt", "Mass vs pt", kFALSE, 1500, 0.0, 30.0, AliReducedVarManager::kPtMC, 200, 2., 4., AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "PtMeas", "", kFALSE, 1500, 0.0, 30.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 1500, 0.0, 10.0, AliReducedVarManager::kPxMC);
      man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 1500, 0.0, 10.0, AliReducedVarManager::kPyMC);
      man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 1500, 0.0, 10.0, AliReducedVarManager::kPzMC);
      man->AddHistogram(classStr.Data(), "P", "", kFALSE, 1500, 0.0, 10.0, AliReducedVarManager::kPMC);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRapMC);
      man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEtaMC);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kThetaMC);
      man->AddHistogram(classStr.Data(), "PDG", "PDG value", kFALSE, 300, -600, 600, AliReducedVarManager::kPdgMC);
      man->AddHistogram(classStr.Data(), "PDGmother", "PDG value jpsi mother", kFALSE, 300, -6000, 6000, AliReducedVarManager::kPdgMC+1);
      man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Number of jpsi vs measured multiplicity", kFALSE,  nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks);
      man->AddHistogram(classStr.Data(), "multDist_evt_gen", "Number of jpsi vs true multiplicity", kFALSE,  nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
      man->AddHistogram(classStr.Data(), "multCorrel_evt", "True multiplicity vs measured multiplicity for Jpsi", kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
      man->AddHistogram(classStr.Data(), "MeanPt_vtxNcontrib", "mean pt vs vertex n contributors", kTRUE,  nMultBins, minMult, maxMult, AliReducedVarManager::kNVtxContributors, 3000, 0, 30., AliReducedVarManager::kPtMC);
      continue;
    }   // end if for Pair classes of histograms 

    for (int icut = 0; icut < task->GetNMeasMultCuts(); icut++) {
      const char* cutName = task->GetMeasMultcutName(icut);
      Double_t ptbins[10] = {0., 1., 2., 3., 4., 5., 6., 8., 12., 100.};
      Double_t* binsPt = ptbins;
      Double_t multbins[nMultBinspp]; for (int b = 0; b <= nMultBinspp; b++) {multbins[b] = minMultpp+b*(maxMultpp-minMultpp)/nMultBinspp;}
      Double_t* binsMult = multbins;
      if(classStr.Contains(Form("JpsiPtMultCorrel_%s", cutName))) {// Counts for all truth particles
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        if (isMC) {
          man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen", "True Pt spectrum (all generated) vs true multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_meas", "True Pt spectrum (all generated) vs meas multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_gen", "True multiplicity distribution", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09 + 2, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_gen_trig", "True multiplicity distribution after trigger", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09 + 2, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "Pt_multCorrel_evt", "True multiplicity vs measured multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, nMultBinspp, binsMult, AliReducedVarManager::kMCNch09+2, 9, binsPt, AliReducedVarManager::kPtMC); 
          man->AddHistogram(classStr.Data(), "MeanPt_multDist_evt_meas", "mean pt vs measured multiplicity distribution", kTRUE,  nMultBinspp, minMult, maxMult, AliReducedVarManager::kNGlobalTracks+icut, 3000, 0, 30., AliReducedVarManager::kPtMC);
        }
        continue;
      }
      if(classStr.Contains(Form("JpsiPtMultCorrelMeas_%s", cutName))) {// Counts, filled only when measured
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, 9, binsPt, AliReducedVarManager::kPt); //2D unfolding
        if(isMC) {
          man->AddHistogram(classStr.Data(), "multCorrel_prim", "True Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut,  9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "ptCorrel_prim", "True Pt spectrum vs measured pt spectrum", kFALSE, 9, binsPt, AliReducedVarManager::kPt, 9, binsPt, AliReducedVarManager::kPtMC); 
          man->AddHistogram(classStr.Data(), "multPtSpec_prim_meas", "True Pt spectrum vs true multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09, 9, binsPt, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "multPtSpec_trk_inter", "Pt measured spectrum vs true multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09, 9, binsPt, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, 9, binsPt, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_gen", "True multiplicity distribution", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09 + 2, 9, binsPt, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "Pt_multDist_evt_gen_trig", "True multiplicity distribution after trigger", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kMCNch09 + 2, 9, binsPt, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "Pt_multCorrel_evt", "True multiplicity vs measured multiplicity", kFALSE,  nMultBinspp, binsMult, AliReducedVarManager::kNGlobalTracks+icut, nMultBinspp, binsMult, AliReducedVarManager::kMCNch09+2, 9, binsPt, AliReducedVarManager::kPt); 
        }
        continue;
      }
    }

    if(classStr.Contains("DeltaPhi_JpsiTruth_JpsiCandidate")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "DeltaPhi", "Delta phi between true Jpsi and best Jpsi candidate", kFALSE,  400, 0., 2*M_PI, AliReducedVarManager::kPhiJpsiMCTruth);
    }

    for(int cutMode = 0; cutMode < 4*nEstimators; cutMode++) {
      if(classStr.Contains(Form("pp_13TeV_Data_cutMode_%d", 100+cutMode))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4);
        //man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
      }
      if(classStr.Contains(Form("pp_13TeV_MC_cutMode_%d",100+cutMode))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4);
        //man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
        man->AddHistogram(classStr.Data(), "multDist_evt_gen", "True multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09 + 2);
        man->AddHistogram(classStr.Data(), "multDist_evt_gen_trig", "True multiplicity distribution after trigger", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09 + 2);
        man->AddHistogram(classStr.Data(), "multCorrel_evt", "True multiplicity vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09+2); 
        /*man->AddHistogram(classStr.Data(), "multCorrel_prim", "True Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09); //2D unfolding - not used - filling is false
        man->AddHistogram(classStr.Data(), "ptCorrel_prim", "True Pt spectrum vs measured pt spectrum", kFALSE, 100, 0., 10., AliReducedVarManager::kPt, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen", "True Pt spectrum (all generated primaries) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_evtloss", "True Pt spectrum (generated primaries from rejected events) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_notrig", "True Pt spectrum (generated primaries from untriggered events) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_meas", "True Pt spectrum vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_prim_meas", "True Pt spectrum primaries vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_sec_meas", "True Pt spectrum secondaries vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas_evtcont", "True Pt spectrum (from contaminated events) vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_inter", "Pt measured spectrum vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
      */
      }

      if(classStr.Contains(Form("pPb_5TeV_Data_cutMode_%d",100+cutMode))) {
         man->AddHistClass(classStr.Data());
         cout << classStr.Data() << endl;
         man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4);
         man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
      }
      if(classStr.Contains(Form("pPb_5TeV_MC_cutMode_%d",100+cutMode))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4);
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
        man->AddHistogram(classStr.Data(), "multDist_evt_gen", "True multiplicity distribution", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09);
        man->AddHistogram(classStr.Data(), "multDist_evt_gen_trig", "True multiplicity distribution after trigger", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09+1);
        man->AddHistogram(classStr.Data(), "multCorrel_evt", "True multiplicity vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09+2); 
        man->AddHistogram(classStr.Data(), "multCorrel_prim", "True Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09); //2D unfolding - not used - filling is false
        man->AddHistogram(classStr.Data(), "ptCorrel_prim", "True Pt spectrum vs measured pt spectrum", kFALSE, 100, 0., 10., AliReducedVarManager::kPt, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen", "True Pt spectrum (all generated primaries) vs true multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_evtloss", "True Pt spectrum (generated primaries from rejected events) vs true multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09+cutMode/4, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_notrig", "True Pt spectrum (generated primaries from untriggered events) vs true multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09+cutMode/4, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_meas", "True Pt spectrum vs true multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_prim_meas", "True Pt spectrum primaries vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
         man->AddHistogram(classStr.Data(), "multPtSpec_trk_sec_meas", "True Pt spectrum secondaries vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
         man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas_evtcont", "True Pt spectrum (from contaminated events) vs measured multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
         man->AddHistogram(classStr.Data(), "multPtSpec_trk_inter", "Pt measured spectrum vs true multiplicity", kFALSE,  nMultBinspPb, minMultpPb, maxMultpPb, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
      }

    }

  }  // end loop over histogram classes
}
