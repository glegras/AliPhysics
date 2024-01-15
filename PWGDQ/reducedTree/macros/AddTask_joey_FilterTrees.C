#include <iostream>
#include <fstream>
#include <string>

#include "TF1.h"


void Setup(AliReducedAnalysisFilterTrees* processor, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/);
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/,  TString prod /*="LHC10h"*/);
void SetupMixingHandler(AliReducedAnalysisFilterTrees* task);
void DefineHistograms(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/);


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


//__________________________________________________________________________________________
AliAnalysisTask* AddTask_glegras_FilterTrees(Bool_t isAliRoot=kTRUE, Int_t runMode=1, Bool_t isMC = kFALSE, Bool_t isInjected = kFALSE, Bool_t ispp = kTRUE, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_glegras_FilterTrees(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisFilterTrees* filterTask = new AliReducedAnalysisFilterTrees("FilterTrees","filter DST trees");
  filterTask->Init();
  filterTask->SetFilteredTreeWritingOption(AliReducedAnalysisTaskSE::kBaseEventsWithBaseTracks);
  filterTask->SetWriteFilteredTracks(kFALSE);
  filterTask->SetWriteFilteredPairs(kFALSE);
  filterTask->SetBuildCandidatePairs(AliReducedPairInfo::kJpsiToEE);
  filterTask->SetBuildCandidateLikePairs(kTRUE);
  filterTask->SetRunCandidatePrefilter(kTRUE);
  filterTask->SetRunCandidatePrefilterOnSameCharge(kFALSE);
  filterTask->SetRejectEmptyEvents(kTRUE);
  //filterTask->SetRunEventMixing(kTRUE);
  //filterTask->SetComputeMult(kTRUE);
  filterTask->SetRunOverMC(isMC);


  if (isMC && isInjected) {
    TF1* fPtJpsi = new TF1("fPtJpsi","[0]*x/pow(1+(x/[1])*(x/[1]),[2])",0,30);
    fPtJpsi->SetParameters(1,4.09,3.04); // arxiv:2108.01906
    TFile* filePtWeights = new TFile("/gluster1/glegras/InjectedJpsiPtWeights.root","read");
    TH1F* hpt; filePtWeights->GetObject("Pt",hpt);
    if(!hpt) {
       filePtWeights = new TFile("~/alice/AliPhysics/PWGDQ/reducedTree/macros/InjectedJpsiPtWeights.root","read");
       filePtWeights->GetObject("Pt",hpt);
    }
    TH1F* hPtWeights = new TH1F("hPtWeights","Weights for rescaling injected Jpsi",hpt->GetNbinsX(),hpt->GetBinLowEdge(1),hpt->GetBinCenter(hpt->GetNbinsX())+hpt->GetBinCenter(hpt->GetNbinsX())/2-hpt->GetBinLowEdge(hpt->GetNbinsX())/2);
    for (int n=1;n<hpt->GetNbinsX()+1;n++) hPtWeights->SetBinContent(n,fPtJpsi->Eval(hPtWeights->GetBinCenter(n))/hpt->GetBinContent(n));
    hPtWeights->Scale(1./hPtWeights->GetMaximum());

    filterTask->SetMCJpsiPtWeights(hPtWeights); //only if MC and injected Jpsi
  }
/*
  TFile* fileJpsiMass = new TFile("/gluster1/glegras/JpsiMassDist.root","read");
  TF1* fJpsiMass; fileJpsiMass->GetObject("signalFitCB",fJpsiMass);
  if(!fJpsiMass) {
    fileJpsiMass = new TFile("~/alice/AliPhysics/PWGDQ/reducedTree/macros/JpsiMassDist.root","read");
    fileJpsiMass->GetObject("signalFitCB",fJpsiMass);
  }
  filterTask->SetJpsiMassDist(fJpsiMass);
  if(!filterTask->GetJpsiMassDist()) cout<<"No jpsi mass distribution setup"<<endl;
*/

  Setup(filterTask, ispp, prod);
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


//_________________________________________________________________
void Setup(AliReducedAnalysisFilterTrees* processor, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  // 
  
  bool isMC = processor->GetRunOverMC();

  TF1* fV0cut = new TF1("fV0cut", V0cut, 2e5, 3e5);

  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, 0., fV0cut, kTRUE, AliReducedVarManager::kRunNo, 0., 0., kTRUE, AliReducedVarManager::kHighMultV0Triggered, 0.5, 1.5);
  //evCut1->AddCut(AliReducedVarManager::kINT7Triggered, 0.1, 2.);
  //evCut1->AddCut(AliReducedVarManager::kHighMultV0Triggered, -0.1, 0.1);
  //evCut1->AddCut(AliReducedVarManager::kHighMultSPDTriggered, -0.1, 0.1);
  processor->AddEventCut(evCut1);

  //Set track cuts for multiplicity info
/*
  // Measured track cuts (max. 8 cuts) - names should not be included one in the other
  AliReducedTrackCut* measMultCut = new AliReducedTrackCut("standardCut","Cut for measured mult");
  measMultCut->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCut->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  if(ispp || isMC) measMultCut->SetTrackFilterBit(8); //SPDAny recquirement
  else measMultCut->SetTrackFilterBit(3); //charged tracks with TPC cluster recquirement
  
  measMultCut->SetTrackFilterBit(6);
  processor->AddMeasuredMultTrackCut(measMultCut);

  AliReducedTrackCut* measMultCutPt = new AliReducedTrackCut("Pt02","Cut for measured mult"); 
  measMultCutPt->AddCut(AliReducedVarManager::kPt, 0.2, 1e5);
  measMultCutPt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutPt->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  if(ispp || isMC) measMultCutPt->SetTrackFilterBit(8); //SPDAny recquirement
  else measMultCutPt->SetTrackFilterBit(3); //charged tracks with TPC cluster recquirement
  processor->AddMeasuredMultTrackCut(measMultCutPt);

  AliReducedTrackCut* measMultCutnoSPDAny = new AliReducedTrackCut("noSPDAny","Cut for measured mult"); 
  measMultCutnoSPDAny->AddCut(AliReducedVarManager::kPt, 0.15, 1e5);
  measMultCutnoSPDAny->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutnoSPDAny->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  if(ispp || isMC) measMultCutnoSPDAny->SetTrackFilterBit(7); 
  else measMultCutnoSPDAny->SetTrackFilterBit(3); //charged tracks with TPC cluster recquirement
  processor->AddMeasuredMultTrackCut(measMultCutnoSPDAny);

  AliReducedTrackCut* measMultCutnoSPDAnyPt = new AliReducedTrackCut("cut4","Cut for measured mult"); 
  measMultCutnoSPDAnyPt->AddCut(AliReducedVarManager::kPt, 0.2, 1e5);
  measMultCutnoSPDAnyPt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  measMultCutnoSPDAnyPt->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  if(ispp || isMC) measMultCutnoSPDAnyPt->SetTrackFilterBit(7);
  else measMultCutnoSPDAnyPt->SetTrackFilterBit(3); //charged tracks with TPC cluster recquirement
  processor->AddMeasuredMultTrackCut(measMultCutnoSPDAnyPt);


  // MC cut on true multiplicity  
  AliReducedTrackCut* trueMultCut = new AliReducedTrackCut("trueMult","Cut for true mult");
  trueMultCut->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueMultCut->AddCut(AliReducedVarManager::kCharge, -0.1, 0.1, kTRUE);
  trueMultCut->SetRejectPureMC(kFALSE);
  trueMultCut->SetMCFilterBit(13); // physical primaries
*/
  processor->AddTrueMultTrackCut(trueMultCut);

  // Set track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kP, 1.,100.);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  //standardCut->SetRequestTOFout();
  standardCut->SetRequestSPDany();
  standardCut->SetRejectPureMC(kTRUE);
  standardCut->SetTrackFilterBit(0);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  processor->AddTrackCut(standardCut); 
  processor->AddCandidateLeg1Cut(standardCut);
  
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt07","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.2,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  //if(ispp) prefTrackCut1->SetTrackFilterBit(2); //electron prefilter
  if(ispp) prefTrackCut1->SetTrackFilterBit(1); //electron prefilter
  prefTrackCut1->SetRejectPureMC(kTRUE);
  processor->AddCandidateLeg1PrefilterCut(prefTrackCut1);  
  
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddCandidateLeg1PairPrefilterCut(prefPairCut);

  
  AliReducedTrackCut* pairCut = new AliReducedTrackCut("JpsiCut","Pt pair selection");
  pairCut->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddCandidatePairCut(pairCut);

  AliReducedTrackCut* jpsiMCCut = new AliReducedTrackCut("jpsiCutMC","Rapidity JpsiSelection");
  jpsiMCCut->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  jpsiMCCut->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC, 442.5, 443.5); 
  //jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC+1, 442.5, 443.5, kTRUE);
  jpsiMCCut->SetRejectPureMC(kFALSE);
  AliReducedTrackCut* electronMCCut = new AliReducedTrackCut("standardMC","Pt electron Selection");
  processor->AddJpsiMotherMCCut(jpsiMCCut,electronMCCut);
  
  AliReducedTrackCut* jpsiMCCutEta = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEta");
  AliReducedTrackCut* electronMCCutEta = new AliReducedTrackCut("standardMCEta","Eta electron Selection");
  electronMCCutEta->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEta,electronMCCutEta);


  AliReducedTrackCut* jpsiMCCutEtaPt = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEtaPt");
  AliReducedTrackCut* electronMCCutEtaPt = new AliReducedTrackCut("standardMCEtaPt","Eta and Pt electron Selection");
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kPt, 1., 1e5);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEtaPt,electronMCCutEtaPt);

  
  //SetupMixingHandler(processor);
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
/*
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
   Float_t zLims[nZbinLimits] = {-10.,  
   // -8.,  -6.,  -4.,  -2., 0.,  2.,  4.,  6.,  8., 
   10. };
   //r(Int_t i=0;i<=10;++i) zLims[nZbinLimits] = i
   
  handler->AddMixingVariable(AliReducedVarManager::kVtxZ, nZbinLimits, zLims); 
  TString histClassNames = handler->GetHistClassNames();
  TObjArray* histClassArr = histClassNames.Tokenize(";");

  task->SetRunEventMixingMult(kFALSE);

  const int nMultBins = 2; Int_t multBins[nMultBins+1] = {0,40,150};
  task->SetMultBinsMixing(nMultBins,multBins);
  for (int i = 0; i<nMultBins;i++){
    AliMixingHandler* handlerMult = (AliMixingHandler*) handler->Clone();
    TString histNamesMult = "";
    for (int j = 0; j<histClassArr->GetEntries(); j++) {
        histNamesMult += histClassArr->At(j)->GetName();
        histNamesMult += Form("_%d_%d;",multBins[i],multBins[i+1]-1);
    }
    handlerMult->SetHistClassNames(histNamesMult);   
    handlerMult->SetHistogramManager(task->GetHistogramManager());
    task->AddMixingHandler(handlerMult);
  }
}
*/

//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisFilterTrees* task, Bool_t ispp/*=kTRUE*/, TString prod /*="LHC10h"*/) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts
  //           are used in the histogram lists
  // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
  //  TODO: make needed changes such that this becomes less prone to mistakes
   /*
  TString runNumbers = ""; const char* fileRuns; const char* fileRuns2;
  int nRuns = 0;
  if(ispp) {
    fileRuns = "/home/glegras/alice/data/runLists/runsMC.txt";
    fileRuns2 = "/gluster1/glegras/runsMC.txt";
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

  AliReducedVarManager::SetRunNumbers(runNumbers);*/

  AliHistogramManager* man = task->GetHistogramManager(); 
   
  TString histClasses = "";
  //histClasses += "Event_BeforeCuts;";
  histClasses += "Event_AfterCuts;";

  //int nEstimators = (task->GetComputeMult() ? task->GetNMeasMultCuts() : 0);
/*
  for (int icut = 0; icut<nEstimators; icut++) {
    histClasses += Form("EventMult_%s_Inclusive;",task->GetMeasMultcutName(icut));
    histClasses += Form("EventMult_%s_MB;",task->GetMeasMultcutName(icut));
    histClasses += Form("EventMult_%s_HM;",task->GetMeasMultcutName(icut));
  }

  histClasses += "Multiplicity_Regions;";
*/
  bool isMC=task->GetRunOverMC();
  
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
      histClasses += Form("Track_LEG1_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,1));
      histClasses += Form("Track_LEG2_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,2));
      if(task->GetRunCandidatePrefilter()) {
        histClasses += Form("Track_LEG1_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,1));
        histClasses += Form("Track_LEG2_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,2));
      }
      histClasses += Form("Pair_Candidate12_%s%s;", 
			  task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      if(task->GetBuildCandidateLikePairs()) {
        	histClasses += Form("Pair_Candidate11_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
	        histClasses += Form("Pair_Candidate22_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      }/*
      if(task->GetRunEventMixing()) {
          histClasses += Form("PairMEPM_%s%s;PairMEPP_%s%s;PairMEMM_%s%s;",
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));

      }
      if(task->GetRunEventMixingMult()) {
        int nBins = task->GetNMultBinsMixing();
        for (int j = 0; j<nBins; j++) {
          int lowMult = task->GetMultBinsMixing(j);int highMult = task->GetMultBinsMixing(j+1)-1;
          histClasses += Form("PairMEPM_%s%s_%d_%d;PairMEPP_%s%s_%d_%d;PairMEMM_%s%s_%d_%d;",
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult,highMult,
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult,highMult,
          task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""),
          lowMult, highMult);
        }
      }*/
    }
    /*if(task->GetRunCandidatePrefilter()) {
        histClasses += "Track_LEG1_PrefilterTrack;";
        histClasses += "Track_LEG2_PrefilterTrack;";
	}*/
  }

  if(isMC) {
    for (int i=0; i<task->GetNJpsiMotherMCCuts();i++) 
      histClasses += Form("PureMCTRUTH_AfterSelection_%s;", task->GetJpsiMotherMCcutName(i));
    
    //histClasses += Form("PureMCTRUTH_DetectedDaughters_%s;", task->GetJpsiMotherMCcutName(task->GetNJpsiMotherMCCuts()-1));
  } 
   
 /*
  for (int cutMode = 0; cutMode < 4*nEstimators; cutMode++) { 
    // 0 for MB multiplicity  unfolding
    // 1 for HM multiplicity unfolding
    // 2 for MB+HM multiplicity unfolding
    // 3 is for Jpsi vs multiplicity unfolding (HM+MB)
    if(ispp){
      if(isMC) histClasses += Form("pp_13TeV_MC_cutMode_%d;",100+cutMode);
      else histClasses += Form("pp_13TeV_Data_cutMode_%d;",100+cutMode);
    }
    else {
      if(isMC) histClasses += Form("pPb_5TeV_MC_cutMode_%d;",100+cutMode);
      else histClasses += Form("pPb_5TeV_Data_cutMode_%d;",100+cutMode);
    }
  }
*/
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};


  int nMultBins; float minMult; float maxMult;
  int nMultBinspp = 201; float minMultpp = -0.5; float maxMultpp = 200.5;
  int nMultBinspPb = 301; float minMultpPb = -0.5; float maxMultpPb = 300.5;

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
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO,
			50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO_TimeFromSOR","Centrality(VZERO) vs time from SOR",kFALSE,
                        90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,30000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksTotal_BeamIntensity0_prof","Number of total tracks per event",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,20000.,AliReducedVarManager::kNtracksTotal);

      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,90, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors_BeamIntensity0","",kTRUE, 100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      Bool_t isCalibrated = kTRUE;
      man->AddHistogram(classStr.Data(), "VZEROA_NEmptyChannels_VtxCent_prof", "No. VZERO-A empty channels per event vs. centrality SPD and vertex Z", kTRUE,
			24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROAemptyChannels);
      man->AddHistogram(classStr.Data(), "VZEROC_NEmptyChannels_VtxCent_prof", "No. VZERO-C empty channels per event vs. centrality SPD and vertex Z", kTRUE,
			24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROCemptyChannels);
      continue;
    }  // end if className contains "Event" 
    /*        
    for (int icut = 0; icut<task->GetNMeasMultCuts(); icut++) {
      const char* cutName = task->GetMeasMultcutName(icut);
      if(classStr.Contains(Form("EventMult_%s",cutName))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(),"NGlobalTracks","Number of global tracks per event",kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"SPDntracklets_V0Mult","SPD ntracklets vs V0 Mult",kTRUE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kSPDntracklets);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_V0Mult","Number of global tracks per event vs V0 Mult",kTRUE, 800, 0, 800, AliReducedVarManager::kVZEROTotalMult, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"Multiplicity","Multiplicity in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09);
        man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks_TProfile","Mean multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kTRUE,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09);
        man->AddHistogram(classStr.Data(),"Multiplicity_NGlobalTracks","Multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kFALSE,nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks + icut,nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_vtxz_TProfile","Mean number of global tracks (in |#eta|<0.9) vs vtxz",kTRUE,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_vtxz","Global tracks (in |#eta|<0.9) vs vtxz",kFALSE,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_RunNumber_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run number",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        //man->AddHistogram(classStr.Data(),"NGlobalTracks_RunNumber_vtxz_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run number vs vtxz",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        //man->AddHistogram(classStr.Data(),"NGlobalTracks_RunNumber_vtxz","Global tracks (in |#eta|<0.9) vs Run number vs vtxz",kFALSE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_RunID_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run ID",kTRUE, nRuns, 0, nRuns, AliReducedVarManager::kRunID,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        man->AddHistogram(classStr.Data(),"NGlobalTracks_RunID_vtxz_TProfile","Mean number of global tracks (in |#eta|<0.9) vs Run ID vs vtxz",kTRUE, nRuns, 0, nRuns, AliReducedVarManager::kRunID,100,-10.,10.,AliReducedVarManager::kVtxZ,nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracks + icut);
        continue;
      }
    }

    if(classStr.Contains("Multiplicity_Regions")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      //Toward region 
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward","Number of Toward global tracks per event",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"MultiplicityToward","Multiplicity Toward in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
      man->AddHistogram(classStr.Data(),"MultiplicityToward_NGlobalTracksToward_TProfile","Mean Toward multiplicity in |#eta|<0.9 vs Toward global tracks (in |#eta|<0.9)",kTRUE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
      man->AddHistogram(classStr.Data(),"MultiplicityToward_NGlobalTrackTowards","Multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Toward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_vtxz_TProfile","Mean number of Toward global tracks (in |#eta|<0.9) vs vtxz",kTRUE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_vtxz","Toward Global tracks (in |#eta|<0.9) vs vtxz",kFALSE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_RunNumber_TProfile","Mean number of Toward global tracks (in |#eta|<0.9) vs Run number",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_RunNumber_vtxz_TProfile","Mean number of Toward global tracks (in |#eta|<0.9) vs Run number vs vtxz",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      //man->AddHistogram(classStr.Data(),"NGlobalTracksToward_RunNumber_vtxz","Global tracks (in |#eta|<0.9) vs Run number vs vtxz",kFALSE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_RunID_TProfile","Mean number of Toward global tracks (in |#eta|<0.9) vs Run ID",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);
      man->AddHistogram(classStr.Data(),"NGlobalTracksToward_RunID_vtxz_TProfile","Mean number of Toward global tracks (in |#eta|<0.9) vs Run ID vs vtxz",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksToward);


      //Transverse region 
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse","Number of Transverse global tracks per event",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      man->AddHistogram(classStr.Data(),"MultiplicityTransverse","Multiplicity Transverse in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
      man->AddHistogram(classStr.Data(),"MultiplicityTransverse_NGlobalTracksTransverse_TProfile","Mean Transverse multiplicity in |#eta|<0.9 vs Transverse global tracks (in |#eta|<0.9)",kTRUE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
      man->AddHistogram(classStr.Data(),"MultiplicityTransverse_NGlobalTrackTransverses","Multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Transverse);
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_vtxz_TProfile","Mean number of Transverse global tracks (in |#eta|<0.9) vs vtxz",kTRUE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_vtxz","Transverse Global tracks (in |#eta|<0.9) vs vtxz",kFALSE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_RunNumber_TProfile","Mean number of Transverse global tracks (in |#eta|<0.9) vs Run number",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_RunNumber_vtxz_TProfile","Mean number of Transverse global tracks (in |#eta|<0.9) vs Run number vs vtxz",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      //man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_RunNumber_vtxz","Global tracks (in |#eta|<0.9) vs Run number vs vtxz",kFALSE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);      
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_RunID_TProfile","Mean number of Transverse global tracks (in |#eta|<0.9) vs Run ID",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      man->AddHistogram(classStr.Data(),"NGlobalTracksTransverse_RunID_vtxz_TProfile","Mean number of Transverse global tracks (in |#eta|<0.9) vs Run ID vs vtxz",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksTransverse);
      
      //Away region 
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway","Number of Away global tracks per event",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"MultiplicityAway","Multiplicity Away in |#eta|<0.9",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
      man->AddHistogram(classStr.Data(),"MultiplicityAway_NGlobalTracksAway_TProfile","Mean Away multiplicity in |#eta|<0.9 vs Away global tracks (in |#eta|<0.9)",kTRUE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
      man->AddHistogram(classStr.Data(),"MultiplicityAway_NGlobalTrackAways","Multiplicity in |#eta|<0.9 vs global tracks (in |#eta|<0.9)",kFALSE, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway, nMultBins, minMult, maxMult,AliReducedVarManager::kMCNch09Away);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_vtxz_TProfile","Mean number of Away global tracks (in |#eta|<0.9) vs vtxz",kTRUE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_vtxz","Away Global tracks (in |#eta|<0.9) vs vtxz",kFALSE,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_RunNumber_TProfile","Mean number of Away global tracks (in |#eta|<0.9) vs Run number",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_RunNumber_vtxz_TProfile","Mean number of Away global tracks (in |#eta|<0.9) vs Run number vs vtxz",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      //man->AddHistogram(classStr.Data(),"NGlobalTracksAway_RunNumber_vtxz","Global tracks (in |#eta|<0.9) vs Run number vs vtxz",kFALSE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_RunID_TProfile","Mean number of Away global tracks (in |#eta|<0.9) vs Run ID",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);
      man->AddHistogram(classStr.Data(),"NGlobalTracksAway_RunID_vtxz_TProfile","Mean number of Away global tracks (in |#eta|<0.9) vs Run ID vs vtxz",kTRUE,nRuns, 0, nRuns, AliReducedVarManager::kRunID,100,-10.,10.,AliReducedVarManager::kVtxZ, nMultBins, minMult, maxMult,AliReducedVarManager::kNGlobalTracksAway);

      continue;
    }*/

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
      man->AddHistogram(classStr.Data(), "Pt_TimeFromSOR", "<p_{T}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_TimeFromSOR", "<#eta> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi_TimeFromSOR", "<#varphi> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_TPCchi2", "DCAxy vs TPC chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAz_TPCchi2", "DCAz vs TPC chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_ITSchi2", "DCAxy vs ITS chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAz_ITSchi2", "DCAz vs ITS chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_goldenChi2", "DCAxy vs golden chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAz_goldenChi2", "DCAz vs golden chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(),"DCAxy_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(),"DCAz_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "NsigmaTPC", "Nsigma TPC", kFALSE, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSig);
      man->AddHistogram(classStr.Data(), "NsigmaTPC_P", "Nsigma TPC vs p", kFALSE, 50, -3.0, 3.0, AliReducedVarManager::kTPCnSig,50,1.,10., AliReducedVarManager::kP);

        
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
       man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTime);
      /* man->AddHistogram(classStr.Data(), "Mass_V0K0s", "Invariant mass, K^{0}_{s} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0);
       man->AddHistogram(classStr.Data(), "Mass_V0Lambda", "Invariant mass, #Lambda^{0} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+1);
       man->AddHistogram(classStr.Data(), "Mass_V0ALambda", "Invariant mass, #bar{#Lambda^{0}} assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+2);
       man->AddHistogram(classStr.Data(), "Mass_V0Gamma", "Invariant mass, #gamma conversion assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+3);
       */ man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
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
       
       continue;
    }   // end if for Pair classes of histograms 
    
    if(classStr.Contains("PairME")){
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE,500, 0.0, 5.0, AliReducedVarManager::kMass);
     // man->AddHistogram(classStr.Data(), "Multip", "Invariant mass vs multiplicity", kFALSE, 150, 0., 150., AliReducedVarManager::kNtracksSelected, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 300, -1.5, 1.5, AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Eta", "Eta", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 630, 0.0, 6.3, AliReducedVarManager::kPhi);
      continue;
    } 

    if(classStr.Contains("MCTRUTH")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassMC);
      man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTimeMC);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 500, 0.0, 30.0, AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPxMC);
      man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPyMC);
      man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPzMC);
      man->AddHistogram(classStr.Data(), "P", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPMC);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRapMC);
      man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEtaMC);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kThetaMC);
      man->AddHistogram(classStr.Data(), "PDG", "PDG value", kFALSE, 300, -600, 600, AliReducedVarManager::kPdgMC);
      man->AddHistogram(classStr.Data(), "PDGmother", "PDG value jpsi mother", kFALSE, 300, -6000, 6000, AliReducedVarManager::kPdgMC+1);
     /* man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Number of jpsi vs measured multiplicity", kFALSE,  nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks);
      man->AddHistogram(classStr.Data(), "multDist_evt_gen", "Number of jpsi vs true multiplicity", kFALSE,  nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
      man->AddHistogram(classStr.Data(), "multCorrel_evt", "True multiplicity vs measured multiplicity for Jpsi", kFALSE, nMultBins, minMult, maxMult, AliReducedVarManager::kNGlobalTracks, nMultBins, minMult, maxMult, AliReducedVarManager::kMCNch09);
      continue;*/
    }   // end if for Pair classes of histograms 
/*
    for(int cutMode=0; cutMode<4*nEstimators;cutMode++) {
      if(classStr.Contains(Form("pp_13TeV_Data_cutMode_%d",100+cutMode))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4);
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
      }
      if(classStr.Contains(Form("pp_13TeV_MC_cutMode_%d",100+cutMode))) {
        man->AddHistClass(classStr.Data());
        cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "multDist_evt_meas", "Measured multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4);
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas", "Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //2D unfolding
        man->AddHistogram(classStr.Data(), "multDist_evt_gen", "True multiplicity distribution", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09);
        man->AddHistogram(classStr.Data(), "multDist_evt_gen_trig", "True multiplicity distribution after trigger", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09+1);
        man->AddHistogram(classStr.Data(), "multCorrel_evt", "True multiplicity vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09+2); 
        man->AddHistogram(classStr.Data(), "multCorrel_prim", "True Pt spectrum vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09); //2D unfolding - not used - filling is false
        man->AddHistogram(classStr.Data(), "ptCorrel_prim", "True Pt spectrum vs measured pt spectrum", kFALSE, 100, 0., 10., AliReducedVarManager::kPt, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen", "True Pt spectrum (all generated primaries) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_evtloss", "True Pt spectrum (generated primaries from rejected events) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_gen_notrig", "True Pt spectrum (generated primaries from untriggered events) vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_prim_meas", "True Pt spectrum vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPtMC); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_prim_meas", "True Pt spectrum primaries vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_sec_meas", "True Pt spectrum secondaries vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_meas_evtcont", "True Pt spectrum (from contaminated events) vs measured multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kNGlobalTracks+cutMode/4, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
        man->AddHistogram(classStr.Data(), "multPtSpec_trk_inter", "Pt measured spectrum vs true multiplicity", kFALSE,  nMultBinspp, minMultpp, maxMultpp, AliReducedVarManager::kMCNch09, 100, 0., 10., AliReducedVarManager::kPt); //pt bins should be changed
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

    }*/

  }  // end loop over histogram classes
}
