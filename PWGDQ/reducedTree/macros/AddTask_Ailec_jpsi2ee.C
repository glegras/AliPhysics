#include "AliReducedAnalysisJpsi2ee.cxx"

 void HistogramManager(AliReducedAnalysisJpsi2ee* task, AliHistogramManager* man, TString prod) ;
 void MixingHandler(AliReducedAnalysisJpsi2ee* task, AliMixingHandler* handler) ;
 void DefineHistograms(AliReducedAnalysisJpsi2ee* task, TString prod) ;
 void Setup(AliReducedAnalysisJpsi2ee* processor,  TString prod) ;

AliAnalysisTask* AddTask_Ailec_jpsi2ee( Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){

 AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
 jpsi2eeAnalysis->Init();
 Setup(jpsi2eeAnalysis, prod);
 AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager");
 task->AddTask(jpsi2eeAnalysis);

 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

 AliAnalysisDataContainer* cReducedEvent = NULL;

 cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");

 mgr->AddTask(task);

 mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());

 AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("Histos", THashList::Class(),
                AliAnalysisManager::kOutputContainer, "AnalysisHistograms_AilecTask.root");

 mgr->ConnectOutput(task, 1, cOutputHist );

 return task;
}
 
 void Setup(AliReducedAnalysisJpsi2ee* processor,  TString prod) {
	 	 
  processor->SetRunEventMixing(kTRUE);
  processor->SetRunPairing(kTRUE);
  //processor->SetRunOverMC(kFALSE);
  processor->SetRunLikeSignPairing(kTRUE);
  
  // Set event cuts
  AliReducedEventCut* evCut2 = new AliReducedEventCut("VertexZ","Vertex selection");
  evCut2->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut2->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  processor->AddEventCut(evCut2);
    
  //Set track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kP, 1,30.0);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  standardCut->SetRequestSPDany();
  processor->AddTrackCut(standardCut);
  
   // set track prefilter cuts
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.7,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // set pair prefilter cuts
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddPrefilterPairCut(prefPairCut);
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut1);
  
  MixingHandler(processor, processor->GetMixingHandler());

  HistogramManager(processor, processor->GetHistogramManager(), prod);
   
}

void HistogramManager(AliReducedAnalysisJpsi2ee* task, AliHistogramManager* man, TString prod) {
  
  // setup the histograms manager
  
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod);
  
  AliReducedVarManager::SetUseVars(man->GetUsedVars());
}


void MixingHandler(AliReducedAnalysisJpsi2ee* task, AliMixingHandler* handler) {
   //
   // setup the mixing handler
   handler = task->GetMixingHandler();  
   handler->SetPoolDepth(50);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1.0);
   handler->SetDownscaleTracks(1);
   handler->SetMixLikeSign(kTRUE);
   const Int_t nZbinLimits = 2;
   Float_t zLims[nZbinLimits] = {-10., /*  -8.,  -6.,  -4.,  -2., 0.,  2.,  4.,  6.,  8., */10. };
   //r(Int_t i=0;i<=10;++i) zLims[nZbinLimits] = i
   
 handler->AddMixingVariable(AliReducedVarManager::kVtxZ, nZbinLimits, zLims);   
}

void DefineHistograms(AliReducedAnalysisJpsi2ee* task, TString prod){
	

 
 AliHistogramManager* man = task->GetHistogramManager(); 
 
 TString histClasses = "";	
 histClasses += "Event_AfterCuts;";

  for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      
      histClasses += Form("Track_%s;", cutName.Data());
      //histClasses += Form("PairTRPM_%s;", cutName.Data());
      histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
      histClasses += Form("PairMEPM_%s;PairMEPP_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
     
  }
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};
  runNBins = 1140;
  runHistRange[0] = 252200.;
  runHistRange[1] = 252400.;
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  int kNMassBins = 500;
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
     // man->AddHistogram(classStr.Data(),"kNtracks2Selected","",  kTRUE, 150, 0., 150., AliReducedVarManager::kNtracks2Selected);
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,3000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksTotal_BeamIntensity0_prof","Number of total tracks per event",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,20000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,300,0.,300.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,90, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors_BeamIntensity0","",kTRUE, 100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
    
  }  
    
    // Histograms for pairs
    if(classStr.Contains("PairME")){
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE,kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
     // man->AddHistogram(classStr.Data(), "Multip", "Invariant mass vs multiplicity", kFALSE, 150, 0., 150., AliReducedVarManager::kNtracksSelected, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 300, -1.5, 1.5, AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Eta", "Eta", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 630, 0.0, 6.3, AliReducedVarManager::kPhi);
     }   
      
      if(classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE")) {
	  man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
	      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Multip", "Invariant mass vs multiplicity", kFALSE, 150, 0., 150., AliReducedVarManager::kNtracksSelected, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass); 
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR", "Invariant mass vs time from SOR", kFALSE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR_prof", "<Invariant mass> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_Cent_TimeFromSOR_prof", "<Invariant mass> vs centrality and time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 18, 0.0, 90., AliReducedVarManager::kCentVZERO, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "Pt_coarse", "", kFALSE, 20, 0.0, 20.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "DCAxy", "", kFALSE, 200, -2.0, 2.0, AliReducedVarManager::kDcaXY);
        man->AddHistogram(classStr.Data(), "DCAz", "", kFALSE, 200, -2.0, 2.0, AliReducedVarManager::kDcaZ);
        man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTime);
        //man->AddHistogram(classStr.Data(), "PseusoproperTimeError", "Error on Pseudoproper decay time", kFALSE,100, 0, 0.1, AliReducedVarManager::kPseudoProperTimeError);
	      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
	      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
        man->AddHistogram(classStr.Data(), "Leg1TPCchi2_Leg2TPCchi2", "", kFALSE, 100, 0.0, 6.0, AliReducedVarManager::kPairLegTPCchi2, 100, 0.0, 6.0, AliReducedVarManager::kPairLegTPCchi2+1);
        man->AddHistogram(classStr.Data(), "Leg1ITSchi2_Leg2ITSchi2", "", kFALSE, 100, 0.0, 6.0, AliReducedVarManager::kPairLegITSchi2, 100, 0.0, 6.0, AliReducedVarManager::kPairLegITSchi2+1);
        man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
        man->AddHistogram(classStr.Data(), "CosThetaStarCS_ptMC", "cos(#theta^{*})_{CS} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 50, 0., 1., AliReducedVarManager::kPtMC);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE_pt", "cos(#theta^{*})_{HE} vs pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 100, 0., 1., AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "CosThetaStarHE_pt_coarse", "cos(#theta^{*})_{HE} vs pt", kFALSE, 10, -1.0, 1.0, AliReducedVarManager::kPairThetaHE, 20, 0., 20., AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "PhiStarCS", "#varphi^{*}_{CS}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiCS);
        man->AddHistogram(classStr.Data(), "PhiStarHE", "#varphi^{*}_{HE}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiHE);         
        Double_t v2PtLims[9] = {0.,0.15,0.3,0.5,1.0,3.0,5.0,7.0,10.0};
        Double_t v2MassLims[8] = {2.4,2.6,2.8,2.92,3.04,3.16,3.50,4.0};
        Double_t v2CentLims[4] = {30.,50.,70.,90.};
        //man->AddHistogram(classStr.Data(), "v2VZEROA_massPtCent", "v_{2}^{VZERO-A} (mass,p_{T})", kTRUE, 7, v2MassLims, AliReducedVarManager::kMass, 8, v2PtLims, AliReducedVarManager::kPt, 3, v2CentLims, AliReducedVarManager::kCentVZERO, "", "", "", AliReducedVarManager::kVZEROFlowVn+0*6+1);
        //man->AddHistogram(classStr.Data(), "v2VZEROC_massPtCent", "v_{2}^{VZERO-C} (mass,p_{T})", kTRUE, 7, v2MassLims, AliReducedVarManager::kMass, 8, v2PtLims, AliReducedVarManager::kPt, 3, v2CentLims, AliReducedVarManager::kCentVZERO, "", "", "", AliReducedVarManager::kVZEROFlowVn+1*6+1);
      }   // end if "QA"
   
  
   
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
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -2.0, 2.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -2.0, 2.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(),"DCAxy_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(),"DCAz_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaZ); 
     // man->AddHistogram(classStr.Data(),"kNtracks2Selected","",  kTRUE, 150, 0., 150., AliReducedVarManager::kNtracks2Selected);

}
}
}

