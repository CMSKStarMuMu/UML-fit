using namespace RooFit;
using namespace std;

double deltaR(double eta1, double phi1, double eta2, double phi2);

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

double PDGB0Mass = 5.2797;
double PDGJpsiMass = 3.0969;
double PDGPsiPrimeMass = 3.6861;
double PDGKstMass = 0.8956;

map<int,vector<double>> swapCut = {
  {6, {0.13, 0.15, 0.07, 0.43, -0.45, 2.0}},
  {7, {0.07, 0.13, 0.11, 0.45, -0.65, 1.9}},
  {8, {0.09, 0.13, 0.07, 0.23, -0.45, 1.0}} };

TCanvas* c [nBins];

void createDataset(int year, int q2Bin = -1, int data = 0, int XGBv = 8, bool plot = false)
{
  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( year<6 || year>8 ) return;

  string XGBstr = "";
  if (!data && (XGBv==8 || (XGBv>0 && XGBv<6))) XGBstr = Form("_XGBv%i",XGBv);

  bool isJpsi = false;
  bool isPsi  = false;
  bool isLMNR = false;
  
  if (q2Bin==4)      isJpsi = true;
  else if (q2Bin==6) isPsi = true;
  else isLMNR = true;

  string dataString = data ? "DATA" : "MC";

  // define angular variables and variable for PU-reweighting
  RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
  RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
  RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (ctK, ctL, phi);
  RooRealVar wei ("weight","weight",1);
  RooRealVar mass("mass","mass", 0,10);
  // random variable [0,1] to keep or reject the event when performing data-like stat. studies
  RooRealVar rand("rand", "rand", 0,1);
  TRandom rand_gen(1029);
  RooArgSet reco_vars (ctK, ctL, phi, mass, rand);
  if (data==0) reco_vars.add(wei);

  // flags to mark which q2 bins should be filled
  bool runBin [nBins];
  string shortString [nBins];
  string longString  [nBins];
  for (int i=0; i<nBins; ++i) {
    runBin [i] = false;
    if ( q2Bin!=-1 && q2Bin!=i ) continue;
    if ( !data && (q2Bin==-1 && ( i==4 || i==6 ) ) ) continue; // b4 and b6 use different MC samples, so don't overwrite them when running the non-resonant for all the bins
    runBin [i] = true;
    shortString [i] = Form("b%i",i);
    longString  [i] = Form("q2 bin %i",i);
  }

  // Load ntuples
  string year_str = Form("201%i", year);
  auto f_num = TFile::Open(Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/%s/201%i.root",data==0?Form("MC-%s%s",isJpsi?"Jpsi":(isPsi?"Psi":"LMNR"),XGBv>0?Form("-XGBv%i",XGBv):""):"data",year));
  auto t_num = (TTree*)f_num->Get("ntuple");
  // if (data==0 && XGBv>9) {
  //   t_num->AddFriend("wTree",Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-%s-TMVAv%i/201%i.root",isJpsi?"Jpsi":(isPsi?"Psi":"LMNR"),XGBv-10,year));
  //   XGBstr = Form("_TMVAv%i",XGBv-10);
  // } else if (data==0 && XGBv>5) {
  //   t_num->AddFriend("wTree",Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-%s-XGBv%i/201%i.root",isJpsi?"Jpsi":(isPsi?"Psi":"LMNR"),XGBv,year));
  //   XGBstr = Form("_XGBv%i",XGBv);
  // }

  // TChain* t_num = new TChain();
  // if (data==0 && isLMNR)
  //   t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/201%iMC_LMNR_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root/ntuple", year, year));
  // else if (data==0 && isJpsi) {
  //   if (XGBv==2) t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/reweightV2/201%iMC_JPSI_withMCw_v2_addxcutvariable.root/ntuple", year, year));
  //   else if (XGBv>2) t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/reweightV%i/201%iMC_JPSI_v%i_MCw_addxcutvariable.root/ntuple",year,XGBv,year,XGBv));
  //   else t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/201%iMC_JPSI_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root/ntuple", year, year));
  // }
  // else if (data==0 && isPsi)
  //   t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/201%iMC_PSI_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root/ntuple", year, year));
  // else
  //   t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/fixBkg/201%idata_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root/ntuple", year, year));
  int numEntries = t_num->GetEntries();
  cout << "Reading tree with nr. of events = " << numEntries << endl;
  int counter;

  t_num->SetBranchStatus("*",0);
  // Import branches from ntuples:
  // angular variables
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_num->SetBranchStatus("cos_theta_k",1);
  t_num->SetBranchStatus("cos_theta_l",1);
  t_num->SetBranchStatus("phi_kst_mumu",1);
  t_num->SetBranchAddress( "cos_theta_k"     , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l"     , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu"    , &recoPhi       );

  // dimuon mass variables
  double recoDimuMass;
  t_num->SetBranchStatus("mumuMass",1);
  t_num->SetBranchAddress( "mumuMass", &recoDimuMass );

  // B0 mass variable
  double recoB0Mass;
  t_num->SetBranchStatus("tagged_mass",1);
  t_num->SetBranchAddress( "tagged_mass", &recoB0Mass );

  int passB0Psi_lmnr, passB0Psi_jpsi, passB0Psi_psip;
  t_num->SetBranchStatus("passB0Psi_lmnr",1);
  t_num->SetBranchStatus("passB0Psi_jpsi",1);
  t_num->SetBranchStatus("passB0Psi_psip",1);
  t_num->SetBranchAddress( "passB0Psi_lmnr", &passB0Psi_lmnr );
  t_num->SetBranchAddress( "passB0Psi_jpsi", &passB0Psi_jpsi );
  t_num->SetBranchAddress( "passB0Psi_psip", &passB0Psi_psip );

  // Anti-swap cut
  double kstTrkmPt		     = 0;
  double kstTrkmEta		     = 0;
  double kstTrkmPhi		     = 0;
  double kstTrkpPt		     = 0;
  double kstTrkpEta		     = 0;
  double kstTrkpPhi		     = 0;
  double mumPt			     = 0;
  double mumEta			     = 0;
  double mumPhi			     = 0;
  double mupPt			     = 0;
  double mupEta			     = 0;
  double mupPhi			     = 0;
  double selected_swapped_b0_mass    = 0;
  double selected_swapped_kstar_mass = 0;
  double selected_swapped_mumu_mass  = 0;
  Long64_t charge_pf_swapped_track     = 0;

  t_num->SetBranchStatus("kstTrkmPt",1);
  t_num->SetBranchStatus("kstTrkmEta",1);
  t_num->SetBranchStatus("kstTrkmPhi",1);
  t_num->SetBranchStatus("kstTrkpPt",1);
  t_num->SetBranchStatus("kstTrkpEta",1);
  t_num->SetBranchStatus("kstTrkpPhi",1);
  t_num->SetBranchStatus("mumPt",1);
  t_num->SetBranchStatus("mumEta",1);
  t_num->SetBranchStatus("mumPhi",1);
  t_num->SetBranchStatus("mupPt",1);
  t_num->SetBranchStatus("mupEta",1);
  t_num->SetBranchStatus("mupPhi",1);
  t_num->SetBranchStatus("selected_swapped_b0_mass",1);
  t_num->SetBranchStatus("selected_swapped_kstar_mass",1);
  t_num->SetBranchStatus("selected_swapped_mumu_mass",1);
  t_num->SetBranchStatus("charge_pf_swapped_track",1);

  t_num->SetBranchAddress("kstTrkmPt"		    , &kstTrkmPt		   );
  t_num->SetBranchAddress("kstTrkmEta"		    , &kstTrkmEta		   );
  t_num->SetBranchAddress("kstTrkmPhi"		    , &kstTrkmPhi		   );
  t_num->SetBranchAddress("kstTrkpPt"		    , &kstTrkpPt		   );
  t_num->SetBranchAddress("kstTrkpEta"		    , &kstTrkpEta		   );
  t_num->SetBranchAddress("kstTrkpPhi"		    , &kstTrkpPhi		   );
  t_num->SetBranchAddress("mumPt"		    , &mumPt			   );
  t_num->SetBranchAddress("mumEta"       	    , &mumEta			   );
  t_num->SetBranchAddress("mumPhi"     		    , &mumPhi			   );
  t_num->SetBranchAddress("mupPt"		    , &mupPt			   );
  t_num->SetBranchAddress("mupEta"		    , &mupEta			   );
  t_num->SetBranchAddress("mupPhi"		    , &mupPhi			   );
  t_num->SetBranchAddress("selected_swapped_b0_mass"   , &selected_swapped_b0_mass   );
  t_num->SetBranchAddress("selected_swapped_kstar_mass", &selected_swapped_kstar_mass);
  t_num->SetBranchAddress("selected_swapped_mumu_mass" , &selected_swapped_mumu_mass );
  t_num->SetBranchAddress("charge_pf_swapped_track"    , &charge_pf_swapped_track    );

  // B0-kinematic variables
  // double recoB0pT, recoB0eta;
  // t_num->SetBranchAddress( "bPt"    , &recoB0pT  );
  // t_num->SetBranchAddress( "bEta"   , &recoB0eta );

  double tagB0;
  t_num->SetBranchStatus("tagB0",1);
  t_num->SetBranchAddress( "tagB0"    , &tagB0     );

  // event number for even/odd splitting
  Long64_t eventN;
  t_num->SetBranchStatus("eventN",1);
  t_num->SetBranchAddress( "eventN", &eventN     );

  int XCut = 0;
  if (isJpsi) {
    t_num->SetBranchStatus("xcut",1);
    t_num->SetBranchAddress( "xcut", &XCut );
  }

  int xBin;

  // event pileup weight
  if (data==0) {
    // flavour tagging variables
    double genSignal;
    t_num->SetBranchStatus("genSignal",1);
    t_num->SetBranchAddress( "genSignal", &genSignal );
  
    std::cout << "is MC"  << std::endl;
    float PUweight = 1;
    t_num->SetBranchStatus("weight",1);
    t_num->SetBranchAddress( "weight", &PUweight );

    // MC weights
    double XGBweight = 1;
    float fXGBweight = 1;
    if (XGBv==8) {
      t_num->SetBranchStatus("XGBv8",1);
      t_num->SetBranchAddress( "XGBv8", &fXGBweight );
    }
    else if (XGBv>9) {
      t_num->SetBranchStatus("MCw",1);
      t_num->SetBranchAddress( "MCw", &XGBweight );
    }
    else if (XGBv>5) {
      t_num->SetBranchStatus("MCw",1);
      t_num->SetBranchAddress( "MCw", &fXGBweight );
    }
    else if (XGBv>0) {
      t_num->SetBranchStatus("sf_to_data",1);
      t_num->SetBranchAddress( "sf_to_data", &fXGBweight );
    }

    // Define datasets for five efficiency terms
    RooDataSet* data_ctRECO_ev [nBins];
    RooDataSet* data_ctRECO_od [nBins];
    RooDataSet* data_wtRECO_ev [nBins];
    RooDataSet* data_wtRECO_od [nBins];
    for (int i=0; i<nBins; ++i) {
      if (runBin[i] ) {
      data_ctRECO_ev [i] = new RooDataSet( ("data_ctRECO_ev_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (even)",
					   reco_vars, "weight" );
      data_ctRECO_od [i] = new RooDataSet( ("data_ctRECO_od_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (odd)",
					   reco_vars, "weight" );
      data_wtRECO_ev [i] = new RooDataSet( ("data_wtRECO_ev_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (even)",
					   reco_vars, "weight" );
      data_wtRECO_od [i] = new RooDataSet( ("data_wtRECO_od_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (odd)",
					   reco_vars, "weight" );
      }
    }
    
    // Prepare numerator dataset
    cout<<"Starting numerator dataset filling..."<<endl;
    counter=0;
    for (int iCand=0; iCand<numEntries; ++iCand) {
      t_num->GetEntry(iCand);
      // anti-radiation cut
      if (isLMNR && passB0Psi_lmnr == 0) continue;
      else if (isJpsi && passB0Psi_jpsi == 0) continue;
      else if (isPsi  && passB0Psi_psip == 0)  continue;
  
      // find q2 bin of current candidate
      xBin=-1;
      for (int i=0; i<nBins; ++i)
        if ( runBin[i] )
  	if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
  	     ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
  	  xBin = i;
  	  break;
  	}
      if (xBin<0) continue;
      
      // apply cut for bin 4 
      if (isJpsi && XCut>0) continue;

      // anti-swap cut
      if (fabs(selected_swapped_mumu_mass-PDGJpsiMass)<swapCut[year][0] &&
	  fabs(selected_swapped_kstar_mass-PDGKstMass)<swapCut[year][1] &&
	  fabs(selected_swapped_b0_mass-PDGB0Mass-selected_swapped_mumu_mass+PDGJpsiMass)<swapCut[year][2] &&
	  ( ( charge_pf_swapped_track<0 && (deltaR(kstTrkmEta,kstTrkmPhi,mumEta,mumPhi)<swapCut[year][3] &&
					    (mumPt-kstTrkmPt)/kstTrkmPt>swapCut[year][4] &&
					    (mumPt-kstTrkmPt)/kstTrkmPt<swapCut[year][5]) ) ||
	    ( charge_pf_swapped_track>0 && (deltaR(kstTrkpEta,kstTrkpPhi,mupEta,mupPhi)<swapCut[year][3] &&
					    (mupPt-kstTrkpPt)/kstTrkpPt>swapCut[year][4] &&
					    (mupPt-kstTrkpPt)/kstTrkpPt<swapCut[year][5]) ) ) ) continue;

      // status display
      if ( iCand > 1.0*counter*numEntries/100 ) {
        cout<<counter<<"%"<<endl;
        counter += 10;
      }
      //generate random variable [0,1]
      // fill dataset
      ctK.setVal(recoCosThetaK);
      ctL.setVal(recoCosThetaL);
      phi.setVal(recoPhi);
      mass.setVal(recoB0Mass);
      rand.setVal(rand_gen.Uniform(1));
      if (genSignal != tagB0+1) { // correctly tagged events
	if (eventN%2==0) data_ctRECO_ev[xBin]->add(reco_vars,PUweight*XGBweight*fXGBweight);
	else data_ctRECO_od[xBin]->add(reco_vars,PUweight*XGBweight*fXGBweight);
      } else { // wrongly tagged events
	if (eventN%2==0) data_wtRECO_ev[xBin]->add(reco_vars,PUweight*XGBweight*fXGBweight);
	else data_wtRECO_od[xBin]->add(reco_vars,PUweight*XGBweight*fXGBweight);
      }
    }
    cout<<"Dataset prepared"<<endl;

    // Save datasets in workspaces
    RooWorkspace *ws_ev [nBins];
    RooWorkspace *ws_od [nBins];
    for (int i=0; i<nBins; ++i) if (runBin[i]) {
      // Skip the creation of a file when the correct-tag efficiency cannot be computed (empty numerators)
      // which usually means that either you are using a resonant MC, which does not fill signal q2 bins,
      // or using a bin too fine, or out of range
      // If correct-tag numerator is filled and wrong-tag is not, a warning is returned
      if ( data_ctRECO_ev[i]->numEntries()==0 || data_ctRECO_od[i]->numEntries()==0 ) {
  	cout<<"Error: ctRECO is empty in q2 bin "<<i<<endl;
  	continue;
      }
      if ( data_wtRECO_ev[i]->numEntries()==0 || data_wtRECO_od[i]->numEntries()==0 ) {
  	cout<<"Warning: wtRECO is empty in q2 bin "<<i<<endl;
      }
      ws_ev[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin even datasets");
      ws_od[i] = new RooWorkspace(("ws_"+shortString[i]+"p1").c_str(),"Workspace with single-bin odd datasets");
      ws_ev[i]->import( *data_ctRECO_ev[i] );
      ws_od[i]->import( *data_ctRECO_od[i] );
      ws_ev[i]->import( *data_wtRECO_ev[i] );
      ws_od[i]->import( *data_wtRECO_od[i] );
      TFile* fout = new TFile( ( "reco"+dataString+"Dataset_"+shortString[i]+ "_" + year_str + XGBstr + ".root" ).c_str(), "RECREATE" );
      ws_ev[i]->Write();
      ws_od[i]->Write();
      fout->Close();
    }
  
    // Plot 1D distributions of datasets
    if (plot ) {
  
      // to keep all distributions visible in the same plot, the ones with higher stats (tipically denominators) need to be rescaled
      double rescFac1 = 1.0/12;
      double rescFac2 = 1.0;
      double rescFac3 = 1.0/25;
  
      TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);
  
      RooPlot* xframe_ev [nBins];
      RooPlot* yframe_ev [nBins];
      RooPlot* zframe_ev [nBins];
      RooPlot* xframe_od [nBins];
      RooPlot* yframe_od [nBins];
      RooPlot* zframe_od [nBins];
  
      bool legFilled = false;
  
      RooDataSet* data_RECO_ev [nBins];
      RooDataSet* data_RECO_od [nBins];
  
      for (int i=0; i<nBins; ++i) if (runBin[i]) {
  	// Create dataset containing both correct-tag and wrong-tag events
  	data_RECO_ev [i] = new RooDataSet( ("data_RECO_ev_"+shortString[i]).c_str(), "Reconstructed candidates after selections (even)", data_ctRECO_ev[i], vars );
  	data_RECO_od [i] = new RooDataSet( ("data_RECO_od_"+shortString[i]).c_str(), "Reconstructed candidates after selections (odd)" , data_ctRECO_od[i], vars );
  	data_RECO_ev[i]->append(*(data_wtRECO_ev[i]));
  	data_RECO_od[i]->append(*(data_wtRECO_od[i]));
  
  	// create frames (one for each bin/parity/variable, but all the six efficiency terms are plotted together)
  	c [i] = new TCanvas(("c_"+shortString[i]).c_str(),("Num and Den 1D projections - "+longString[i]).c_str(),2000,1400);
  	xframe_ev [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (even)").c_str()));
  	yframe_ev [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (even)").c_str()));
  	zframe_ev [i] = phi.frame(Title((longString[i]+" #phi distributions (even)").c_str()));
  	xframe_od [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (odd)").c_str()));
  	yframe_od [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (odd)").c_str()));
  	zframe_od [i] = phi.frame(Title((longString[i]+" #phi distributions (odd)").c_str()));
  
  	// plot datasets on frames
  	if (!legFilled) { // the first time assign names to tag them in the legend
  	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40),Name("plCTrecoNum"));
  	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40),Name("plWTrecoNum"));
  	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40),Name("plRecoNum"));
  	} else {
  	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  	}
  	data_ctRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	data_wtRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	data_RECO_od  [i]->plotOn(xframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  	data_ctRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	data_wtRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	data_RECO_ev  [i]->plotOn(yframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  	data_ctRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	data_wtRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	data_RECO_od  [i]->plotOn(yframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  	data_ctRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	data_wtRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	data_RECO_ev  [i]->plotOn(zframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  	data_ctRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
  	data_wtRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
  	data_RECO_od  [i]->plotOn(zframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
  
  	xframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
  	yframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
  	zframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
  	xframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
  	yframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
  	zframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
  	xframe_ev[i]->SetMaximum(xframe_ev[i]->GetMaximum()*rescFac1);
  	yframe_ev[i]->SetMaximum(yframe_ev[i]->GetMaximum()*rescFac1);
  	zframe_ev[i]->SetMaximum(zframe_ev[i]->GetMaximum()*rescFac1);
  	xframe_od[i]->SetMaximum(xframe_od[i]->GetMaximum()*rescFac1);
  	yframe_od[i]->SetMaximum(yframe_od[i]->GetMaximum()*rescFac1);
  	zframe_od[i]->SetMaximum(zframe_od[i]->GetMaximum()*rescFac1);
  
  	// Fill legend (only one time)
  	if (!legFilled) {
  	  string strRescFac1 = (rescFac1<1?Form(" (*%1.2f)",rescFac1):"");
  	  string strRescFac2 = (rescFac2<1?Form(" (*%1.2f)",rescFac2):"");
  	  string strRescFac3 = (rescFac3<1?Form(" (*%1.2f)",rescFac3):"");
  	  leg->AddEntry(xframe_ev[i]->findObject("plRecoNum"  ),"Post-selection distribution","lep");
  	  leg->AddEntry(xframe_ev[i]->findObject("plCTrecoNum"),"Post-selection correct-tag distribution","lep");
  	  leg->AddEntry(xframe_ev[i]->findObject("plWTrecoNum"),"Post-selection wrong-tag distribution","lep");
  	  legFilled = true;
  	}
  
  	// Plot even distributions in the top row and odd ones in the bottom row
  	c[i]->Divide(3,2);
  	c[i]->cd(1);
  	gPad->SetLeftMargin(0.17); 
  	xframe_ev[i]->Draw();
  	leg->Draw("same");
  	c[i]->cd(2);
  	gPad->SetLeftMargin(0.17); 
  	yframe_ev[i]->Draw();
  	leg->Draw("same");
  	c[i]->cd(3);
  	gPad->SetLeftMargin(0.17); 
  	zframe_ev[i]->Draw();
  	leg->Draw("same");
  	c[i]->cd(4);
  	gPad->SetLeftMargin(0.17); 
  	xframe_od[i]->Draw();
  	leg->Draw("same");
  	c[i]->cd(5);
  	gPad->SetLeftMargin(0.17); 
  	yframe_od[i]->Draw();
  	leg->Draw("same");
  	c[i]->cd(6);
  	gPad->SetLeftMargin(0.17); 
  	zframe_od[i]->Draw();
  	leg->Draw("same");
  
  	c[i]->SaveAs( ("plotDist_d/dist_RECO_"+dataString+"_"+shortString[i]+".pdf").c_str() );
        }
    }
  
  }
  
  
  else{
    RooDataSet* data [nBins];
    for (int i=0; i<nBins; ++i) {
      if (runBin[i] ){
        data [i] = new RooDataSet( ("data_"+shortString[i]).c_str(), "Reconstructed candidates after selections",
					   reco_vars);
      }
    }
  
    // Prepare numerator dataset
    cout<<"Starting numerator dataset filling..."<<endl;
    counter=0;
    for (int iCand=0; iCand<numEntries; ++iCand) {
      t_num->GetEntry(iCand);
      // anti-radiation cut
      if (isLMNR && passB0Psi_lmnr == 0) continue;
      else if (isJpsi && passB0Psi_jpsi == 0) continue;
      else if (isPsi  && passB0Psi_psip == 0)  continue;
  
      // find q2 bin of current candidate
      xBin=-1;
      for (int i=0; i<nBins; ++i)
        if ( runBin[i] )
  	if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
  	     ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
  	  xBin = i;
  	  break;
  	}
      if (xBin<0) continue;

      // apply cut for bin 4 
      if (isJpsi && XCut>0) continue;

      // anti-swap cut
      if (fabs(selected_swapped_mumu_mass-PDGJpsiMass)<swapCut[year][0] &&
	  fabs(selected_swapped_kstar_mass-PDGKstMass)<swapCut[year][1] &&
	  fabs(selected_swapped_b0_mass-PDGB0Mass-selected_swapped_mumu_mass+PDGJpsiMass)<swapCut[year][2] &&
	  ( ( charge_pf_swapped_track<0 && (deltaR(kstTrkmEta,kstTrkmPhi,mumEta,mumPhi)<swapCut[year][3] &&
					    (mumPt-kstTrkmPt)/kstTrkmPt>swapCut[year][4] &&
					    (mumPt-kstTrkmPt)/kstTrkmPt<swapCut[year][5]) ) ||
	    ( charge_pf_swapped_track>0 && (deltaR(kstTrkpEta,kstTrkpPhi,mupEta,mupPhi)<swapCut[year][3] &&
					    (mupPt-kstTrkpPt)/kstTrkpPt>swapCut[year][4] &&
					    (mupPt-kstTrkpPt)/kstTrkpPt<swapCut[year][5]) ) ) ) continue;

      // status display
      if ( iCand > 1.0*counter*numEntries/100 ) {
        cout<<counter<<"%"<<endl;
        counter += 10;
      }
      // fill dataset
      ctK.setVal(recoCosThetaK);
      ctL.setVal(recoCosThetaL);
      phi.setVal(recoPhi);
      mass.setVal(recoB0Mass);
      rand.setVal(rand_gen.Uniform(1));
      data[xBin]->add(reco_vars);
    }
    cout<<"Dataset prepared"<<endl;

    // Save datasets in workspaces
    RooWorkspace *ws [nBins];
    for (int i=0; i<nBins; ++i) if (runBin[i]) {
      if ( data[i]->numEntries()==0  ) {
        cout<<"Error: RECO data is empty in q2 bin "<<i<<endl;
  	continue;
      }
      ws[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin data datasets");
      ws[i]->import( *data[i] );
      TFile* fout = new TFile( ( "reco"+dataString+"Dataset_"+shortString[i]+ "_" + year_str + XGBstr + ".root" ).c_str(), "RECREATE" );
      ws[i]->Write();
      fout->Close();
    }

  }


}

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  float deltaPhi = phi1 - phi2;
  if (fabs(deltaPhi)>TMath::Pi()) {
    if (deltaPhi>0) deltaPhi = deltaPhi - 2*TMath::Pi();
    else deltaPhi = deltaPhi + 2*TMath::Pi();
  }
  float deltaEta = eta1 - eta2;
  return sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
}
