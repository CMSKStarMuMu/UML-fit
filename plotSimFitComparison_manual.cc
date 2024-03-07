#include "interface/utils.h"

string year[3] = {"2016","2017","2018"};

void plotSimFitComparison_manual(int q2Bin = 4)
{
  gROOT->SetBatch(true);
  
  bool verbose = false;

  string shortString = Form("b%ip0",q2Bin);

  static const int nPars = 8;
  string parName[nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

  static const int nMassPars = 2;
  string massParName[nMassPars] = {"mean_{RT}^{%s}",
				   "fsig_"+shortString+"_%s"};
				   // "slope^{%s}"};

  static const int nMassConstPars = 14;
  string massConstParName[nMassConstPars] = {"f_{M}^{%s}",
					     "#alpha_{RT1}^{%s}",
					     "#alpha_{RT2}^{%s}",
					     "#sigma_{RT1}^{%s}",
					     "#sigma_{RT2}^{%s}",
					     "f^{RT%s}",
					     "n_{RT1}^{%s}",
					     "n_{RT2}^{%s}",
					     "#alpha_{WT1}^{%s}",
					     "#alpha_{WT2}^{%s}",
					     "#sigma_{WT1}^{%s}",
					     "deltaPeakVar^{%s}",
					     "n_{WT1}^{%s}",
					     "n_{WT2}^{%s}"};

  static const int nSPars = 6;
  string sParName[nSPars] = {"Fs","As","A4s","A5s","A7s","A8s"};

  double PDG_Fl  [2] = {0.571,0.463};
  double PDG_Fleh[2] = {0.007,0.028};
  double PDG_Flel[2] = {0.007,0.040};

  double bestParSim[nPars];
  double lowParSim [nPars];
  double highParSim[nPars];
  double bestPar[3][nPars];
  double lowPar [3][nPars];
  double highPar[3][nPars];
  
  double bestSParSim[nSPars];
  double lowSParSim [nSPars];
  double highSParSim[nSPars];
  double bestSPar[3][nSPars];
  double lowSPar [3][nSPars];
  double highSPar[3][nSPars];
  
  double bestMassParSim[3][nMassPars];
  double errMassParSim [3][nMassPars];
  double bestMassConstParSim[3][nMassConstPars];
  double errMassConstParSim [3][nMassConstPars];
  double bestMassPar[3][nMassPars];
  double errMassPar [3][nMassPars];
  double bestMassConstPar[3][nMassConstPars];
  double errMassConstPar [3][nMassConstPars];

  double massConstVal[3][nMassConstPars];
  double massConstErr[3][nMassConstPars];


  // First reads values of parameters which are not P_i or F_L from the log file of the simultaneous fit
  // for AN v9
  // string finLogName = "/afs/cern.ch/work/d/dini/public/ResonantXGBv8/Corr_frac_sigma/Bin%i-allYears-Bern.log";
  // for WP 90
  string finLogName = "/afs/cern.ch/work/d/dini/public/ResonantSwapCut/Bin4SimBern/logs_simFit4d/simfit_data_fullAngularMass_Swave_%i_0_0_0_0_1_8_0_1_2_2016_2017_2018_unbl4.log";
//   string finLogName = "/afs/cern.ch/work/d/dini/public/ResonantSwapCut/Bin6SimNoZed/logs_simFit4d/simfit_data_fullAngularMass_Swave_%i_0_0_0_0_0_8_0_1_2_2016_2017_2018_unbl4.log";
  // string finLogName = "logs_parSub/simfit_data_fullAngularMass_Swave_0_%i_0.out";
  // string finLogName = "logs_simFit4d/simfit_data_fullAngularMass_Swave_%i_1_0_0_0_0_1_2016_2017_2018.log";
  ifstream fin_log(Form(finLogName.c_str(),q2Bin));
  if ( !fin_log ) {
    cout<<"Log file is problematic!"<<endl;
    return;
  }
  string fin_line;
  while (getline(fin_log, fin_line)) {
    stringstream ss(fin_line);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> vec(begin, end);
    if (vec.size()<1) continue;
    // cout<<vec[0]<<endl;
    for (int iPar = 0; iPar<nSPars; ++iPar)
      if ( vec[0].compare(sParName[iPar]) == 0 ) {
	bestSParSim[iPar] = atof(vec[2].c_str());
	highSParSim[iPar] = atof(vec[4].c_str());
	lowSParSim [iPar] = -1. * atof(vec[4].c_str());
      }
    for (int iPar = 0; iPar<3*nMassPars; ++iPar)
      if ( vec[0].compare(Form(massParName[iPar/3].c_str(),year[iPar%3].c_str())) == 0 ) {
	bestMassParSim[iPar%3][iPar/3] = atof(vec[2].c_str()); // [year][parIndex]
	errMassParSim [iPar%3][iPar/3] = atof(vec[4].c_str());
	if (verbose) std::cout << "sim results (free): \t" << Form(massParName[iPar/3].c_str(),year[iPar%3].c_str()) << " = " << bestMassParSim[iPar%3][iPar/3] << std::endl;
      }
    for (int iPar = 0; iPar<3*nMassConstPars; ++iPar)
      if ( vec[0].compare(Form(massConstParName[iPar/3].c_str(),year[iPar%3].c_str())) == 0 ) {
	bestMassConstParSim[iPar%3][iPar/3] = atof(vec[2].c_str());
	errMassConstParSim [iPar%3][iPar/3] = atof(vec[4].c_str());
	if (verbose) std::cout << "sim results (constr): \t" << Form(massConstParName[iPar/3].c_str(),year[iPar%3].c_str()) << " = " << bestMassConstParSim[iPar%3][iPar/3] << std::endl;
      }
  }
  fin_log.close();

  // INPUTS
  // Now retrieves the angular parameters from the root file (they are saved in the fitResultsTree)
  string finName = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/simFitResults/simFitResult_data_fullAngularMass_Swave_%s_b%ip1_XGBv8.root";
  /// tmp to change back
  string finNameSingleYear = "/afs/cern.ch/work/d/dini/public/ResonantSwapCut/Bin4SingleBern/simFitResults4d/simFitResult_data_fullAngularMass_Swave_%s_b%i-XGBv8_unbl4.root";
  // which single year fits are available
  std::vector<int> haveSingleFits = {2016,2017,2018};
  int nYears = (int)haveSingleFits.size();

  // Here the files used to get the constraints
  string finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-LMNR-XGBv8/%s.root";
  if (q2Bin==4) finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-Jpsi-XGBv8/%s.root";
  else if (q2Bin==6) finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-Psi-XGBv8/%s.root";


  // Get the angular pars from the root file 
  auto finSim = TFile::Open(Form(finName.c_str(),"201620172018",q2Bin));
  if ( !finSim || finSim->IsZombie() ) {
    cout<<"Sim file is problematic!"<<endl;
    return;
  }

  auto tinSim = (TTree*)finSim->Get("fitResultsTree");
  if ( !tinSim || tinSim->IsZombie() || tinSim->GetEntries() < 1 ) {
    cout<<"Sim tree is problematic!"<<endl;
    finSim->Close();
    return;
  }
  for (int iPar=0; iPar<nPars; ++iPar) {
    tinSim->SetBranchAddress((parName[iPar]+"_best").c_str(),&bestParSim[iPar]);
    tinSim->SetBranchAddress((parName[iPar]+"_low" ).c_str(),&lowParSim [iPar]);
    tinSim->SetBranchAddress((parName[iPar]+"_high").c_str(),&highParSim[iPar]);
  }
  tinSim->GetEntry(0);

  // Get pars from the single year fits, if any 
  for (int iYear=0; iYear<nYears; ++iYear) {
    // temp to change back
    finName = finNameSingleYear;
//     auto fin = TFile::Open(Form(finNameSingleYear.c_str(),year[iYear].c_str(),q2Bin));
    auto fin = TFile::Open(Form(finName.c_str(),year[iYear].c_str(),q2Bin));
    if ( !fin || fin->IsZombie() ) {
      cout<<year[iYear]<<" file is problematic!"<<endl;
      return;
    }
    auto tin= (TTree*)fin->Get("fitResultsTree");
    if ( !tin || tin->IsZombie() || tin->GetEntries() < 1 ) {
      cout<<year[iYear]<<" tree is problematic!"<<endl;
      fin->Close();
      return;
    }
    // Get P-wave angular pars 
    for (int iPar=0; iPar<nPars; ++iPar) {
      tin->SetBranchAddress((parName[iPar]+"_best").c_str(),&bestPar[iYear][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_low" ).c_str(),&lowPar [iYear][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_high").c_str(),&highPar[iYear][iPar]);
    }
    tin->GetEntry(0);

    // Get S-wave and mass pars 
    auto ws = (RooWorkspace*)fin->Get("wsp_out");
    if ( !ws || ws->IsZombie() ) {
      cout<<year[iYear]<<" workspace is problematic!"<<endl;
      return;
    }

    for (int iPar=0; iPar<nSPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(sParName[iPar].c_str());
      bestSPar[iYear][iPar] = par->getValV();
      lowSPar [iYear][iPar] = par->getErrorLo();
      highSPar[iYear][iPar] = par->getErrorHi();
    }
    for (int iPar=0; iPar<nMassPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massParName[iPar].c_str(),year[iYear].c_str()));
      if (!par) {
	cout<<"Error: "<<massParName[iPar]<<" not found in file "<<Form(finName.c_str(),year[iYear].c_str(),q2Bin)<<endl;
	return;
      }
      bestMassPar[iYear][iPar] = par->getValV();
      errMassPar [iYear][iPar] = par->getError();
      if (verbose) std::cout << "single fit (free): " << Form(massParName[iPar].c_str(),year[iYear].c_str()) << " = " << bestMassPar[iYear][iPar]  << std::endl;
    }
    for (int iPar=0; iPar<nMassConstPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massConstParName[iPar].c_str(),year[iYear].c_str()));
      if (!par) {
	cout<<"Error: "<<Form(massConstParName[iPar].c_str(),year[iYear].c_str()) <<" not found in file "<<Form(finName.c_str(),year[iYear].c_str(),q2Bin)<<endl;
	return;
      }
      bestMassConstPar[iYear][iPar] = par->getValV();
      errMassConstPar [iYear][iPar] = par->getError();
      if (verbose) std::cout << "single fit (constr): " << Form(massConstParName[iPar].c_str(),year[iYear].c_str()) << " = " << bestMassConstPar[iYear][iPar]  << std::endl;
    }

    fin->Close();
  }

  // Now read constraints from the mass-pdf files
  for (int iYear=0; iYear<3; ++iYear) {
    auto fin2 = TFile::Open(Form(finName2.c_str(),year[iYear].c_str()));
    if ( !fin2 || fin2->IsZombie() ) {
      cout<<year[iYear]<<" file2 is problematic!"<<endl;
      return;
    }
    auto ws2 = (RooWorkspace*)fin2->Get("w");
    if ( !ws2 || ws2->IsZombie() ) {
      cout<<year[iYear]<<" workspace2 is problematic!"<<endl;
      return;
    }
    ws2->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    double mean_rt_val = ws2->var(Form("mean_{RT}^{%i}", q2Bin))->getVal(); 
    double mean_rt_err = ws2->var(Form("mean_{RT}^{%i}", q2Bin))->getError(); 
    for (int iPar=0; iPar<nMassConstPars; ++iPar) {
      if (iPar==8) ws2->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
      if (iPar==11) {
	double mean_wt_val = ws2->var(Form("mean_{WT}^{%i}", q2Bin))->getVal(); 
	double mean_wt_err = ws2->var(Form("mean_{WT}^{%i}", q2Bin))->getError(); 
	massConstVal[iYear][iPar] = mean_rt_val-mean_wt_val;
	massConstErr[iYear][iPar] = sqrt(mean_rt_err*mean_rt_err+mean_wt_err*mean_wt_err); 
      } else if (iPar>0) {
	auto par2 = (RooRealVar*)ws2->var(Form(massConstParName[iPar].c_str(),Form("%i",q2Bin)));
	if ( !par2 || par2->IsZombie() ) {
	  cout<<Form(massConstParName[iPar].c_str(),year[iYear].c_str())<<" constraint not found"<<endl;
	  return;
	}
	massConstVal[iYear][iPar] = par2->getValV();
	massConstErr[iYear][iPar] = par2->getError();
      } else {
	double nrt_mc = ws2->var(Form("nRT_%i",q2Bin))->getVal();
	double nwt_mc = ws2->var(Form("nWT_%i",q2Bin))->getVal();
	massConstVal[iYear][iPar] = 1.;
        double fraction = nwt_mc / (nrt_mc + nwt_mc);
	massConstErr[iYear][iPar] = fM_sigmas[atoi(year[iYear].c_str())][q2Bin] / fraction/(1-fraction); 
      }
      cout<< "constr val " << Form(massConstParName[iPar].c_str(),year[iYear].c_str()) << " = " << massConstVal[iYear][iPar]<<endl;
      
    }
    fin2->Close();
  }


  // Plot P-wave parameters
  double x[nYears][nPars];
  double xe[nPars];
  double xSim[nPars];
  double xeSim[nPars];

  double ySim [nPars];
  double ylSim[nPars];
  double yhSim[nPars];
  double y [nYears][nPars];
  double yl[nYears][nPars];
  double yh[nYears][nPars];

  for (int iPar=0; iPar<nPars; ++iPar) {
    xSim[iPar]  = 0.5+iPar;
    xeSim[iPar] = 0.5;
    xe[iPar]    = 1.0/(2*nYears);
    ySim [iPar] = 0;
    ylSim[iPar] = -1 * lowParSim [iPar];
    yhSim[iPar] = highParSim[iPar];
    for (int iYear=0; iYear < nYears; ++iYear) {
      x [iYear][iPar] = (2.0*iYear+1)/(2*nYears)+iPar;
      y [iYear][iPar] = bestPar[iYear][iPar] - bestParSim[iPar];
      yl[iYear][iPar] = -1 * lowPar [iYear][iPar];
      yh[iYear][iPar] = highPar[iYear][iPar];
    }
  }

  TCanvas canv("canv","canv",1500,1000);
  canv.cd();

  gStyle->SetOptStat(0);
  // cCover.cd()->SetTicky(2);
  string plotTitle = Form("Comparison of simultaneous fit results - q2 bin %i;;(single-year result) - (simultaneous result)",q2Bin);
  auto hLab = new TH1S ("hLab",plotTitle.c_str(),nPars,0,nPars);
  for (int iBin=0; iBin<nPars; ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,parName[iBin].c_str());
  double yRange = 0.1;
  if (q2Bin==4) yRange = 0.05;
  hLab->SetMinimum(-1*yRange);
  hLab->SetMaximum(yRange);
  // hLab->GetYaxis()->SetTickLength(0.006);
  hLab->Draw();

  auto grSim = new TGraphAsymmErrors(nPars,xSim,ySim,xeSim,xeSim,ylSim,yhSim);
  grSim->SetName("grSim");
  grSim->SetFillColor(30);
  grSim->Draw("2");

  TGraphAsymmErrors* grPDG;
  double yPDG[1];
  double ylPDG[1];
  double yhPDG[1];
  if (q2Bin==4 || q2Bin==6) {
    uint iCR = (q2Bin-4) / 2;
    yPDG [0] = PDG_Fl[iCR]-bestParSim[0];
    ylPDG[0] = PDG_Flel[iCR];
    yhPDG[0] = PDG_Fleh[iCR];
    grPDG = new TGraphAsymmErrors(1,xSim,yPDG,xeSim,xeSim,ylPDG,yhPDG);
    grPDG->SetName("grPDG");
    grPDG->SetLineColor(15);
    grPDG->SetLineWidth(2);
    grPDG->Draw("P");
  }

  TLine line (0,0,nPars,0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.SetLineColor(13);
  line.Draw();

  TGraphAsymmErrors* gr [nYears];
  int colours [3] = {2,kViolet+7,4};
  for (int iYear=0; iYear< nYears ; ++iYear) {
    gr[iYear] = new TGraphAsymmErrors(nPars,x[iYear],y[iYear],xe,xe,yl[iYear],yh[iYear]);
    gr[iYear]->SetName(Form("gr%i",iYear));
    gr[iYear]->SetLineColor(colours[iYear]);
    gr[iYear]->SetLineWidth(2);
    gr[iYear]->Draw("P");
  }

  TLegend leg (0.15,0.66,0.4,0.9);
  leg.AddEntry(grSim,"Simultaneous result","f");
  if (q2Bin==4 || q2Bin==6) leg.AddEntry(grPDG,"PDG value","lep");
  for (int iYear=0; iYear< nYears; ++iYear)
    leg.AddEntry(gr[iYear], Form("Result %i", haveSingleFits[iYear]),"lep");
  leg.Draw();

  canv.SaveAs(("plotSimFit4d_d/comparisonSimFit_result_"+shortString+"_XGBv8.pdf").c_str());

  // Plot S-wave parameters
  double xS[nYears][nSPars];
  double xSe[nSPars];
  double xSSim[nSPars];
  double xSeSim[nSPars];

  double ySSim [nSPars];
  double ySlSim[nSPars];
  double yShSim[nSPars];
  double yS [nYears][nSPars];
  double ySl[nYears][nSPars];
  double ySh[nYears][nSPars];

  for (int iPar=0; iPar<nSPars; ++iPar) {
    xSSim[iPar]  = 0.5+iPar;
    xSeSim[iPar] = 0.5;
    xSe[iPar]    = 1.0/(2*nYears);
    ySSim [iPar] = 0;
    ySlSim[iPar] = -1 * lowSParSim [iPar];
    yShSim[iPar] = highSParSim[iPar];
    for (int iYear=0; iYear<nYears; ++iYear) {
      xS [iYear][iPar] = (2.0*iYear+1)/(2*nYears)+iPar;
      yS [iYear][iPar] = bestSPar[iYear][iPar] - bestSParSim[iPar];
      ySl[iYear][iPar] = -1 * lowSPar [iYear][iPar];
      ySh[iYear][iPar] = highSPar[iYear][iPar];
    }
  }

  TCanvas canvS("canvS","canvS",1500,1000);
  canvS.cd();

  string plotTitleS = Form("Comparison of simultaneous fit S-wave parameters - q2 bin %i;;(single-year result) - (simultaneous result)",q2Bin);
  auto hLabS = new TH1S ("hLabS",plotTitleS.c_str(),nSPars,0,nSPars);
  for (int iBin=0; iBin<nSPars; ++iBin)
    hLabS->GetXaxis()->SetBinLabel(iBin+1,sParName[iBin].c_str());
  double yRangeS = 0.2;
  if (q2Bin==4) yRangeS = 0.1;
  hLabS->SetMinimum(-1*yRangeS);
  hLabS->SetMaximum(yRangeS);
  hLabS->Draw();

  auto grSimS = new TGraphAsymmErrors(nSPars,xSSim,ySSim,xSeSim,xSeSim,ySlSim,yShSim);
  grSimS->SetName("grSimS");
  grSimS->SetFillColor(30);
  grSimS->Draw("2");

  TLine lineS (0,0,nSPars,0);
  lineS.SetLineWidth(2);
  lineS.SetLineStyle(2);
  lineS.SetLineColor(13);
  lineS.Draw();

  TGraphAsymmErrors* grS [nYears];
  for (int iYear=0; iYear<nYears; ++iYear) {
    grS[iYear] = new TGraphAsymmErrors(nSPars,xS[iYear],yS[iYear],xSe,xSe,ySl[iYear],ySh[iYear]);
    grS[iYear]->SetName(Form("grS%i",iYear));
    grS[iYear]->SetLineColor(colours[iYear]);
    grS[iYear]->SetLineWidth(2);
    grS[iYear]->Draw("P");
  }

  TLegend legS (0.65,0.73,0.9,0.9);
  legS.AddEntry(grSimS,"Simultaneous result","f");
  for (int iYear=0; iYear<nYears; ++iYear)
    legS.AddEntry(grS[iYear], Form("Result %i", haveSingleFits[iYear]),"lep");
  legS.Draw();

  canvS.SaveAs(("plotSimFit4d_d/comparisonSimFit_Swave_"+shortString+"_XGBv8.pdf").c_str());


  // Plot mass parameters
  double xMass[nYears*(nMassPars+nMassConstPars)];
  double xMasse[nYears*(nMassPars+nMassConstPars)];

  double yMassSim [nYears*(nMassPars+nMassConstPars)];
  double yMasseSim[nYears*(nMassPars+nMassConstPars)];
  double yMass [nYears*(nMassPars+nMassConstPars)];
  double yMasse[nYears*(nMassPars+nMassConstPars)];
  double yMassConst [nYears*(nMassPars+nMassConstPars)];
  double yMasseConst[nYears*(nMassPars+nMassConstPars)];

  for (int iPar=0; iPar<nYears*(nMassPars+nMassConstPars); ++iPar) {
    xMass[iPar]  = 0.5+iPar;
    xMasse[iPar] = 0.5;
    yMassSim [iPar] = 0.;
    yMasseSim[iPar] = 1.;
  }

  for (int iPar=0; iPar<nMassPars; ++iPar) {
    for (int iYear = 0; iYear < nYears; iYear++){
      int this_index = iPar*nYears + iYear;
      int year_index = haveSingleFits[iYear] - 2016;
      yMass [this_index] = (bestMassPar[iYear][iPar] - bestMassParSim[year_index][iPar]) / errMassParSim[year_index][iPar];
      yMasse[this_index] = errMassPar[iYear][iPar] / errMassParSim[year_index][iPar];
      yMassConst [this_index] = 1e99;
      yMasseConst[this_index] = 0;
    }
  }  


  for (int iPar=0; iPar<nMassConstPars; ++iPar) {
    for (int iYear = 0; iYear < nYears; iYear++){
      int this_index = nMassPars*nYears + iPar*nYears + iYear;
      int year_index = haveSingleFits[iYear] - 2016;
      yMass [this_index] = (bestMassConstPar[iYear][iPar] - bestMassConstParSim[year_index][iPar]) / errMassConstParSim[year_index][iPar];
      yMasse[this_index] =  errMassConstPar[iYear][iPar] / errMassConstParSim[year_index][iPar];
      yMassConst [this_index] = (massConstVal[year_index][iPar] - bestMassConstParSim[year_index][iPar]) / errMassConstParSim[year_index][iPar];
      yMasseConst[this_index] = massConstErr[year_index][iPar] / errMassConstParSim[year_index][iPar];
    }
  }

  TCanvas canvM("canvM","canvM",1500,1000);
  canvM.cd();

  string plotTitleM = Form("Comparison of simultaneous fit mass parameters - q2 bin %i;;pulls",q2Bin);
  auto hLabM = new TH1S ("hLabM",plotTitleM.c_str(),nYears*(nMassPars+nMassConstPars),0,nYears*(nMassPars+nMassConstPars));
  for (int iBin=0; iBin<nYears*nMassPars; ++iBin){
    string this_year = Form("%i", haveSingleFits[iBin%nYears]);
    hLabM->GetXaxis()->SetBinLabel(iBin+1,Form(massParName[iBin/nYears].c_str(), this_year.c_str()));
  }
  for (int iBin=0; iBin<nYears; ++iBin){
    string this_year = Form("%i", haveSingleFits[iBin%nYears]);
    hLabM->GetXaxis()->SetBinLabel(iBin+1+nYears*nMassPars,Form("R^{%s}", this_year.c_str()));
  }  
  for (int iBin=nYears; iBin<nYears*nMassConstPars; ++iBin){
    string this_year = Form("%i", haveSingleFits[iBin%nYears]);
    hLabM->GetXaxis()->SetBinLabel(iBin+1+nYears*nMassPars,Form(massConstParName[iBin/nYears].c_str(), this_year.c_str()));
    if (massConstParName[iBin/nYears].find("eltaPeak") !=std::string::npos)
      hLabM->GetXaxis()->SetBinLabel(iBin+1+nYears*nMassPars,Form("#Deltam^{%s}", this_year.c_str()));      
  
  }
  double yRangeM = 9;
  if (q2Bin==4) yRangeM = 30;//30;
  hLabM->SetMinimum(-1*yRangeM);
  hLabM->SetMaximum(yRangeM);
  hLabM->Draw();
  hLabM->GetXaxis()->LabelsOption("v");


  auto grSimM = new TGraphErrors(nYears*(nMassPars+nMassConstPars),xMass,yMassSim,xMasse,yMasseSim);
  grSimM->SetName("grSimM");
  grSimM->SetFillColor(30);
  grSimM->Draw("2");

  TLine lineM (0,0,nYears*(nMassPars+nMassConstPars),0);
  lineM.SetLineWidth(2);
  lineM.SetLineStyle(2);
  lineM.SetLineColor(13);
  lineM.Draw();

  auto grM = new TGraphErrors(nYears*(nMassPars+nMassConstPars),xMass,yMass,xMasse,yMasse);
  grM->SetName("grM");
  grM->SetLineColor(2);
  grM->SetLineWidth(2);
  grM->Draw("P");

  auto grConst = new TGraphErrors(nYears*(nMassPars+nMassConstPars),xMass,yMassConst,xMasse,yMasseConst);
  grConst->SetName("grConst");
  grConst->SetLineColor(1);\
  grConst->SetLineWidth(2);
  grConst->SetFillStyle(0);
  grConst->SetMarkerStyle(20);
  grConst->Draw("5p");

  TLegend legM (0.1,0.73,0.35,0.9);
  legM.AddEntry(grSimM,"Simultaneous result","f");
  legM.AddEntry(grM,"Single-year result","lep");
  legM.AddEntry(grConst,"Constraint","fp");
  legM.Draw();

  canvM.SaveAs(("plotSimFit4d_d/comparisonSimFit_massParameters_"+shortString+"_XGBv8.pdf").c_str());
  
  // Plot contrained mass parameters with traditional pull plots
  {
    double xMass[3*nMassConstPars];
    double xMasse[3*nMassConstPars];

    double yMassSim [3*nMassConstPars];
    double yMasseSim[3*nMassConstPars];
    double yMass [3*nMassConstPars];
    double yMasse[3*nMassConstPars];
    double yMassConst [3*nMassConstPars];
    double yMasseConst[3*nMassConstPars];

    for (int iPar=0; iPar<3*nMassConstPars; ++iPar) {
      xMass[iPar]  = 0.5+iPar;
      xMasse[iPar] = 0.5;
      yMassConst [iPar] = 0.;
      yMasseConst[iPar] = 1.;
    }
    for (int iPar=0; iPar<nMassConstPars; ++iPar) {
      for (int iYear=0; iYear<3; ++iYear) {
        int this_index = iPar*3 + iYear;
        yMassSim [this_index] = (bestMassConstParSim[iYear][iPar] - massConstVal[iYear][iPar]) / massConstErr[iYear][iPar];
        yMasseSim[this_index] = errMassConstParSim[iYear][iPar] / massConstErr[iYear][iPar];
      }
    }

    TCanvas canvM("canvPull","canvPull",1000,1500);
    auto vp = canvM.cd();
    vp->SetGridx();

    string plotTitleM = Form("Pulls of constrained fit parameters - q2 bin %i;(#theta_{fit} - #theta_{constr})/#sigma(#theta)_{constr};",q2Bin);
    double yRangeM = 9;
    if (q2Bin==4) yRangeM = 20;
    auto hLabM = new TH2S ("hLabPull",plotTitleM.c_str(),
			   1,-1*yRangeM,yRangeM,
			   3*nMassConstPars,0,3*nMassConstPars);
    for (int iBin=0; iBin<3; ++iBin)
      hLabM->GetYaxis()->SetBinLabel(iBin+1,Form("R^{%s}",year[iBin%3].c_str()));
    for (int iBin=3; iBin<3*nMassConstPars; ++iBin){
      hLabM->GetYaxis()->SetBinLabel(iBin+1,Form(massConstParName[iBin/3].c_str(),year[iBin%3].c_str()));
      if (massConstParName[iBin/3].find("eltaPeak") !=std::string::npos)
        hLabM->GetYaxis()->SetBinLabel(iBin+1,Form("#Deltam^{%s}",year[iBin%3].c_str()));
    }    
    hLabM->Draw();

    auto grSimM = new TGraphErrors(3*nMassConstPars,yMassSim,xMass,yMasseSim,xMasse);
    grSimM->SetName("grSimPull");
    grSimM->SetLineColor(30);
    grSimM->SetLineWidth(2);
    grSimM->Draw("P");

    canvM.SaveAs(("plotSimFit4d_d/comparisonSimFit_massParametersPulls_"+shortString+"_XGBv8.pdf").c_str());

  }

}
