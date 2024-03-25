#include "interface/utils.h"

string year[3] = {"2016","2017","2018"};

void plotSimFitComparison(int q2Bin = 6, int q2Stat = -1)
{
  string shortString = Form("b%ip1",q2Bin);
  string statString = q2Stat<0 ? "" : Form("_b%istat-0",q2Stat);

  static const int nPars = 8;
  string parName[nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

  static const int nMassPars = 3;
  string massParName[nMassPars] = {"mean_{RT}^{%s}",
				   "fsig_"+shortString+"_%s",
				   "slope^{%s}"};

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
  
  double bestMassParSim[3*nMassPars];
  double errMassParSim [3*nMassPars];
  double bestMassConstParSim[3*nMassConstPars];
  double errMassConstParSim [3*nMassConstPars];
  double bestMassPar[3*nMassPars];
  double errMassPar [3*nMassPars];
  double bestMassConstPar[3*nMassConstPars];
  double errMassConstPar [3*nMassConstPars];

  double massConstVal[3*nMassConstPars];
  double massConstErr[3*nMassConstPars];

  // to be cross checked if exists (and if the parity is correct)
  string finName = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/simFitResults/simFitResult_data_fullAngularMass_Swave_%s_b%ip1_XGBv8.root";
  //// use the following for plotting results from sub-samples of the data control channel
  // string finName = "simFitResults4d/simFitResult_data_fullAngularMass_Swave_%s" + statString + "_b%i-XGBv8_wp90_unbl4.root";
  string finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-LMNR-XGBv8/%s.root";
  if (q2Bin==4) finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-Jpsi-XGBv8/%s.root";
  else if (q2Bin==6) finName2 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-Psi-XGBv8/%s.root";

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

  auto wsSim = (RooWorkspace*)finSim->Get("wsp_out");
  if ( !wsSim || wsSim->IsZombie() ) {
    cout<<"Sim workspace is problematic!"<<endl;
    return;
  }
  for (int iPar=0; iPar<nSPars; ++iPar) {
    auto par = (RooRealVar*)wsSim->var(sParName[iPar].c_str());
    bestSParSim[iPar] = par->getValV();
    lowSParSim [iPar] = par->getErrorLo();
    highSParSim[iPar] = par->getErrorHi();
  }
  for (int iPar=0; iPar<3*nMassPars; ++iPar) {
    auto par = (RooRealVar*)wsSim->var(Form(massParName[iPar/3].c_str(),year[iPar%3].c_str()));
    bestMassParSim[iPar] = par->getValV();
    errMassParSim [iPar] = par->getError();
  }
  for (int iPar=0; iPar<3*nMassConstPars; ++iPar) {
    auto par = (RooRealVar*)wsSim->var(Form(massConstParName[iPar/3].c_str(),year[iPar%3].c_str()));
    bestMassConstParSim[iPar] = par->getValV();
    errMassConstParSim [iPar] = par->getError();
  }
  
  finSim->Close();


  for (int iy=0; iy<3; ++iy) {
    auto fin = TFile::Open(Form(finName.c_str(),year[iy].c_str(),q2Bin));
    if ( !fin || fin->IsZombie() ) {
      cout<<year[iy]<<" file is problematic!"<<endl;
      return;
    }
    auto fin2 = TFile::Open(Form(finName2.c_str(),year[iy].c_str()));
    if ( !fin2 || fin2->IsZombie() ) {
      cout<<year[iy]<<" file2 is problematic!"<<endl;
      return;
    }

    auto tin= (TTree*)fin->Get("fitResultsTree");
    if ( !tin || tin->IsZombie() || tin->GetEntries() < 1 ) {
      cout<<year[iy]<<" tree is problematic!"<<endl;
      fin->Close();
      return;
    }
    for (int iPar=0; iPar<nPars; ++iPar) {
      tin->SetBranchAddress((parName[iPar]+"_best").c_str(),&bestPar[iy][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_low" ).c_str(),&lowPar [iy][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_high").c_str(),&highPar[iy][iPar]);
    }
    tin->GetEntry(0);

    auto ws = (RooWorkspace*)fin->Get("wsp_out");
    if ( !ws || ws->IsZombie() ) {
      cout<<year[iy]<<" workspace is problematic!"<<endl;
      return;
    }
    auto ws2 = (RooWorkspace*)fin2->Get("w");
    if ( !ws2 || ws2->IsZombie() ) {
      cout<<year[iy]<<" workspace2 is problematic!"<<endl;
      return;
    }

    for (int iPar=0; iPar<nSPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(sParName[iPar].c_str());
      bestSPar[iy][iPar] = par->getValV();
      lowSPar [iy][iPar] = par->getErrorLo();
      highSPar[iy][iPar] = par->getErrorHi();
    }
    for (int iPar=0; iPar<nMassPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massParName[iPar].c_str(),year[iy].c_str()));
      bestMassPar[iPar*3+iy] = par->getValV();
      errMassPar [iPar*3+iy] = par->getError();
    }
    ws2->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    double mean_rt_val = ws2->var(Form("mean_{RT}^{%i}", q2Bin))->getVal(); 
    double mean_rt_err = ws2->var(Form("mean_{RT}^{%i}", q2Bin))->getError(); 
    for (int iPar=0; iPar<nMassConstPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massConstParName[iPar].c_str(),year[iy].c_str()));
      bestMassConstPar[iPar*3+iy] = par->getValV();
      errMassConstPar [iPar*3+iy] = par->getError();
      if (iPar==8) ws2->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
      if (iPar==11) {
	double mean_wt_val = ws2->var(Form("mean_{WT}^{%i}", q2Bin))->getVal(); 
	double mean_wt_err = ws2->var(Form("mean_{WT}^{%i}", q2Bin))->getError(); 
	massConstVal[iPar*3+iy] = mean_rt_val-mean_wt_val;
	massConstErr[iPar*3+iy] = sqrt(mean_rt_err*mean_rt_err+mean_wt_err*mean_wt_err); 
      } else if (iPar>0) {
	auto par2 = (RooRealVar*)ws2->var(Form(massConstParName[iPar].c_str(),Form("%i",q2Bin)));
	if ( !par2 || par2->IsZombie() ) {
	  cout<<Form(massConstParName[iPar].c_str(),year[iy].c_str())<<" constraint not found"<<endl;
	  return;
	}
	massConstVal[iPar*3+iy] = par2->getValV();
	massConstErr[iPar*3+iy] = par2->getError();
	if (q2Stat>=0) massConstErr[iPar*3+iy] *= constr_width_fac[atoi(year[iy].c_str())][2*q2Stat+iPar/8];
      } else {
	double nrt_mc = ws2->var(Form("nRT_%i",q2Bin))->getVal();
	double nwt_mc = ws2->var(Form("nWT_%i",q2Bin))->getVal();
	massConstVal[iy] = 1.;
	massConstErr[iy] = fM_sigmas[atoi(year[iy].c_str())][q2Stat<0?q2Bin:q2Stat] / ( nwt_mc / (nrt_mc + nwt_mc) );
	cout<<"==========="
	    <<iy<<"\t"<<nrt_mc<<"\t"<<nwt_mc<<"\t"
	    <<massConstErr[iy]<<"\t"
	    <<bestMassConstPar[iy]<<"\t"
	    <<errMassConstPar [iPar*3+iy]<<endl;
      }
    }

    fin->Close();
    fin2->Close();

  }

  // Plot P-wave parameters
  double x[3][nPars];
  double xe[nPars];
  double xSim[nPars];
  double xeSim[nPars];

  double ySim [nPars];
  double ylSim[nPars];
  double yhSim[nPars];
  double y [3][nPars];
  double yl[3][nPars];
  double yh[3][nPars];

  for (int iPar=0; iPar<nPars; ++iPar) {
    xSim[iPar]  = 0.5+iPar;
    xeSim[iPar] = 0.5;
    xe[iPar]    = 1.0/6;
    ySim [iPar] = 0;
    ylSim[iPar] = -1 * lowParSim [iPar];
    yhSim[iPar] = highParSim[iPar];
    for (int iy=0; iy<3; ++iy) {
      x [iy][iPar] = (2.0*iy+1)/6+iPar;
      y [iy][iPar] = bestPar[iy][iPar] - bestParSim[iPar];
      yl[iy][iPar] = -1 * lowPar [iy][iPar];
      yh[iy][iPar] = highPar[iy][iPar];
    }
  }

  TCanvas canv("canv","canv",1500,1000);
  canv.cd();

  gStyle->SetOptStat(0);
  // cCover.cd()->SetTicky(2);
  string plotTitle = Form("Fit results - q2 bin %i%s;;(single-year result) - (simultaneous result)",q2Bin,q2Stat>=0?Form(" subsample (q2 bin %i stat)",q2Stat):"");
  auto hLab = new TH1S ("hLab",plotTitle.c_str(),nPars,0,nPars);
  for (int iBin=0; iBin<nPars; ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,parName[iBin].c_str());
  double yRange = 0.1;
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

  TGraphAsymmErrors* gr [3];
  int colours [3] = {2,kViolet+7,4};
  for (int iy=0; iy<3; ++iy) {
    gr[iy] = new TGraphAsymmErrors(nPars,x[iy],y[iy],xe,xe,yl[iy],yh[iy]);
    gr[iy]->SetName(Form("gr%i",iy));
    gr[iy]->SetLineColor(colours[iy]);
    gr[iy]->SetLineWidth(2);
    gr[iy]->Draw("P");
  }

  TLegend leg (0.15,0.66,0.4,0.9);
  leg.AddEntry(grSim,"Simultaneous result","f");
  if (q2Bin==4 || q2Bin==6) leg.AddEntry(grPDG,"PDG value","lep");
  for (int iy=0; iy<3; ++iy)
    leg.AddEntry(gr[iy],("Result "+year[iy]).c_str(),"lep");
  leg.Draw();

  canv.SaveAs(("plotSimFit4d_d/comparisonSimFit_result_"+shortString+statString+"_XGBv8.pdf").c_str());

  // Plot S-wave parameters
  double xS[3][nSPars];
  double xSe[nSPars];
  double xSSim[nSPars];
  double xSeSim[nSPars];

  double ySSim [nSPars];
  double ySlSim[nSPars];
  double yShSim[nSPars];
  double yS [3][nSPars];
  double ySl[3][nSPars];
  double ySh[3][nSPars];

  for (int iPar=0; iPar<nSPars; ++iPar) {
    xSSim[iPar]  = 0.5+iPar;
    xSeSim[iPar] = 0.5;
    xSe[iPar]    = 1.0/6;
    ySSim [iPar] = 0;
    ySlSim[iPar] = -1 * lowSParSim [iPar];
    yShSim[iPar] = highSParSim[iPar];
    for (int iy=0; iy<3; ++iy) {
      xS [iy][iPar] = (2.0*iy+1)/6+iPar;
      yS [iy][iPar] = bestSPar[iy][iPar] - bestSParSim[iPar];
      ySl[iy][iPar] = -1 * lowSPar [iy][iPar];
      ySh[iy][iPar] = highSPar[iy][iPar];
    }
  }

  TCanvas canvS("canvS","canvS",1500,1000);
  canvS.cd();

  string plotTitleS = Form("Fit S-wave parameters - q2 bin %i%s;;(single-year result) - (simultaneous result)",q2Bin,q2Stat>=0?Form(" subsample (q2 bin %i stat)",q2Stat):"");
  auto hLabS = new TH1S ("hLabS",plotTitleS.c_str(),nSPars,0,nSPars);
  for (int iBin=0; iBin<nSPars; ++iBin)
    hLabS->GetXaxis()->SetBinLabel(iBin+1,sParName[iBin].c_str());
  double yRangeS = 0.2;
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

  TGraphAsymmErrors* grS [3];
  for (int iy=0; iy<3; ++iy) {
    grS[iy] = new TGraphAsymmErrors(nSPars,xS[iy],yS[iy],xSe,xSe,ySl[iy],ySh[iy]);
    grS[iy]->SetName(Form("grS%i",iy));
    grS[iy]->SetLineColor(colours[iy]);
    grS[iy]->SetLineWidth(2);
    grS[iy]->Draw("P");
  }

  TLegend legS (0.65,0.73,0.9,0.9);
  legS.AddEntry(grSimS,"Simultaneous result","f");
  for (int iy=0; iy<3; ++iy)
    legS.AddEntry(grS[iy],("Result "+year[iy]).c_str(),"lep");
  legS.Draw();

  canvS.SaveAs(("plotSimFit4d_d/comparisonSimFit_Swave_"+shortString+statString+"_XGBv8.pdf").c_str());

  // Plot mass parameters
  double xMass[3*(nMassPars+nMassConstPars)];
  double xMasse[3*(nMassPars+nMassConstPars)];

  double yMassSim [3*(nMassPars+nMassConstPars)];
  double yMasseSim[3*(nMassPars+nMassConstPars)];
  double yMass [3*(nMassPars+nMassConstPars)];
  double yMasse[3*(nMassPars+nMassConstPars)];
  double yMassConst [3*(nMassPars+nMassConstPars)];
  double yMasseConst[3*(nMassPars+nMassConstPars)];

  for (int iPar=0; iPar<3*(nMassPars+nMassConstPars); ++iPar) {
    xMass[iPar]  = 0.5+iPar;
    xMasse[iPar] = 0.5;
    yMassSim [iPar] = 0.;
    yMasseSim[iPar] = 1.;
  }
  for (int iPar=0; iPar<3*nMassPars; ++iPar) {
    yMass [iPar] = (bestMassPar[iPar] - bestMassParSim[iPar]) / errMassParSim[iPar];
    yMasse[iPar] = errMassPar[iPar] / errMassParSim[iPar];
    yMassConst [iPar] = 1e99;
    yMasseConst[iPar] = 0;
  }
  for (int iPar=0; iPar<3*nMassConstPars; ++iPar) {
    yMass [iPar+3*nMassPars] = (bestMassConstPar[iPar] - bestMassConstParSim[iPar]) / errMassConstParSim[iPar];
    yMasse[iPar+3*nMassPars] = errMassConstPar[iPar] / errMassConstParSim[iPar];
    yMassConst [iPar+3*nMassPars] = (massConstVal[iPar] - bestMassConstParSim[iPar]) / errMassConstParSim[iPar];
    yMasseConst[iPar+3*nMassPars] = massConstErr[iPar] / errMassConstParSim[iPar];
  }

  TCanvas canvM("canvM","canvM",1500,1000);
  canvM.cd();

  string plotTitleM = Form("Fit mass parameters - q2 bin %i%s;;pulls",q2Bin,q2Stat>=0?Form(" subsample (q2 bin %i stat)",q2Stat):"");
  auto hLabM = new TH1S ("hLabM",plotTitleM.c_str(),3*(nMassPars+nMassConstPars),0,3*(nMassPars+nMassConstPars));
  for (int iBin=0; iBin<3*nMassPars; ++iBin)
    hLabM->GetXaxis()->SetBinLabel(iBin+1,Form(massParName[iBin/3].c_str(),year[iBin%3].c_str()));
  for (int iBin=0; iBin<3; ++iBin)
    hLabM->GetXaxis()->SetBinLabel(iBin+1+3*nMassPars,Form("R^{%s}",year[iBin%3].c_str()));
  for (int iBin=3; iBin<3*nMassConstPars; ++iBin)
    hLabM->GetXaxis()->SetBinLabel(iBin+1+3*nMassPars,Form(massConstParName[iBin/3].c_str(),year[iBin%3].c_str()));
  double yRangeM = q2Stat>=0 ? 3 : 10.9;
  hLabM->SetMinimum(-1*yRangeM);
  hLabM->SetMaximum(yRangeM);
  hLabM->Draw();

  auto grSimM = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMassSim,xMasse,yMasseSim);
  grSimM->SetName("grSimM");
  grSimM->SetFillColor(30);
  grSimM->Draw("2");

  TLine lineM (0,0,3*(nMassPars+nMassConstPars),0);
  lineM.SetLineWidth(2);
  lineM.SetLineStyle(2);
  lineM.SetLineColor(13);
  lineM.Draw();

  auto grM = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMass,xMasse,yMasse);
  grM->SetName("grM");
  grM->SetLineColor(2);
  grM->SetLineWidth(2);
  grM->Draw("P");

  auto grConst = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMassConst,xMasse,yMasseConst);
  grConst->SetName("grConst");
  grConst->SetLineColor(1);
  grConst->SetLineWidth(2);
  grConst->SetFillStyle(0);
  grConst->SetMarkerStyle(20);
  grConst->Draw("5p");

  TLegend legM (0.1,0.73,0.35,0.9);
  legM.AddEntry(grSimM,"Simultaneous result","f");
  legM.AddEntry(grM,"Single-year result","lep");
  legM.AddEntry(grConst,"Constraint","fp");
  legM.Draw();

  canvM.SaveAs(("plotSimFit4d_d/comparisonSimFit_massParameters_"+shortString+statString+"_XGBv8.pdf").c_str());

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
    for (int iPar=0; iPar<3*nMassConstPars; ++iPar) {
      yMass [iPar] = (bestMassConstPar[iPar] - massConstVal[iPar]) / massConstErr[iPar];
      yMasse[iPar] = errMassConstPar[iPar] / massConstErr[iPar];
      yMassSim [iPar] = (bestMassConstParSim[iPar] - massConstVal[iPar]) / massConstErr[iPar];
      yMasseSim[iPar] = errMassConstParSim[iPar] / massConstErr[iPar];
    }

    TCanvas canvM("canvPull","canvPull",1000,1500);
    auto vp = canvM.cd();
    vp->SetGridx();

    string plotTitleM = Form("Pulls of contrained fit parameters - q2 bin %i%s;(#theta_{fit} - #theta_{constr})/#sigma(#theta)_{constr};",q2Bin,q2Stat>=0?Form(" subsample (q2 bin %i stat)",q2Stat):"");
    auto hLabM = new TH2S ("hLabPull",plotTitleM.c_str(),
			   1,-1*yRangeM,yRangeM,
			   3*nMassConstPars,0,3*nMassConstPars);
    for (int iBin=0; iBin<3; ++iBin)
      hLabM->GetYaxis()->SetBinLabel(iBin+1,Form("R^{%s}",year[iBin%3].c_str()));
    for (int iBin=3; iBin<3*nMassConstPars; ++iBin)
      hLabM->GetYaxis()->SetBinLabel(iBin+1,Form(massConstParName[iBin/3].c_str(),year[iBin%3].c_str()));
    hLabM->Draw();

    auto grSimM = new TGraphErrors(3*nMassConstPars,yMassSim,xMass,yMasseSim,xMasse);
    grSimM->SetName("grSimPull");
    grSimM->SetLineColor(30);
    grSimM->SetLineWidth(2);
    grSimM->Draw("P");

    // TLine lineM (0,0,0,3*nMassConstPars);
    // lineM.SetLineWidth(2);
    // lineM.SetLineStyle(2);
    // lineM.SetLineColor(13);
    // lineM.Draw();

    // auto grM = new TGraphErrors(3*nMassConstPars,yMass,xMass,yMasse,xMasse);
    // grM->SetName("grPull");
    // grM->SetLineColor(2);
    // grM->SetLineWidth(2);
    // grM->Draw("P");

    // auto grConst = new TGraphErrors(3*nMassConstPars,xMass,yMassConst,xMasse,yMasseConst);
    // grConst->SetName("grConst");
    // grConst->SetLineColor(1);
    // grConst->SetLineWidth(2);
    // grConst->SetFillStyle(0);
    // grConst->SetMarkerStyle(20);
    // grConst->Draw("5p");

    // TLegend legM (0.1,0.73,0.35,0.9);
    // legM.AddEntry(grSimM,"Simultaneous result","lep");
    // legM.AddEntry(grM,"Single-year result","lep");
    // // legM.AddEntry(grConst,"Constraint","fp");
    // legM.Draw();

    canvM.SaveAs(("plotSimFit4d_d/comparisonSimFit_massParametersPulls_"+shortString+statString+"_XGBv8.pdf").c_str());

  }

}
