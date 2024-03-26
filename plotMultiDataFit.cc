#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

static const int nBins = 8;
float binBorders [nBins+1] = { 1.1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};

static const int nUncBins = 100;
double binsUnc [nUncBins+1];

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};
string parTitle[nPars] = {"F_{L}","P_{1}","P_{2}","P_{3}","P'_{4}","P'_{5}","P'_{6}","P'_{8}"};
double parMin  [nPars] = {0,-1,-0.5,-0.5,-1*sqrt(2),-1*sqrt(2),-1*sqrt(2),-1*sqrt(2)};
double parMax  [nPars] = {1, 1, 0.5, 0.5,   sqrt(2),   sqrt(2),   sqrt(2),   sqrt(2)};

static const int nQuant = 4;
double quantPerc [nQuant] = {0.025,0.16,0.84,0.975};

int colors [12] = { 633, 417, 879, 857, 839, 801, 921, 607, 807, 419, 907, 402 };
// int colors [13] = { 633, 417, 879, 857, 839, 887, 801, 921, 607, 807, 419, 907, 402 };

static const int nStats = 2;
int aq2stat [nStats] = { 0, 5 };

double diffMax = 0.0499;

void plotMultiDataFit ()
{

  gROOT->SetBatch(true);

  vector< vector<TH1D*> > vHistBest (nPars);
  vector< vector<TH1D*> > vHistPull (nPars);
  vector< vector<TH1D*> > vHistErrH (nPars);
  vector< vector<TH1D*> > vHistErrL (nPars);

  vector< vector<double> > vHistBestRECO (nPars);
  vector< vector<double> > vHistErrHRECO (nPars);
  vector< vector<double> > vHistErrLRECO (nPars);

  vector< vector<double> > vMean (nPars);
  vector< vector<double> > vRMS  (nPars);
  vector< vector<double> > vBias (nPars);
  vector< vector<double> > vMeanErr (nPars);

  vector<int> vq2Bins (0);

  binsUnc[0]=0.006;
  for (int i=0; i<nUncBins; ++i)
    binsUnc[i+1] = binsUnc[i] * pow(0.4/binsUnc[0],1./nUncBins); // to reach 0.4 as largest error

  double q2Val [nBins];
  double q2Err [nBins];

  for (int i=0; i<nBins; ++i) {
    q2Val[i] = 0.5 * (binBorders[i+1]+binBorders[i]);
    q2Err[i] = 0.5 * (binBorders[i+1]-binBorders[i]);
  }

  for (int iColor = 0; iColor<nStats; ++iColor) {

    int q2stat = aq2stat[iColor];

    int q2Bin = 4;
    vq2Bins.push_back(q2stat);

    TChain fitResultsTree ("fitResultsTree","");
    string filename = Form("/eos/cms/store/user/fiorendi/p5prime/angularFits/dataCR/signal-stat/simFitResult_data_fullAngularMass_Swave_201620172018_b%istat-*_b%i-XGBv8_wp90_unbl4.root",q2stat,q2Bin);
    fitResultsTree.Add(filename.c_str());

//     string filename_fR = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/simFitResults/simFitResult_data_fullAngularMass_Swave_201620172018_b%ip1_XGBv8.root",q2Bin); // link will be updated
    string filename_fR = Form("/afs/cern.ch/work/d/dini/public/ResonantSwapCut/Bin4SimBern/simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_unbl4.root",q2Bin);
    TFile* filein_fR = TFile::Open(filename_fR.c_str());
    TTree* fitResultsTree_fR = (TTree*)filein_fR->Get("fitResultsTree");
    if (!fitResultsTree_fR || fitResultsTree_fR->GetEntries() != 1) {
      cout<<"Error, unexpected number of entries in fitResultsTree in file: "<<filename_fR<<endl;
      return;
    }

    int nSamp = fitResultsTree.GetEntries();
    cout<<"Number of samples: "<<nSamp<<endl;

    vector<double> vBest(nPars);
    vector<double> vHigh(nPars);
    vector<double> vLow (nPars);

    for (int iPar=0; iPar<nPars; ++iPar) {
      
      fitResultsTree.SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vBest[iPar]);
      fitResultsTree.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh[iPar]);
      fitResultsTree.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow [iPar]);
      
      fitResultsTree_fR->SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vBest[iPar]);
      fitResultsTree_fR->SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh[iPar]);
      fitResultsTree_fR->SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow [iPar]);

      vHistBest[iPar].push_back( new TH1D(Form("hBest%i%i",iColor,iPar),Form("%s results of fits to control-region sub-samples with signal q2-bin statistics;%s;# of results",parTitle[iPar].c_str(),parTitle[iPar].c_str()),100,parMin[iPar],parMax[iPar]) );
      vHistPull[iPar].push_back( new TH1D(Form("hPull%i%i",iColor,iPar),Form("%s pulls in fits to control-region sub-samples with signal q2-bin statistics;%s pulls;# of results",parTitle[iPar].c_str(),parTitle[iPar].c_str()),30,-5,5) );
      vHistErrH[iPar].push_back( new TH1D(Form("hErrH%i%i",iColor,iPar),Form("%s MINOS uncertainties of fits to control-region sub-samples with signal q2-bin statistics;#sigma(%s);# of results",parTitle[iPar].c_str(),parTitle[iPar].c_str()),nUncBins,binsUnc) );
      vHistErrL[iPar].push_back( new TH1D(Form("hErrL%i%i",iColor,iPar),Form("%s MINOS uncertainties of fits to control-region sub-samples with signal q2-bin statistics;#sigma(%s);# of results",parTitle[iPar].c_str(),parTitle[iPar].c_str()),nUncBins,binsUnc) );

      vHistBest[iPar].back()->SetLineColor(colors[iColor]);
      vHistPull[iPar].back()->SetLineColor(colors[iColor]);
      vHistErrH[iPar].back()->SetLineColor(colors[iColor]);
      vHistErrL[iPar].back()->SetLineColor(colors[iColor]);
      vHistErrL[iPar].back()->SetFillColor(colors[iColor]);
      vHistErrL[iPar].back()->SetFillStyle(3345);

      vMean[iPar].push_back( 0 );
      vRMS [iPar].push_back( 0 );
      vBias[iPar].push_back( 0 );
      vMeanErr[iPar].push_back( 0 );

    }

    fitResultsTree_fR->GetEntry(0);
    for (int iPar=0; iPar<nPars; ++iPar) {
      vHistBestRECO[iPar].push_back( vBest[iPar] );
      vHistErrHRECO[iPar].push_back( vHigh[iPar] );
      vHistErrLRECO[iPar].push_back( vLow [iPar] );
    }

    for (int iEn=0; iEn<fitResultsTree.GetEntries(); ++iEn) {
      fitResultsTree.GetEntry(iEn);
      for (int iPar=0; iPar<nPars; ++iPar) {
	vHistBest[iPar].back()->Fill(vBest[iPar]);
	vHistPull[iPar].back()->Fill((vBest[iPar]-vHistBestRECO[iPar].back())/(vBest[iPar]>vHistBestRECO[iPar].back()?abs(vLow[iPar]):abs(vHigh[iPar])));
	vHistErrH[iPar].back()->Fill(vHigh[iPar]);
	vHistErrH[iPar].back()->Fill(abs(vLow [iPar])); // To create a stacked histogram
	vHistErrL[iPar].back()->Fill(abs(vLow [iPar]));

	vMean[iPar].back() += vBest[iPar];
	vRMS[iPar].back() += vBest[iPar]*vBest[iPar];
      }
    }

    for (int iPar=0; iPar<nPars; ++iPar) {
      // cout<<vRMS[iPar].back()<<", "<<vBias[iPar].back()<<"("<<vBias[iPar].back() * vBias[iPar].back()<<") -> "<<( vRMS[iPar].back() - vBias[iPar].back() * vBias[iPar].back() ) / ( fitResultsTree.GetEntries() - 1 )<<endl;
      vRMS[iPar].back() = sqrt( ( vRMS[iPar].back() - vMean[iPar].back() * vMean[iPar].back() / fitResultsTree.GetEntries() ) / ( fitResultsTree.GetEntries() - 1 ) );
      vMean[iPar].back() = vMean[iPar].back() / fitResultsTree.GetEntries();
      vBias[iPar].back() = vMean[iPar].back() - vHistBestRECO[iPar].back();
      vMeanErr[iPar].back() = vRMS[iPar].back() / sqrt( fitResultsTree.GetEntries() );

      printf("%s:\tBias (wrt full Jpsi fit result) = %.5f\tRMS deviation: %.5f\n",parName[iPar].c_str(),vBias[iPar].back(),vRMS[iPar].back());
    }

  }

  int nPlotBins = vHistBestRECO[0].size();
  if (nPlotBins<1) {
    cout<<"ERROR, no q2 bins processed!"<<endl;
    return;
  }

  gStyle->SetOptStat(0);

  vector<TCanvas*> cDistr (nPars);
  vector<TCanvas*> cUncert (nPars);
  vector<TCanvas*> cPull (nPars);
  vector<TCanvas*> cResult (nPars);
  vector< vector<TLine*> > lineRECO (nPars);
  vector< vector<TLine*> > lineRMS (nPars);

  for (int iPar=0; iPar<nPars; ++iPar) {

    cDistr[iPar] = new TCanvas(Form("cDistr%i",iPar),Form("%s distribution",parTitle[iPar].c_str()),2000,1000);
    cUncert[iPar] = new TCanvas(Form("cUncert%i",iPar),Form("%s uncertainty",parTitle[iPar].c_str()),2000,1000);
    cPull[iPar] = new TCanvas(Form("cPull%i",iPar),Form("%s pulls",parTitle[iPar].c_str()),2000,1000);
    cResult[iPar] = new TCanvas(Form("cResult%i",iPar),Form("%s results",parTitle[iPar].c_str()),1000,1000);

    cDistr[iPar]->cd();
    vHistBest[iPar][0]->Draw();
    double ymax = vHistBest[iPar][0]->GetMaximum();

    cPull[iPar]->cd();
    vHistPull[iPar][0]->Draw();
    double ymaxPull = vHistPull[iPar][0]->GetMaximum();

    cUncert[iPar]->cd()->SetLogx();
    // copy underflow and overflow in first and last bins
    vHistErrH[iPar][0]->AddBinContent(1,vHistErrH[iPar][0]->GetBinContent(0));
    vHistErrL[iPar][0]->AddBinContent(1,vHistErrL[iPar][0]->GetBinContent(0));
    vHistErrH[iPar][0]->AddBinContent(nUncBins,vHistErrH[iPar][0]->GetBinContent(nUncBins+1));
    vHistErrL[iPar][0]->AddBinContent(nUncBins,vHistErrL[iPar][0]->GetBinContent(nUncBins+1));
    vHistErrH[iPar][0]->GetXaxis()->SetMoreLogLabels();
    vHistErrH[iPar][0]->Draw();
    vHistErrL[iPar][0]->Draw("same");
    double ymaxUnc = vHistErrH[iPar][0]->GetMaximum();

    TLegend* leg;
    TLegend* legPull;
    TLegend* legUnc;
    if ( parName[iPar].compare("P4p")==0 || parName[iPar].compare("P5p")==0 || parName[iPar].compare("P1")==0 )
      leg = new TLegend(0.67,0.57,0.87,0.87,"q^{2} bin");
    else
      leg = new TLegend(0.15,0.57,0.35,0.87,"q^{2} bin");
    if ( parName[iPar].compare("P4p")==0 || parName[iPar].compare("P8p")==0 || parName[iPar].compare("P3")==0 )
      legUnc = new TLegend(0.15,0.62,0.35,0.87,"q^{2} bin");
    else
      legUnc = new TLegend(0.67,0.62,0.87,0.87,"q^{2} bin");
    legUnc->SetNColumns(2);
    legPull = new TLegend(0.67,0.57,0.87,0.87,"q^{2} bin");
    if (nPlotBins>1) {
      // vHistBest[iPar][0]->SetTitle( Form("%s results of data-like MC sample fits",parTitle[iPar].c_str()) );
      // vHistErrH[iPar][0]->SetTitle( Form("%s MINOS uncertainties of data-like MC sample fits",parTitle[iPar].c_str()) );
      leg->SetBorderSize(0);
      leg->AddEntry(vHistBest[iPar][0],Form("%i [Bias:%.3f RMS:%.3f]",vq2Bins[0],vBias[iPar][0],vRMS[iPar][0]),"l");
      legPull->SetBorderSize(0);
      legPull->AddEntry(vHistPull[iPar][0],Form("%i [Mean:%.3f RMS:%.3f]",vq2Bins[0],vHistPull[iPar][0]->GetMean(),vHistPull[iPar][0]->GetRMS()),"l");
    }

    legUnc->SetBorderSize(0);
    legUnc->AddEntry(vHistErrH[iPar][0],Form("%i higher, mean: %.3f",vq2Bins[0], vHistErrH[iPar][0]->GetMean()),"f");
    legUnc->AddEntry(vHistErrL[iPar][0],Form("%i lower, mean: %.3f",vq2Bins[0], vHistErrL[iPar][0]->GetMean()),"f");

    printf("%s:\t Bin 0 Mean (high/low) = %.5f\t\n",parName[iPar].c_str(),(vHistErrH[iPar][0]->GetMean()+vHistErrL[iPar][0]->GetMean())/2);

    for (int iBin=1; iBin<nPlotBins; ++iBin) {
      cDistr[iPar]->cd();
      vHistBest[iPar][iBin]->Draw("same");
      if (ymax < vHistBest[iPar][iBin]->GetMaximum()) ymax = vHistBest[iPar][iBin]->GetMaximum();

      cPull[iPar]->cd();
      vHistPull[iPar][iBin]->Draw("same");
      if (ymaxPull < vHistPull[iPar][iBin]->GetMaximum()) ymaxPull = vHistPull[iPar][iBin]->GetMaximum();

      cUncert[iPar]->cd();
      vHistErrH[iPar][iBin]->AddBinContent(1,vHistErrH[iPar][iBin]->GetBinContent(0));
      vHistErrL[iPar][iBin]->AddBinContent(1,vHistErrL[iPar][iBin]->GetBinContent(0));
      vHistErrH[iPar][iBin]->AddBinContent(nUncBins,vHistErrH[iPar][iBin]->GetBinContent(nUncBins+1));
      vHistErrL[iPar][iBin]->AddBinContent(nUncBins,vHistErrL[iPar][iBin]->GetBinContent(nUncBins+1));
      vHistErrH[iPar][iBin]->Draw("same");
      vHistErrL[iPar][iBin]->Draw("same");
      if (ymaxUnc < vHistErrH[iPar][iBin]->GetMaximum()) ymaxUnc = vHistErrH[iPar][iBin]->GetMaximum();

      leg->AddEntry(vHistBest[iPar][iBin],Form("%i [Bias:%.3f RMS:%.3f]",vq2Bins[iBin],vBias[iPar][iBin],vRMS[iPar][iBin]),"l");
      legPull->AddEntry(vHistPull[iPar][iBin],Form("%i [Mean:%.3f RMS:%.3f]",vq2Bins[iBin],vHistPull[iPar][iBin]->GetMean(),vHistPull[iPar][iBin]->GetRMS()),"l");
      legUnc->AddEntry(vHistErrH[iPar][iBin],Form("%i higher, mean: %.3f",vq2Bins[iBin], vHistErrH[iPar][iBin]->GetMean()),"f");
      legUnc->AddEntry(vHistErrL[iPar][iBin],Form("%i lower, mean: %.3f",vq2Bins[iBin], vHistErrL[iPar][iBin]->GetMean()),"f");

      printf("%s:\tBin 5 Mean (high/low) = %.5f\t\n",parName[iPar].c_str(),(vHistErrH[iPar][iBin]->GetMean()+vHistErrL[iPar][iBin]->GetMean())/2);
    }

    vHistBest[iPar][0]->GetYaxis()->SetRangeUser(0,1.1*ymax);
    vHistPull[iPar][0]->GetYaxis()->SetRangeUser(0,1.1*ymaxPull);
    vHistErrH[iPar][0]->GetYaxis()->SetRangeUser(0,1.1*ymaxUnc);

    for (int iBin=0; iBin<nPlotBins; ++iBin) {
      cDistr[iPar]->cd();
      lineRECO[iPar].push_back( new TLine(vHistBestRECO[iPar][iBin],0,vHistBestRECO[iPar][iBin],1.1*ymax) );
      lineRECO[iPar].back()->SetLineWidth(2);
      lineRECO[iPar].back()->SetLineColor(colors[iBin]);
      lineRECO[iPar].back()->Draw();

      cUncert[iPar]->cd();
      lineRMS[iPar].push_back( new TLine(vRMS[iPar][iBin],0,vRMS[iPar][iBin],1.1*ymaxUnc) );
      lineRMS[iPar].back()->SetLineWidth(2);
      lineRMS[iPar].back()->SetLineColor(colors[iBin]);
      // lineRMS[iPar].back()->Draw();
    }

    cDistr[iPar]->cd();
    if (nPlotBins>1) leg->Draw();
    cPull[iPar]->cd();
    if (nPlotBins>1) legPull->Draw();
    cUncert[iPar]->cd();
    legUnc->Draw();

    cDistr[iPar]->SaveAs(Form("plotSimFit_d/simfit_dataCRsub_%s_dist_wp90.pdf",parName[iPar].c_str()));
    cPull[iPar]->SaveAs(Form("plotSimFit_d/simfit_dataCRsub_%s_pulls_wp90.pdf",parName[iPar].c_str()));
    cUncert[iPar]->SaveAs(Form("plotSimFit_d/simfit_dataCRsub_%s_uncert_wp90.pdf",parName[iPar].c_str()));

  }
  
  // Plot P-wave parameters
  
  double x[nStats][nPars];
  double xe[nPars];
  double xSim[nPars];
  double xeSim[nPars];

  double ySim [nPars];
  double ylSim[nPars];
  double yhSim[nPars];
  double y [nStats][nPars];
  double yl[nStats][nPars];
  double yh[nStats][nPars];

  for (int iPar=0; iPar<nPars; ++iPar) {
    xSim[iPar]  = 0.5+iPar;
    xeSim[iPar] = 0.5;
    xe[iPar]    = 1.0/(2*nStats);
    ySim [iPar] = 0;
    ylSim[iPar] = -1 * vHistErrLRECO[iPar][0];
    yhSim[iPar] = vHistErrHRECO[iPar][0];
    for (int iStat=0; iStat<nStats; ++iStat) {
      x [iStat][iPar] = (2.0*iStat+1)/(2*nStats)+iPar;
      y [iStat][iPar] = vMean[iPar][iStat] - vHistBestRECO[iPar][0];
      yl[iStat][iPar] = yh[iStat][iPar] = vRMS[iPar][iStat];
    }
  }

  TCanvas canv("canv","canv",1500,1000);
  canv.cd();

  gStyle->SetOptStat(0);
  // cCover.cd()->SetTicky(2);
  string plotTitle = "Control-region subsample fit results;;(subsample result) - (full-sample result)";
  auto hLab = new TH1S ("hLab",plotTitle.c_str(),nPars,0,nPars);
  for (int iBin=0; iBin<nPars; ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,parName[iBin].c_str());
  double yRange = 0.25;
  hLab->SetMinimum(-1*yRange);
  hLab->SetMaximum(yRange);
  // hLab->GetYaxis()->SetTickLength(0.006);
  hLab->Draw();

  auto grSim = new TGraphAsymmErrors(nPars,xSim,ySim,xeSim,xeSim,ylSim,yhSim);
  grSim->SetName("grSim");
  grSim->SetFillColor(429);
  grSim->Draw("2");

  TLine line (0,0,nPars,0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.SetLineColor(13);
  line.Draw();

  TGraphAsymmErrors* gr [3];
  for (int iStat=0; iStat<nStats; ++iStat) {
    gr[iStat] = new TGraphAsymmErrors(nPars,x[iStat],y[iStat],xe,xe,yl[iStat],yh[iStat]);
    gr[iStat]->SetName(Form("gr%i",iStat));
    gr[iStat]->SetLineColor(colors[iStat]);
    gr[iStat]->SetLineWidth(2);
    gr[iStat]->Draw("P");
  }

  TLegend leg (0.15,0.70,0.4,0.9);
  leg.AddEntry(grSim,"Control-region fit","f");
  for (int iStat=0; iStat<nStats; ++iStat)
    leg.AddEntry(gr[iStat],Form("q2-bin %i stat",aq2stat[iStat]),"lep");
  leg.Draw();

  canv.SaveAs("plotSimFit_d/simfit_dataCRsub_results_xgbv8_wp90.pdf");

}
