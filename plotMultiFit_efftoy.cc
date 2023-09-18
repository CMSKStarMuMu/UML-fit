#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

static const int nBins = 8;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};

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

double diffMax = 0.0999;
// double diffMax = 0.0499;

void plotMultiFit_efftoy (int binIndex=-1, int parity=0, int whichSamples = 3, bool ref4dFit = true)
{

  // whichSamples=0 -> plot fit results to 3D MC subsamples
  // whichSamples=1 -> plot fit results to 4D MC subsamples
  // whichSamples=2 -> plot fit results to 4D MC subsamples + toy background
  // whichSamples=3 -> plot fit results to data
  
  vector< vector<TH1D*> > vHistBest (nPars);
  vector< vector<TH1D*> > vHistErrH (nPars);
  vector< vector<TH1D*> > vHistErrL (nPars);

  vector< vector<TH1D*> > vHistBestDense (nPars);

  vector< vector<double> > vHistBestRECO (nPars);
  vector< vector<double> > vHistErrHRECO (nPars);
  vector< vector<double> > vHistErrLRECO (nPars);

  vector< vector<double> > vMean (nPars);
  vector< vector<double> > vRMS  (nPars);
  vector< vector<double> > vBias (nPars);
  vector< vector<double> > vMeanErr (nPars);

  vector<int> vq2Bins (0);

  vector<string> table (0);
  table.push_back("\\textbf{$q^2$-bin}");
  for (int iPar=0; iPar<nPars; ++iPar)
    table[0] = table[0] + " & \\textbf{$"+ parTitle[iPar]+"$}";

  binsUnc[0]=0.006;
  for (int i=0; i<nUncBins; ++i)
    binsUnc[i+1] = binsUnc[i] * pow(0.4/binsUnc[0],1./nUncBins); // to reach 0.4 as largest error

  double q2Val [nBins];
  double q2Err [nBins];

  for (int i=0; i<nBins; ++i) {
    q2Val[i] = 0.5 * (binBorders[i+1]+binBorders[i]);
    q2Err[i] = 0.5 * (binBorders[i+1]-binBorders[i]);
  }

  fstream fs ("../confSF/KDE_SF.list", fstream::in);
  string fullLine = "";

  int iColor = 0;

  do {
    getline(fs,fullLine);
    if (fullLine.empty()) break;
    cout<<"Running "<<fullLine<<endl;

    stringstream ss(fullLine);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> splitLine(begin, end);

    int q2Bin = atoi(splitLine[0].c_str());
    if (binIndex>=0 && binIndex!=q2Bin) continue;

    TChain fitResultsTree ("fitResultsTree","");
    // string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/simFitResults4d/simFitResult_recoMC_fullAngularMass_toyeff201620172018_dataStat-0_b%i_good.root",q2Bin);
    // string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/MCstat/simFitResults4d/simFitResult_recoMC_fullAngularMass_toybkg201620172018_dataStat-0_b%i_toyeff*.root",q2Bin);
    string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/MCstat/simFitResults4d/xgbv8/simFitResult_recoMC_fullAngularMass_toybkg201620172018_dataStat-0_p%i_b%i_toyeff*.root",parity,q2Bin);
    if (whichSamples==3)
      filename = Form("simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_toyeff*.root",q2Bin);
    // string filename = Form("simFitResults4d/simFitResult_recoMC_fullAngularMass_toybkg201620172018_dataStat-*_b%i.root",q2Bin);
    // if (plot4dFit==1) filename = Form("simFitResults4d/simFitResult_recoMC_fullAngularMass201620172018_dataStat-*_b%i.root",q2Bin);
    // if (plot4dFit==0) filename = Form("simFitResults/simFitResult_recoMC_fullAngular201620172018_dataStat-*_b%i.root",q2Bin);
    fitResultsTree.Add(filename.c_str());

    // string filename_fR = Form("simFitResults4d/simFitResult_recoMC_fullAngularMass_toybkg201620172018_dataStat-0_b%i.root",q2Bin);
    string filename_fR = Form("simFitResults4d/xgbv8/simFitResult_recoMC_fullAngularMass_toybkg201620172018_dataStat-0_b%i.root",q2Bin);
    if (whichSamples==3)
      filename_fR = Form("./simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_unbl4.root",q2Bin);
    // string filename_fR = Form("../UML-fit-integration/simFitResults4d/simFitResult_recoMC_fullAngularMass201620172018_dataStat-0_b%i.root",q2Bin);
    // string filename_fR = Form("simFitResults4d/simFitResult_recoMC_fullAngularMass201620172018_MCStat_b%i.root",q2Bin);
    // if (!ref4dFit) filename_fR = Form("simFitResults/simFitResult_recoMC_fullAngular201620172018_MCStat_b%i.root",q2Bin);
    TFile* filein_fR = TFile::Open(filename_fR.c_str());
    TTree* fitResultsTree_fR = (TTree*)filein_fR->Get("fitResultsTree");
    if (!fitResultsTree_fR || fitResultsTree_fR->GetEntries() != 1) {
      cout<<"Error, unexpected numebr of entries in fitResultsTree in file: "<<filename_fR<<endl;
      iColor--;
      continue;
    }

    int nSamp = fitResultsTree.GetEntries();
    cout<<"Number of samples: "<<nSamp<<endl;
    if (nSamp<2) {iColor--; continue;}

    vq2Bins.push_back(q2Bin);

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

      vHistBest[iPar].push_back( new TH1D(Form("hBest%i%i",q2Bin,iPar),Form("%s results of data fit with toy efficiencies - q2 bin %i;%s;# of results",parTitle[iPar].c_str(),q2Bin,parTitle[iPar].c_str()),100,parMin[iPar],parMax[iPar]) );
      vHistBestDense[iPar].push_back( new TH1D(Form("hBestDense%i%i",q2Bin,iPar),"",
					       (parMax[iPar]-parMin[iPar])*300,
					       parMin[iPar],parMax[iPar]) );
      vHistErrH[iPar].push_back( new TH1D(Form("hErrH%i%i",q2Bin,iPar),Form("%s MINOS uncertainties of data fit with toy efficiencies - q2 bin %i;#sigma(%s);# of results",parTitle[iPar].c_str(),q2Bin,parTitle[iPar].c_str()),nUncBins,binsUnc) );
      vHistErrL[iPar].push_back( new TH1D(Form("hErrL%i%i",q2Bin,iPar),Form("%s MINOS uncertainties of data fit with toy efficiencies - q2 bin %i;#sigma(%s);# of results",parTitle[iPar].c_str(),q2Bin,parTitle[iPar].c_str()),nUncBins,binsUnc) );

      vHistBest[iPar].back()->SetLineColor(colors[iColor]);
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
      vHistErrLRECO[iPar].push_back( (whichSamples==3?-1:1)*vLow [iPar] );
    }

    // Pre-cycle
    for (int iEn=0; iEn<fitResultsTree.GetEntries(); ++iEn) {
      fitResultsTree.GetEntry(iEn);
      for (int iPar=0; iPar<nPars; ++iPar)
	vHistBestDense[iPar].back()->Fill(vBest[iPar]);
    }

    // Identification of a ~10 sigma range in which to evaluate the RMS
    // results of problematic fits that fall outside it will not be used
    double p[3]={0.16,0.5,0.84};
    double q[3];
    vector<double> vRangeMin (nPars);
    vector<double> vRangeMax (nPars);
    int nSigma = 10;
    for (int iPar=0; iPar<nPars; ++iPar) {
      vHistBestDense[iPar].back()->GetQuantiles(3,q,p);
      vRangeMin[iPar] = nSigma*q[0] - (nSigma-1)*q[1];
      vRangeMax[iPar] = nSigma*q[2] - (nSigma-1)*q[1];
      // cout<<vRangeMin[iPar]<<"\t"<<vRangeMax[iPar]<<endl;
    }

    vector<int> vNevents (nPars,0);
    for (int iEn=0; iEn<fitResultsTree.GetEntries(); ++iEn) {
      fitResultsTree.GetEntry(iEn);
      for (int iPar=0; iPar<nPars; ++iPar) {
	if (vBest[iPar]<vRangeMin[iPar] || vBest[iPar]>vRangeMax[iPar]) continue;
	vHistBest[iPar].back()->Fill(vBest[iPar]);
	vHistErrH[iPar].back()->Fill(vHigh[iPar]-vBest[iPar]);
	vHistErrH[iPar].back()->Fill(vBest[iPar]-vLow [iPar]); // To create a stacked histogram
	vHistErrL[iPar].back()->Fill(vBest[iPar]-vLow [iPar]);

	vMean[iPar].back() += vBest[iPar];
	vRMS[iPar].back() += vBest[iPar]*vBest[iPar];
	++vNevents[iPar];
      }
    }

    string line = Form("%i",q2Bin);
    string firstLine = "\\textbf{$q^2$-bin}";

    for (int iPar=0; iPar<nPars; ++iPar) {
      // cout<<vRMS[iPar].back()<<", "<<vBias[iPar].back()<<"("<<vBias[iPar].back() * vBias[iPar].back()<<") -> "<<( vRMS[iPar].back() - vBias[iPar].back() * vBias[iPar].back() ) / ( fitResultsTree.GetEntries() - 1 )<<endl;
      vRMS[iPar].back() = sqrt( ( vRMS[iPar].back() - vMean[iPar].back() * vMean[iPar].back() / vNevents[iPar] ) / ( vNevents[iPar] - 1 ) );
      vMean[iPar].back() = vMean[iPar].back() / vNevents[iPar];
      vBias[iPar].back() = vMean[iPar].back() - vHistBestRECO[iPar].back();
      vMeanErr[iPar].back() = vRMS[iPar].back() / sqrt( vNevents[iPar] );

      printf("%s:\tBias (wrt RECO result) = %.5f\tRMS deviation: %.5f\n",parName[iPar].c_str(),vBias[iPar].back(),vRMS[iPar].back());

      line = line + Form(" & %.3f",vRMS[iPar].back());

    }

    table.push_back(line);

  } while (++iColor);

  int nPlotBins = vHistBestRECO[0].size();
  if (nPlotBins<1) {
    cout<<"ERROR, no q2 bins processed!"<<endl;
    return;
  }

  gStyle->SetOptStat(0);

  vector<TCanvas*> cDistr (nPars);
  vector<TCanvas*> cUncert (nPars);
  vector<TCanvas*> cResult (nPars);
  vector< vector<TLine*> > lineRECO (nPars);
  vector< vector<TLine*> > lineRMS (nPars);

  for (int iPar=0; iPar<nPars; ++iPar) {

    cDistr[iPar] = new TCanvas(Form("cDistr%i",iPar),Form("%s distribution",parTitle[iPar].c_str()),2000,1000);
    cUncert[iPar] = new TCanvas(Form("cUncert%i",iPar),Form("%s uncertainty",parTitle[iPar].c_str()),2000,1000);
    cResult[iPar] = new TCanvas(Form("cResult%i",iPar),Form("%s results",parTitle[iPar].c_str()),1000,1000);

    cDistr[iPar]->cd();
    vHistBest[iPar][0]->Draw();
    double ymax = vHistBest[iPar][0]->GetMaximum();

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
    TLegend* legUnc;
    if ( parName[iPar].compare("P4p")==0 || parName[iPar].compare("P5p")==0 || parName[iPar].compare("P1")==0 )
      leg = new TLegend(0.67,0.57,0.87,0.87,"q^{2} bin");
    else
      leg = new TLegend(0.15,0.57,0.35,0.87,"q^{2} bin");
    if ( parName[iPar].compare("P4p")==0 || parName[iPar].compare("P8p")==0 || parName[iPar].compare("P3")==0 )
      legUnc = new TLegend(0.15,0.57,0.35,0.87,"q^{2} bin");
    else
      legUnc = new TLegend(0.67,0.57,0.87,0.87,"q^{2} bin");
    legUnc->SetNColumns(2);
    if (nPlotBins>1) {
      vHistBest[iPar][0]->SetTitle( Form("%s results of data fit with toy efficiencies",parTitle[iPar].c_str()) );
      vHistErrH[iPar][0]->SetTitle( Form("%s MINOS uncertainties of data fit with toy efficiencies",parTitle[iPar].c_str()) );
      leg->SetBorderSize(0);
      leg->AddEntry(vHistBest[iPar][0],Form("%i [Bias:%.3f RMS:%.3f]",vq2Bins[0],vBias[iPar][0],vRMS[iPar][0]),"l");
    }

    legUnc->SetBorderSize(0);
    legUnc->AddEntry(vHistErrH[iPar][0],Form("%i higher",vq2Bins[0]),"f");
    legUnc->AddEntry(vHistErrL[iPar][0],Form("%i lower",vq2Bins[0]),"f");

    for (int iBin=1; iBin<nPlotBins; ++iBin) {
      cDistr[iPar]->cd();
      vHistBest[iPar][iBin]->Draw("same");
      if (ymax < vHistBest[iPar][iBin]->GetMaximum()) ymax = vHistBest[iPar][iBin]->GetMaximum();

      cUncert[iPar]->cd();
      vHistErrH[iPar][iBin]->AddBinContent(1,vHistErrH[iPar][iBin]->GetBinContent(0));
      vHistErrL[iPar][iBin]->AddBinContent(1,vHistErrL[iPar][iBin]->GetBinContent(0));
      vHistErrH[iPar][iBin]->AddBinContent(nUncBins,vHistErrH[iPar][iBin]->GetBinContent(nUncBins+1));
      vHistErrL[iPar][iBin]->AddBinContent(nUncBins,vHistErrL[iPar][iBin]->GetBinContent(nUncBins+1));
      vHistErrH[iPar][iBin]->Draw("same");
      vHistErrL[iPar][iBin]->Draw("same");
      if (ymaxUnc < vHistErrH[iPar][iBin]->GetMaximum()) ymaxUnc = vHistErrH[iPar][iBin]->GetMaximum();

      leg   ->AddEntry(vHistBest[iPar][iBin],Form("%i [Bias:%.3f RMS:%.3f]",vq2Bins[iBin],vBias[iPar][iBin],vRMS[iPar][iBin]),"l");
      legUnc->AddEntry(vHistErrH[iPar][iBin],Form("%i higher",vq2Bins[iBin]),"f");
      legUnc->AddEntry(vHistErrL[iPar][iBin],Form("%i lower",vq2Bins[iBin]),"f");
    }

    vHistBest[iPar][0]->GetYaxis()->SetRangeUser(0,1.1*ymax);
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
      lineRMS[iPar].back()->Draw();
    }

    cDistr[iPar]->cd();
    if (nPlotBins>1) leg->Draw();
    cUncert[iPar]->cd();
    legUnc->Draw();

    string confString = "recoMC_toybkg";
    if (whichSamples==3)
      confString = "data";
    // if (plot4dFit) toyConfString = toyConfString + "_4Dbkg";
    // if (!ref4dFit) toyConfString = toyConfString + "_vs3DfullMC";

    cDistr[iPar]->SaveAs(Form("plotSimFit4d_d/xgbv8/simfit_%s_efftoy_%s_dist_XGBv8.pdf",confString.c_str(),parName[iPar].c_str()));
    cUncert[iPar]->SaveAs(Form("plotSimFit4d_d/xgbv8/simfit_%s_efftoy_%s_uncert_XGBv8.pdf",confString.c_str(),parName[iPar].c_str()));

    // Plot resutls vs q2

    double aMean [nBins];
    double aMeanErr [nBins];
    double aReco [nBins];
    double aRecoErrH [nBins];
    double aRecoErrL [nBins];
    double aBias [nBins];
    double aBiasErrH [nBins];
    double aBiasErrL [nBins];
    double aQuantInnerCenter [nBins];
    double aQuantOuterCenter [nBins];
    double aQuantInnerError [nBins];
    double aQuantOuterError [nBins];

    if (nQuant<4) {
      cout<<"Too few quantile values provided: "<<nQuant<<endl;
      return;
    }

    for (int iBin=0; iBin<nBins; ++iBin) {
     aMean[iBin] = -9;
     aMeanErr[iBin] = 1;

     aReco[iBin] = -9;
     aRecoErrH[iBin] = 1;
     aRecoErrL[iBin] = 1;

     aBias[iBin] = -9;
     aBiasErrH[iBin] = 1;
     aBiasErrL[iBin] = 1;

     aQuantInnerCenter[iBin] = -9;
     aQuantOuterCenter[iBin] = -9;
     aQuantInnerError[iBin] = 1;
     aQuantOuterError[iBin] = 1;
    }
    for (int iBin=0; iBin<nPlotBins; ++iBin) {
      aMean[vq2Bins[iBin]] = vMean[iPar][iBin];
      aMeanErr[vq2Bins[iBin]] = vMeanErr[iPar][iBin];

      aReco[vq2Bins[iBin]] = vHistBestRECO[iPar][iBin];
      aRecoErrH[vq2Bins[iBin]] = vHistErrHRECO[iPar][iBin];
      aRecoErrL[vq2Bins[iBin]] = vHistErrLRECO[iPar][iBin];

      aBias[vq2Bins[iBin]] = vBias[iPar][iBin];
      aBiasErrH[vq2Bins[iBin]] = sqrt( vMeanErr[iPar][iBin]*vMeanErr[iPar][iBin] + vHistErrLRECO[iPar][iBin]*vHistErrLRECO[iPar][iBin] );
      aBiasErrL[vq2Bins[iBin]] = sqrt( vMeanErr[iPar][iBin]*vMeanErr[iPar][iBin] + vHistErrHRECO[iPar][iBin]*vHistErrHRECO[iPar][iBin] );

      double quantVal[nQuant];
      vHistBest[iPar][iBin]->GetQuantiles(nQuant,quantVal,quantPerc);
      aQuantInnerCenter[vq2Bins[iBin]] = 0.5 * ( quantVal[1] + quantVal[2] );
      aQuantOuterCenter[vq2Bins[iBin]] = 0.5 * ( quantVal[0] + quantVal[3] );
      aQuantInnerError[vq2Bins[iBin]] = 0.5 * fabs( quantVal[1] - quantVal[2] );
      aQuantOuterError[vq2Bins[iBin]] = 0.5 * fabs( quantVal[0] - quantVal[3] );
    }


    auto GrReco = new TGraphAsymmErrors(nBins, q2Val, aReco, q2Err, q2Err, aRecoErrL, aRecoErrH);
    GrReco->SetName(Form("GrReco%i",iPar));
    GrReco->SetTitle(Form("%s results from data fit with toy efficiencies",parTitle[iPar].c_str()));
    GrReco->GetYaxis()->SetTitle(parTitle[iPar].c_str());
    auto Gr = new TGraphErrors(nBins, q2Val, aMean, q2Err, aMeanErr);
    Gr->SetName(Form("Gr%i",iPar));    
    auto GrDiff = new TGraphAsymmErrors(nBins, q2Val, aBias, q2Err, q2Err, aBiasErrL, aBiasErrH);
    GrDiff->SetName(Form("GrDiff%i",iPar));
    auto GrQuantIn = new TGraphErrors(nBins, q2Val, aQuantInnerCenter, q2Err, aQuantInnerError);
    GrQuantIn->SetName(Form("GrQuantIn%i",iPar));    
    auto GrQuantOut = new TGraphErrors(nBins, q2Val, aQuantOuterCenter, q2Err, aQuantOuterError);
    GrQuantOut->SetName(Form("GrQuantOut%i",iPar));    

    Gr->SetLineColor(kRed+1);
    Gr->SetMarkerColor(kRed+1);
    GrReco->SetLineColor(1);
    GrReco->SetMarkerColor(1);
    GrDiff->SetLineColor(1);
    GrDiff->SetMarkerColor(1);
    Gr->SetLineWidth(2);
    GrReco->SetLineWidth(2);
    GrDiff->SetLineWidth(2);

    GrQuantIn->SetFillColor(38);
    GrQuantOut->SetFillColor(38);
    GrQuantOut->SetLineColor(38);
    GrQuantOut->SetLineWidth(2);
    GrQuantIn->SetFillStyle(3345);
    GrQuantOut->SetFillStyle(0);

    // Grey bands for resonant regions
    double ResX [2] = {0.5*(binBorders[5]+binBorders[4]),0.5*(binBorders[7]+binBorders[6])};
    double ResXe[2] = {0.5*(binBorders[5]-binBorders[4]),0.5*(binBorders[7]-binBorders[6])};
    double ResY [2] = {0.5*(parMax[iPar]+parMin[iPar]),0.5*(parMax[iPar]+parMin[iPar])};
    double ResYe[2] = {0.498*(parMax[iPar]-parMin[iPar]),0.498*(parMax[iPar]-parMin[iPar])};
    double ResD [2] = {0,0};
    double ResDe[2] = { 0.98*diffMax, 0.98*diffMax};
    TGraphErrors *resCover = new TGraphErrors(2,ResX,ResY,ResXe,ResYe);
    resCover->SetName(Form("resCover%i",iPar));
    resCover->SetFillColor(18);
    resCover->SetFillStyle(1001);
    TGraphErrors *resDiffCover = new TGraphErrors(2,ResX,ResD,ResXe,ResDe);
    resDiffCover->SetName(Form("resDiffCover%i",iPar));
    resDiffCover->SetFillColor(18);
    resDiffCover->SetFillStyle(1001);

    // Legend
    TLegend *legRes;
    if (iPar==4 || iPar==5) legRes = new TLegend(0.15,0.65,0.4,0.85);
    // else if (iPar==2) legRes = new TLegend(0.48,0.1,0.9,0.3);
    else legRes = new TLegend(0.4,0.05,0.9,0.25);
    legRes->SetName(Form("legRes%i",iPar));
    legRes->SetBorderSize(0);
    legRes->SetFillColor(1);
    legRes->SetFillStyle(0);
    legRes->SetTextSize(0.032);
    legRes->AddEntry(GrReco,"Fit with nominal efficiency","lep");
    legRes->AddEntry(Gr,"Mean of toy-efficiency fit results","lep");
    legRes->AddEntry(GrQuantIn,"Central 68\% of the toy-efficiency results","f");
    legRes->AddEntry(GrQuantOut,"Central 95\% of the toy-efficiency results","f");

    // Zero line
    TLine *line = new TLine(GrReco->GetXaxis()->GetXmin(),0,GrReco->GetXaxis()->GetXmax(),0);
    line->SetLineColor(14);
    line->SetLineStyle(7);

    cResult[iPar]->cd();
    TPad *pad1 = new TPad(Form("pad1_%i",iPar), "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
  
     
    GrReco->Draw("AP");

    GrReco->GetYaxis()->SetLabelSize(0.);
    GrReco->GetYaxis()->SetTitleSize(0.);
    GrReco->GetYaxis()->SetRangeUser(parMin[iPar],parMax[iPar]);
    TGaxis *axis = new TGaxis( GrReco->GetXaxis()->GetXmin(), parMin[iPar]+0.01,
    			       GrReco->GetXaxis()->GetXmin(), parMax[iPar], 
    			       parMin[iPar]+0.01,parMax[iPar], 
    			       510, "");
    axis->SetName(Form("axis%i",iPar));
    axis->SetTitle(parTitle[iPar].c_str());
    axis->SetLabelFont(43);
    axis->SetLabelSize(20);
    axis->Draw();

    resCover->Draw("e2");
    GrQuantOut->Draw("5");
    GrQuantIn->Draw("2");
    GrReco->Draw("P");
    Gr-> Draw("P");    
    legRes->Draw();

    // plot difference wrt RECO results
    cResult[iPar]->cd();
    TPad *pad2 = new TPad(Form("pad2_%i",iPar), "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();
    pad2->cd();

    // first create axis
    TH1F* auxE2 = new TH1F(Form("auxE2%i",iPar), "", nBins, GrReco->GetXaxis()->GetXmin(), GrReco->GetXaxis()->GetXmax());
    auxE2->SetStats(kFALSE);
    auxE2->SetLineColor(1);
    auxE2->GetXaxis()->SetTitle("q^{2} (GeV^{2})");
    auxE2->GetXaxis()->SetTitleSize(0.12);
    auxE2->GetXaxis()->SetTitleOffset(0.95);
    auxE2->GetXaxis()->SetLabelSize( 0.10);
    auxE2->GetXaxis()->SetTickLength(0.1);
    auxE2->GetYaxis()->SetTitle("Bias");
    auxE2->GetYaxis()->SetTitleSize(0.10);
    auxE2->GetYaxis()->SetTitleOffset(0.45);
    auxE2->GetYaxis()->SetRangeUser(-1*diffMax,diffMax);
    auxE2->GetYaxis()->SetLabelSize(0.07);
    auxE2->GetYaxis()->SetNdivisions(505);
    auxE2->Draw();

    // Bin lines
    std::vector<TLine*> lines;
    for (int i=0; i<nBins-1; ++i) {
      lines.push_back( new TLine(binBorders[i+1], -1*diffMax, binBorders[i+1], diffMax));
      lines[i]->SetLineStyle(3);
      lines[i]->SetLineColor(kGray);
      lines[i]->Draw();
    }
  
    GrDiff -> Draw("P");    
    line->Draw();
    resDiffCover->Draw("e2");

    cResult[iPar]->SaveAs(Form("plotSimFit4d_d/xgbv8/simfit_%s_efftoy_%s_results_XGBv8.pdf",confString.c_str(),parName[iPar].c_str()));

  }

  cout<<"===== Formatted RMS table ========"<<endl;
  for (int iLine=0; iLine<table.size(); ++iLine)
    cout<<table[iLine]<<endl;

}
