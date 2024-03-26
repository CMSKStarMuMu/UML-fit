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

static const int nVariations = 1;
int aq2stat [nVariations] = { 0 };

double diffMax = 0.0499;

void computeSystematicsFromDataFit ()
{

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  vector<int> vq2Bins (0);
  vector<string> table (0);

  binsUnc[0]=0.006;
  for (int i=0; i<nUncBins; ++i)
    binsUnc[i+1] = binsUnc[i] * pow(0.4/binsUnc[0],1./nUncBins); // to reach 0.4 as largest error

  double q2Val [nBins];
  double q2Err [nBins];

  for (int i=0; i<nBins; ++i) {
    q2Val[i] = 0.5 * (binBorders[i+1]+binBorders[i]);
    q2Err[i] = 0.5 * (binBorders[i+1]-binBorders[i]);
  }
  
  for (int q2Bin = 0; q2Bin < nBins; q2Bin++){
  
    vector<double> vNominal(nPars);
    vector<double> vHigh(nPars);
    vector<double> vLow (nPars);

    vector<double> vSyst(nPars);
    vector<double> vSystH(nPars);
    vector<double> vSystL(nPars);

    if (q2Bin==4 || q2Bin==6) {
      continue;
    }

    vector<double> vBias(nPars);
    
    int q2stat = aq2stat[0];
    vq2Bins.push_back(q2stat);
  
  
    TChain fitResultsTree ("fitResultsTree","");
    string filename = Form("/afs/cern.ch/work/d/dini/public/Bin0AndCircularityStudies/FitCircularity/DataBins_SB_Circularity/simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_unbl4.root",q2Bin);
//     string filename = Form("/afs/cern.ch/work/d/dini/public/systPDG-DataBins/simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_unbl4.root",q2Bin);
    //  mass shape syst
    //  string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/massShapeSyst/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_f1.06_unbl4.root",q2Bin);
    //  mistag syst
    //  string filename = Form("simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_wp90_unbl4_mistagSyst_exp.root",q2Bin);
    fitResultsTree.Add(filename.c_str());
  
    string filename_nominal = Form("/afs/cern.ch/work/d/dini/public/Bin0AndCircularityStudies/FitCircularity/DataBins_SB_Effi_Circularity/simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_unbl4.root",q2Bin);
//     string filename_nominal = Form("simFitResults4d/simFitResult_data_fullAngularMass_Swave_201620172018_b%i-XGBv8_wp90_unbl4.root",q2Bin);
    TFile* filein_nominal = TFile::Open(filename_nominal.c_str());
    TTree* fitResultsTree_nominal = (TTree*)filein_nominal->Get("fitResultsTree");
    if (!fitResultsTree_nominal || fitResultsTree_nominal->GetEntries() != 1) {
      cout<<"Error, unexpected number of entries in fitResultsTree in file: "<<filename_nominal<<endl;
      return;
    }
  
    int nSamp = fitResultsTree.GetEntries();
    cout<<"Number of samples: "<<nSamp<<endl;
    
    for (int iPar=0; iPar<nPars; ++iPar) {
      fitResultsTree_nominal->SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vNominal[iPar]);
      fitResultsTree_nominal->SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh[iPar]);
      fitResultsTree_nominal->SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow [iPar]);

      fitResultsTree.SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vSyst[iPar]);
      fitResultsTree.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vSystH[iPar]);
      fitResultsTree.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vSystL[iPar]);
    }
  
    fitResultsTree_nominal->GetEntry(0);
    fitResultsTree.GetEntry(0);
//     for (int iPar=0; iPar<nPars; ++iPar) {
//       std::cout <<  "vNominal:  ipar:" << iPar << "  " << vNominal[iPar]<< std::endl;
//       std::cout <<  "vSyst:  ipar:" << iPar << "  " << vSyst[iPar] << " +- " << vSystH[iPar] << " " << vSystL[iPar] << std::endl;
//     }

    string firstLine = "\\textbf{$q^2$-bin}";
    string line = Form("%i",q2Bin);

      for (int iPar=0; iPar<nPars; ++iPar) {
        vBias[iPar] = vSyst[iPar] - vNominal[iPar];
        printf("%s:\tBias (wrt nominal result) = %.5f\t\n",parName[iPar].c_str(),vBias[iPar]);

        firstLine = firstLine +  " & \\textbf{$" + parTitle[iPar] + "$}";
        if (fabs(vBias[iPar])<0.0001)
  	  line = line + " & $<$ 0.0001";
//         if (fabs(vBias[iPar])<0.0005)
//   	  line = line + " & $<$ 0.001";
        else
//   	  line = line + Form(" & $\\pm %.4f$",(vBias[iPar]));
  	  line = line + Form(" & $\\pm %.4f$",fabs(vBias[iPar]));
      }
      if (table.size()==0) {
        firstLine = firstLine + " \\\\ \\hline";
        table.push_back(firstLine);
      }
      line = line + " \\\\";
      table.push_back(line);
      
    // Plot P-wave parameters
    
    double x[nPars];
    double xe[nPars];
    double xNom[nPars];
    double xeNom[nPars];
  
    double yNominal [nPars];
    double ylNominal[nPars];
    double yhNominal[nPars];
    double y [nPars];
    double yl[nPars];
    double yh[nPars];
  
    for (int iPar=0; iPar<nPars; ++iPar) {
      xNom[iPar]  = 0.5+iPar;
      xeNom[iPar] = 0.5;
      xe[iPar]    = 1.0/(2*nVariations);
      yNominal [iPar] = 0;
      ylNominal[iPar] = -1 * vLow[iPar];
      yhNominal[iPar] = vHigh[iPar];
      x [iPar] = (1)/(2*nVariations)+iPar;
      y [iPar] = vBias[iPar];
      yh[iPar] = vSystH[iPar];
      yl[iPar] = -1 *vSystL[iPar];
    }
  
    TCanvas canv("canv","canv",1500,1000);
    canv.cd();
  
    string plotTitle = Form("Bin %d;;(alternative result) - (nominal result)", q2Bin);
    auto hLab = new TH1S ("hLab",plotTitle.c_str(),nPars,0,nPars);
    for (int iBin=0; iBin<nPars; ++iBin)
      hLab->GetXaxis()->SetBinLabel(iBin+1,parName[iBin].c_str());
    double yRange = 0.25;
    hLab->SetMinimum(-1*yRange);
    hLab->SetMaximum(yRange);
    hLab->Draw();
  
    auto grNominal = new TGraphAsymmErrors(nPars,xNom,yNominal,xeNom,xeNom,ylNominal,yhNominal);
    grNominal->SetName("grNominal");
    grNominal->SetFillColor(429);
    grNominal->Draw("2");
  
    TLine linex (0,0,nPars,0);
    linex.SetLineWidth(2);
    linex.SetLineStyle(2);
    linex.SetLineColor(13);
    linex.Draw();
  
    TGraphAsymmErrors* gr= new TGraphAsymmErrors(nPars,xNom,y,xeNom,xeNom,yl,yh);
    gr->SetLineColor(1);
    gr->SetLineWidth(2);
    gr->Draw("P");

    TLegend leg (0.15,0.75,0.4,0.9);
    leg.AddEntry(grNominal,"Nominal fit","f");
    leg.AddEntry(gr,"sideband stat syst","lep");
//     leg.AddEntry(gr,"mass shape syst","lep");
//     leg.AddEntry(gr,"mistag fraction syst","lep");
    leg.Draw();
  
//     canv.SaveAs(Form("plotSyst/simfit_syst_mistag_results_xgbv8_wp90_bin%d.pdf", q2Bin));
//     canv.SaveAs(Form("plotSyst/simfit_syst_sidebandStat_results_xgbv8_wp90_bin%d.pdf", q2Bin));

  } // end n bins  

  cout<<"===== Formatted bias table ========"<<endl;
  for (int iLine=0; iLine<table.size(); ++iLine)
    cout<<table[iLine]<<endl;

}
