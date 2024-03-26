#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>
#include <TMatrixT.h>
#include <TMatrixTSym.h>

#include <RooAbsPdf.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooAddition.h>
#include <RooRandom.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooBernstein.h"

#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"
#include "PdfSigAngMass.h"
#include "ShapeSigAng.h"

#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"
#include "Fitter.h"
#include "RooBernsteinSideband.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [4*nBins];

double power = 1.0;
bool evalMistagSys = false;


void simfit_data_fullAngularMass_SwaveBin(int q2Bin, int parity, bool multiSample, uint nSample, uint q2stat, int fitOption, int XGBv, int unblind, bool localFiles, bool plot, int save, std::vector<int> years)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string effCString = Form("effCHist_b%ip%i",q2Bin,parity);
  string effWString = Form("effWHist_b%ip%i",q2Bin,parity);
  string intCHistString = "MCint_"+shortString + "t1";
  string intWHistString = "MCint_"+shortString + "t0";
  string all_years = "";
  string year = ""; 
  string isample = ""; 
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;
  string stat = "";
  if (nSample>0)   stat = Form("_b%istat-%i",q2stat,firstSample);
  if (multiSample) stat = stat + Form("-%i",lastSample);

  string sigpdfname = "PDF_sig_ang_mass_"+shortString+"_%i";
  string bkgpdfname = "bkg_pdf_%i";

  string XGBstr = "";
  if (XGBv>9) XGBstr = Form("-TMVAv%i",XGBv-10);
  else if (XGBv>0) XGBstr = Form("-XGBv%i",XGBv);

  std::vector<TFile*> fin_eff;
  std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_sb, wsp_Z4430;
  std::vector<std::vector<RooDataSet*>> data;
  std::vector<RooAbsReal*> effC, effW;
  std::vector<TH3D*> effCHist, effWHist;
  std::vector<TH1D*> intCHist, intWHist;
  std::vector< std::vector<double> > intCVec(years.size(), std::vector<double>(0));
  std::vector< std::vector<double> > intWVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_mass_bkg(0);
  std::vector<RooAbsPdf*> PDF_sig_ang_mass_bkg_penalty(0);
  std::vector<RooGaussian*> c_deltaPeaks, c_fm;
  RooArgSet c_vars_rt, c_pdfs_rt;
  RooArgSet c_vars_wt, c_pdfs_wt;
  RooArgSet c_vars; 

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* ctK = new RooRealVar("ctK", "cos(#theta_{K})", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "cos(#theta_{l})", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "#phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6,"GeV");
  RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
  RooArgSet reco_vars (*ctK, *ctL, *phi, *mass, *rand);
  RooArgSet observables (*ctK, *ctL, *phi, *mass);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl    = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);   
  RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p   = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p   = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p   = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p   = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  RooRealVar* Fs  = new RooRealVar("Fs","F_{S}",0.05,0,1);
  RooRealVar* As  = new RooRealVar("As","A^{S}",0,-1,1);
  RooRealVar* A4s = new RooRealVar("A4s","A_{4}^{S}",0,-1,1);
  RooRealVar* A5s = new RooRealVar("A5s","A_{5}^{S}",0,-1,1);
  RooRealVar* A7s = new RooRealVar("A7s","A_{7}^{S}",0,-1,1);
  RooRealVar* A8s = new RooRealVar("A8s","A_{8}^{S}",0,-1,1);

  // only used if Z is considered
  RooRealVar* AZ  = new RooRealVar("AZ","A_{Z}",0.05,0.,1.);

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    for (uint is = firstSample; is <= lastSample; is++) {
      isample.clear(); isample.assign( Form("%i",is) );
      sample.defineType(("data"+year+"_subs"+isample).c_str());
    }
  } 

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);
  RooSimultaneous* simPdf_penalty = new RooSimultaneous("simPdf_penalty", "simultaneous pdf with penalty term", sample);

  // Define boundary check (returning 0 in physical region and 1 outside)
  BoundCheck* boundary = new BoundCheck("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // Define boundary distance calculator
  BoundDist* bound_dist = new BoundDist("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,true,0,false);

  // Define penalty term (parameters set to zero and will be set sample-by-sample)
  Penalty* penTerm = new Penalty("penTerm","Penalty term",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,0,0,0,0);

  // Random generators
  RooRandom::randomGenerator()->SetSeed(1);

  //
  // retrieve Z4430 model workspace from file
  // only a single year sample of Z is available, therefore the same parametrisation is used for the 3 years of data taking
  //
  RooAbsPdf* Z4430_ang_pdf = 0;
  RooArgSet* Z4430_ang_params = 0;
  RooAbsPdf* Z4430_mass_pdf = 0;
  RooArgSet* Z4430_mass_params = 0;
  if (q2Bin==6 && fitOption>0){
    string filename_Z4430 = "HistZ4430.root";
    if(!localFiles) filename_Z4430 = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/Zmodel/"+filename_Z4430;
    if (!(retrieveWorkspace(filename_Z4430, wsp_Z4430, "wZ4430"))) return;

    // retrieve angular component of the Z pdf
    Z4430_ang_pdf = (RooAbsPdf*) wsp_Z4430[0]->pdf("Z4430_ang_pdf");
    if(!Z4430_ang_pdf){
      std::cout << "Z4430_ang_pdf not found!!!\n" << std::endl;
      wsp_Z4430[0]->allPdfs().Print();
      exit(1);
    }
    Z4430_ang_params = (RooArgSet*) Z4430_ang_pdf->getParameters(RooArgSet(*ctK,*ctL,*phi));
    auto iter_z = Z4430_ang_params->createIterator();
    RooRealVar* ivar_z = (RooRealVar*)iter_z->Next();
    int iii=1;
    while (ivar_z) {
      ivar_z->setConstant(true);
      cout<<Form("ang par. %d) %s = %f\n",iii,ivar_z->GetName(),ivar_z->getVal())<<endl;
      ivar_z = (RooRealVar*) iter_z->Next();
      iii++;
    }
    // retrieve mass component of the Z pdf
    Z4430_mass_pdf = (RooAbsPdf*) wsp_Z4430[0]->pdf("Z4430_mass_pdf");
    if (!Z4430_mass_pdf){
      std::cout << "Z4430_mass_pdf not found!!!\n" << std::endl;
      wsp_Z4430[0]->allPdfs().Print();
      exit(1);
    }
    Z4430_mass_params = (RooArgSet*) Z4430_mass_pdf->getParameters(RooArgSet(*mass));
    auto iter_z_mass = Z4430_mass_params->createIterator();
    RooRealVar* ivar_z_mass =  (RooRealVar*)iter_z_mass->Next();
    iii=1;
    while (ivar_z_mass) {
      ivar_z_mass->setConstant(true);
      if (iii==3 && fitOption==2) ivar_z_mass->setConstant(false);
      cout << Form("mass par %d) %s = %f\n", iii, ivar_z_mass->GetName(), ivar_z_mass->getVal()) << endl;
      ivar_z_mass = (RooRealVar*) iter_z_mass->Next();
      iii++;
    }
  }

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("recoDATADataset_b%i_%i.root", q2Bin, years[iy]);
    if (!localFiles) filename_data = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/data-datasets/" + filename_data;

    // import data (or MC as data proxy)
    if (!(retrieveWorkspace(filename_data, wsp, Form("ws_b%ip0", q2Bin )))) return;

    // import KDE efficiency histograms and partial integral histograms
    string filename = Form((parity==0 ? "KDEeff_b%i_ev_%i.root" : "KDEeff_b%i_od_%i.root"),q2Bin,years[iy]);
    if (!localFiles) {
      if (XGBv<1) filename = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/eff/" + filename;
      else filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/eff-XGBv%i/",XGBv) + filename;
    }
    fin_eff.push_back( new TFile( filename.c_str(), "READ" ));
    if ( !fin_eff[iy] || !fin_eff[iy]->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;   
    }

    effCHist.push_back( (TH3D*)fin_eff[iy]->Get(effCString.c_str()));
    effWHist.push_back( (TH3D*)fin_eff[iy]->Get(effWString.c_str()));
    if ( !effCHist[iy] || effCHist[iy]->IsZombie() || !effWHist[iy] || effWHist[iy]->IsZombie() ) {
      cout<<"Efficiency histogram "<< effCString <<" or " << effWString << " not found in file: "<< filename <<endl;
      return;
    }else{
      cout<<"Efficiency histogram "<< effCString <<" or " << effWString << " found in file: "<< filename <<endl;
    }

    // create efficiency functions
    RooDataHist* effCData = new RooDataHist(("effCData_"+shortString+"_"+year).c_str(),"effCData",vars,effCHist[iy]);
    RooDataHist* effWData = new RooDataHist(("effWData_"+shortString+"_"+year).c_str(),"effWData",vars,effWHist[iy]);
    effC.push_back( new RooHistFunc(("effC_"+shortString+"_"+year).c_str(),
                                    ("effC"+year).c_str() ,
                                    vars,
                                    *effCData,
                                    1));
    effW.push_back( new RooHistFunc(("effW_"+shortString+"_"+year).c_str(),
                                    ("effW"+year).c_str() ,
                                    vars,
                                    *effWData,
                                    1));

    // import precomputed integrals and fill a std::vector
    intCHist.push_back( (TH1D*)fin_eff[iy]->Get(intCHistString.c_str()));
    intWHist.push_back( (TH1D*)fin_eff[iy]->Get(intWHistString.c_str()));
    intCVec.push_back (vector<double> (0));
    intWVec.push_back (vector<double> (0));
    if ( !intCHist[iy] || intCHist[iy]->IsZombie() || !intWHist[iy] || intWHist[iy]->IsZombie() ) {
      cout << "Integral histogram " << intCHistString <<" or " << intWHistString << " not found in file: "<< filename << endl << "Abort" << endl;
      return;
    } else if ( strcmp( intCHist[iy]->GetTitle(), effCHist[iy]->GetTitle() ) || strcmp( intWHist[iy]->GetTitle(), effWHist[iy]->GetTitle() )) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
      cout << "Integral histograms are incoherent with efficiency in file: " << filename << endl;
      cout << "Efficiency (CT) conf: " << effCHist[iy]->GetTitle() <<endl;
      cout << "Integral (CT) conf: "   << intCHist[iy]->GetTitle() <<endl;
      cout << "Efficiency (WT) conf: " << effWHist[iy]->GetTitle() <<endl;
      cout << "Integral (WT) conf: "   << intWHist[iy]->GetTitle() <<endl;
      cout << "Abort"<<endl;
      return;
    } 
    else {
      for (int i=1; i<=intCHist[iy]->GetNbinsX(); ++i) {
        intCVec[iy].push_back(intCHist[iy]->GetBinContent(i));
      }
      for (int i=1; i<=intWHist[iy]->GetNbinsX(); ++i) {
        intWVec[iy].push_back(intWHist[iy]->GetBinContent(i));
      }
    }


    // create roodataset (in case data-like option is selected, only import the correct % of data)
    if (nSample==0)
      data.push_back( createDatasetInData( wsp[iy],  q2Bin,  observables,  shortString  ));
    else
      data.push_back( createDatasetInData( nSample,  firstSample,  lastSample, wsp[iy],
					   q2Bin,  years[iy], observables,  shortString, q2stat ));



    //// Signal part of the pdf ////
    // Mass Component - import mass PDF from fits to the MC
    string filename_mc_mass = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-%s%s/%i.root",q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"LMNR"),XGBstr.c_str(), years[iy]);
    if (localFiles)
      filename_mc_mass = Form("results_fits_%i_fM_%s%s.root",years[iy], q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"lmnr"), XGBv>0 ? Form("_MCw_xgbv%i",XGBv):"");
    if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w")) {
      cout<<"Workspace not found in mass-fit file: "<<filename_mc_mass<<endl;
      return;
    }

    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    RooRealVar* mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
    RooRealVar* sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
    RooRealVar* alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.01,  100.);

    RooAbsPdf* dcb_rt;
    RooRealVar* alpha_rt2 = new RooRealVar (Form("#alpha_{RT2}^{%i}",years[iy] ), "alphart2"    , 0,    -10,   10 );
    RooRealVar* n_rt2     = new RooRealVar (Form("n_{RT2}^{%i}",years[iy])      , "nrt2"        , 0.01,   0.01,  100.);
    RooRealVar* sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,   0 , 0,   0.12, "GeV");
    RooRealVar* f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   0 , 0.,  1.);
    if (q2Bin != 7){
      alpha_rt2 -> setVal(wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal() );
      n_rt2     -> setVal(wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin)     )->getVal() );
    }

    if (q2Bin >= 4){
      sigma_rt2-> setVal(wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal() );
      f1rt     -> setVal(wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal() );
      if (q2Bin < 7) {
        if (nSample==0)  dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
        else             dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
      } else {
        if (nSample==0)  dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
        else             dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
      }
    } 
    else{
        alpha_rt2->setRange(0,10);
        if (nSample==0) dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt  );
        else dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
    }
    
    // create constrained PDF for RT mass
    RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
    constr_rt_list.add(*dcb_rt);
    RooProdPdf * c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list );
    c_vars.add(c_vars_rt);
   
    /// create WT component
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));

    RooRealVar* sigma_wt    = new RooRealVar (Form("#sigma_{WT1}^{%i}",years[iy])   , "sigmawt"    ,  wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal() ,      0,    1, "GeV");
    RooRealVar* alpha_wt1   = new RooRealVar (Form("#alpha_{WT1}^{%i}",years[iy] )  , "alphawt1"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* alpha_wt2   = new RooRealVar (Form("#alpha_{WT2}^{%i}",years[iy] )  , "alphawt2"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);
    RooRealVar* n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);

    double mean_wt_val = wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal(); 
    double mean_wt_err = wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getError(); 

    RooAbsPdf* dcb_wt;
    double deltaPeakValue = mean_rt->getVal()-mean_wt_val;
    double mean_rt_err = wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}", q2Bin))->getError(); 
    double deltaPeakError = sqrt(mean_rt_err*mean_rt_err+mean_wt_err*mean_wt_err); 

    RooRealVar* deltaPeakVar = new RooRealVar ( Form("deltaPeakVar^{%i}", years[iy]),Form("deltaPeakVar^{%i}", years[iy]), deltaPeakValue, 0., 0.2) ;
    RooGaussian* c_deltaPeaks = new RooGaussian(Form("deltaPeaks^{%i}"  , years[iy]) , "c_deltaPeaks", *deltaPeakVar, RooConst( deltaPeakValue ), RooConst(deltaPeakError )); // value to be checked
    RooFormulaVar*mWT_data = new  RooFormulaVar(Form("mWT_data^{%i}",years[iy]), "@0 + @1", RooArgList(*mean_rt, *deltaPeakVar));
    RooRealVar* mean_wt = (RooRealVar*) mWT_data;
    cout << "deltaPeak constraint built, of value " << deltaPeakValue << " +/- " << deltaPeakError << endl;

    c_pdfs_wt.add(*c_deltaPeaks);
    c_vars_wt.add(*deltaPeakVar);

    if (nSample==0) dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );
    else dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt, q2stat );

    // create constrained PDF for WT mass
    RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
    constr_wt_list.add(*dcb_wt);
    RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
    c_vars.add(c_vars_wt);

    // As the signal PDF is written as [ CT + mFrac * MT ] (see the PdfSigAngMass class),
    // the mFrac parameter represents the fraction between mistagged and correctly-tagged events (n_MT/n_CT)
    // Also, since the integral of the efficiencies contains the information on the mistag fraction in MC
    // this parameter represents the ratio between the fitted mFrac and the MC-based one ( n_MT_data / n_CT_data * n_CT_MC / n_MT_MC )
    RooRealVar* mFrac = new RooRealVar(Form("f_{M}^{%i}",years[iy]),"mistag fraction",1, 0, 15);

    // The values of fM_sigmas are computed on data-like MC subsamples as the fluctuation of the fraction of mistagged events ( n_MT/(n_MT+n_CT) )
    // this fluctuation needs to be propagated on the quantity of mFrac, to define a Gaussian contraint
    double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
    double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
    double fraction = nwt_mc / (nrt_mc + nwt_mc);
    // Propagation: sigma(mFrac) = sigma(n_MT/n_CT) * n_CT/n_MT
    //                           = sigma(fraction)/(1-fraction)^2 * (1-fraction)/fraction
    //                           = sigma(fraction) / fraction / (1-fraction)
    double frac_sigma = fM_sigmas[years[iy]][q2Bin]/fraction/(1-fraction);
    if (nSample!=0) frac_sigma = fM_sigmas[years[iy]][q2stat]/fraction/(1-fraction);

    double c_fm_val = 1;
    if (evalMistagSys){
      if (iy==0) c_fm_val = 1.07;
      else if (iy==1) c_fm_val = 1.02;
      else if (iy==2) c_fm_val = 1.03;
    }  
    RooGaussian* c_fm = new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                        RooConst(c_fm_val) , 
                                        RooConst(frac_sigma)
                                        );

    cout << "mFrac = " << fraction << " +/- " << fM_sigmas[years[iy]][q2Bin] << " ---> R = 1 +/- " << frac_sigma << endl;
    c_vars.add(*mFrac); 

    // Angular Component
    RooAbsReal* ang_rt = new ShapeSigAng(("PDF_sig_ang_rt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_rt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
					 *Fs,*As,*A4s,*A5s,*A7s,*A8s,
         		                 *effC[iy], intCVec[iy],
         		                 true
         		                 );
    
    RooAbsReal* ang_wt = new ShapeSigAng(("PDF_sig_ang_wt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_wt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
					 *Fs,*As,*A4s,*A5s,*A7s,*A8s,
         		                 *effW[iy], intWVec[iy],
         		                 false 
         		                 );
    
    // Angular * mass component
    PdfSigAngMass* pdf_sig_ang_mass = nullptr;
    PdfSigAngMass* pdf_sig_ang_mass_penalty = nullptr;
    if (q2Bin < 5)  {
      pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
					    Form(sigpdfname.c_str(),years[iy]),
					    *ctK,*ctL,*phi,*mass,
					    *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
					    *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
					    *mFrac, *c_fm,
					    *ang_rt, *ang_wt,
					    *c_dcb_rt, *c_dcb_wt
					    );
    
      pdf_sig_ang_mass_penalty = new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
						    Form((sigpdfname+"_penalty").c_str(),years[iy]),
						    *ctK,*ctL,*phi,*mass,
						    *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
						    *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
						    *mFrac, *c_fm,
						    *penTerm,
						    *ang_rt, *ang_wt,
						    *c_dcb_rt, *c_dcb_wt
						    );
    }      		                                           
    else if (q2Bin==7){
        pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
                                              Form(sigpdfname.c_str(),years[iy]),
         		                      *ctK,*ctL,*phi,*mass,
                                              *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                              *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
         		                      *mFrac, *c_fm,
         		                      *ang_rt, *ang_wt,
         		                      *c_dcb_rt, *c_dcb_wt
         		                    );
    
        pdf_sig_ang_mass_penalty =  new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                                       Form((sigpdfname+"_penalty").c_str(),years[iy]),
      		                                       *ctK,*ctL,*phi,*mass,
                                                       *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                                       *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
      		                                       *mFrac, *c_fm,
                          		               *penTerm,
            		                               *ang_rt, *ang_wt,
      		                                       *c_dcb_rt, *c_dcb_wt
      		                                     );
    }      		                                           
    else {
      pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
					    Form(sigpdfname.c_str(),years[iy]),
					    *ctK,*ctL,*phi,*mass,
					    *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
					    *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
					    *mFrac, *c_fm,
					    *ang_rt, *ang_wt,
					    *c_dcb_rt, *c_dcb_wt
					    );
    
      pdf_sig_ang_mass_penalty = new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
						    Form((sigpdfname+"_penalty").c_str(),years[iy]),
						    *ctK,*ctL,*phi,*mass,
						    *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
						    *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
						    *mFrac, *c_fm,
						    *penTerm,
						    *ang_rt, *ang_wt,
						    *c_dcb_rt, *c_dcb_wt
						    );
    } 

    // Add constraint on mFrac to the pdf
    auto pdf_sig_ang_mass_mfc = new RooProdPdf(("PDF_sig_ang_mass_mfc_"+shortString+"_"+year).c_str(),
					       ("PDF_sig_ang_mass_mfc_"+year).c_str(),
					       *pdf_sig_ang_mass,
					       *c_fm);
    auto pdf_sig_ang_mass_penalty_mfc = new RooProdPdf(("PDF_sig_ang_mass_penalty_mfc_"+shortString+"_"+year).c_str(),
					       ("PDF_sig_ang_mass_penalty_mfc_"+year).c_str(),
					       *pdf_sig_ang_mass_penalty,
					       *c_fm);
    

    //// Background components ////
    // Read angular pdf for sidebands from external file 
    string filename_sb = Form("savesb_%i_b%i.root", years[iy], q2Bin);
    if (!localFiles) filename_sb = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/" + filename_sb;
    if (!(retrieveWorkspace(filename_sb, wsp_sb, "wsb"))) return;

    RooBernsteinSideband* bkg_ang_pdf = (RooBernsteinSideband*) wsp_sb[iy]->pdf(Form("BernSideBand_bin%i_%i", q2Bin, years[iy]));
    RooArgSet* bkg_ang_params = (RooArgSet*)bkg_ang_pdf->getParameters(observables);
    auto iter = bkg_ang_params->createIterator();
    RooRealVar* ivar =  (RooRealVar*)iter->Next();
    while (ivar) {
      ivar->setConstant(true);
      ivar = (RooRealVar*) iter->Next();
    }

    // create exponential mass pdf for background
    // in Jpsi bin, if fitOption > 0 -> bernstein polynbomial
    RooRealVar* slope = new RooRealVar(Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , -5., -10., 0.);
    RooAbsPdf* bkg_mass_pdf = 0;
    RooExponential* bkg_exp_pdf = 0;
    RooGenericPdf* pdfCut = 0;
    double pol_bmax =1.;
    RooRealVar* b0_bkg_mass = new RooRealVar(Form("b0_bkg_mass-%i",years[iy]) , Form("b0_bkg_mass-%i",years[iy])  ,  pol_bmax  );
    RooRealVar* b1_bkg_mass = new RooRealVar(Form("b1_bkg_mass-%i",years[iy]) , Form("b1_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
    RooRealVar* b2_bkg_mass = new RooRealVar(Form("b2_bkg_mass-%i",years[iy]) , Form("b2_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
    RooRealVar* b3_bkg_mass = new RooRealVar(Form("b3_bkg_mass-%i",years[iy]) , Form("b3_bkg_mass-%i",years[iy])  ,  0.0 );
    RooRealVar* b4_bkg_mass = new RooRealVar(Form("b4_bkg_mass-%i",years[iy]) , Form("b4_bkg_mass-%i",years[iy])  ,  0.1 , 0., pol_bmax);
    if (q2Bin==4 && fitOption>0){
      b0_bkg_mass->setConstant(true);
      b3_bkg_mass->setConstant(true);
      bkg_mass_pdf = new RooBernstein(Form("bkg_mass1_%i",years[iy]),Form("bkg_mass1_%i",years[iy]),  *mass, RooArgList(*b0_bkg_mass, *b1_bkg_mass, *b2_bkg_mass, *b3_bkg_mass,* b4_bkg_mass));
    } else { 
      bkg_exp_pdf = new RooExponential(Form("bkg_mass1_%i",years[iy]),  Form("bkg_mass1_%i",years[iy]) ,  *mass,   *slope  );
      if (q2Bin==3||q2Bin==5) {
       pdfCut = createPdfCuts(q2Bin,years[iy], mass);
       bkg_mass_pdf= new RooProdPdf(Form("bkg_mass_pdf_%i",years[iy]),Form("bkg_mass_pdf_%i",years[iy]),RooArgSet(*bkg_exp_pdf,*pdfCut));
      }else{
       bkg_mass_pdf=bkg_exp_pdf;
      } 
       
    } 

    RooProdPdf* bkg_pdf = new RooProdPdf(Form(bkgpdfname.c_str(),years[iy]),
					 Form(bkgpdfname.c_str(),years[iy]),
					 RooArgList(*bkg_ang_pdf,*bkg_mass_pdf));


    //// Finally build full pdf, including or not Z component ////
    RooAddPdf* full_pdf = 0;
    RooAddPdf* full_pdf_penalty = 0;

    // fraction of signal and bkg pdfs
    RooRealVar *fsig = new RooRealVar( ("fsig_"+shortString+"_"+year).c_str(), ("fsig_"+shortString+"_"+year).c_str(),0,1 );


    // special case of bin 6, include Z component in the signal model
    if(q2Bin==6 && fitOption>0){
      RooProdPdf* Z4430_pdf = 0;
      if(fitOption==3){
        std::cout << Form("Warning! Z(4430) Mass model from Signal MC Fits for Year=%i", years[iy]) << std::endl;
        RooConstVar* FracZ4430WT = new RooConstVar(Form("FracZ4430WT^{%i}",years[iy]),"FracZ4430WT", fraction);
        RooProduct* mass_Frac = new RooProduct( Form("mass_Frac_%i", years[iy]), Form("mass_Frac_%i", years[iy]), RooArgList(*FracZ4430WT, *mFrac));
        RooAddPdf* Mass_All = new RooAddPdf( Form("Mass_All_%i", years[iy]), Form("Mass_All_%i", years[iy]), RooArgList(*c_dcb_wt,*c_dcb_rt), *mass_Frac);
      	Z4430_pdf = new RooProdPdf( Form("Z4430_pdf_%i", years[iy]), Form("Z4430_pdf_%i", years[iy]), RooArgList(*Z4430_ang_pdf, *Mass_All, *c_fm) );
      }else{	
        std::cout << Form("Warning! Z(4430) Mass model from Z(4430) MC Fit for Year=%i", years[iy]) << std::endl;
      	Z4430_pdf = new RooProdPdf( Form("Z4430_pdf_%i",years[iy]), Form("Z4430_pdf_%i", years[iy]), RooArgList(*Z4430_ang_pdf, *Z4430_mass_pdf) );
      }	

      RooAddPdf* pdf_z_sig_ang_mass_mfc = new RooAddPdf(("PDF1_sig_ang_mass_"+shortString+"_"+year).c_str(),
					                ("PDF1_sig_ang_mass_"+year).c_str(),
   					                RooArgList(*Z4430_pdf,*pdf_sig_ang_mass_mfc),
 					                RooArgList( *AZ)
					               );
      
      RooAddPdf* pdf_z_sig_ang_mass_mfc_penalty = new RooAddPdf(("PDF1_sig_ang_mass_penalty _"+shortString+"_"+year).c_str(),
					                        ("PDF1_sig_ang_mass_penalty _"+year).c_str(),
   					                        RooArgList(*Z4430_pdf,*pdf_sig_ang_mass_penalty_mfc),
 					                        RooArgList( *AZ)
              					               );

      full_pdf = new RooAddPdf( ("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
 	     			("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
 	     			RooArgList(*pdf_z_sig_ang_mass_mfc, *bkg_pdf),
 	     			RooArgList(*fsig)
      	     		      );
      full_pdf_penalty = new RooAddPdf(("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
        			       ("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
        			       RooArgList(*pdf_z_sig_ang_mass_mfc_penalty, *bkg_pdf),
        			       RooArgList(*fsig)
        		              );

    }
    else {

      full_pdf = new RooAddPdf(("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
      	     		       ("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
      	     			RooArgList(*pdf_sig_ang_mass_mfc, *bkg_pdf),
      	     			RooArgList(*fsig)
      	     		      );
  
      full_pdf_penalty = new RooAddPdf(("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
      	     			       ("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
      	     				RooArgList(*pdf_sig_ang_mass_penalty_mfc, *bkg_pdf),
      	     				RooArgList(*fsig)
      	     		              );
    }
    PDF_sig_ang_mass_bkg.push_back(full_pdf);
    PDF_sig_ang_mass_bkg_penalty.push_back(full_pdf_penalty);

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if (multiSample) for (uint is = firstSample; is <= lastSample; is++) {
	if ( !data[iy][is] || data[iy][is]->IsZombie() ) {
	  cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
	  return;
	}
	map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
	simPdf        -> addPdf(*PDF_sig_ang_mass_bkg[iy],         ("data"+year+Form("_subs%d",is)).c_str());
	simPdf_penalty-> addPdf(*PDF_sig_ang_mass_bkg_penalty[iy], ("data"+year+Form("_subs%d",is)).c_str());
      }
    else {
      if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
	cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
	return;
      }
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
      simPdf        ->addPdf(*PDF_sig_ang_mass_bkg[iy],         ("data"+year+Form("_subs%d",firstSample)).c_str());
      simPdf_penalty->addPdf(*PDF_sig_ang_mass_bkg_penalty[iy], ("data"+year+Form("_subs%d",firstSample)).c_str());
    }
  
  }

  TFile* fout = 0;
  std::string addInfo = evalMistagSys ? "_mistagSyst_exp" : "";

  if (save>0) fout = new TFile(("simFitResults4d/simFitResult_data_fullAngularMass_Swave_" + all_years + stat + Form("_b%i%s_wp90_unbl%i%s.root", q2Bin, XGBstr.c_str(), unblind, addInfo.c_str())).c_str(),"RECREATE");
  RooWorkspace* wsp_out = 0;

  // save initial par values    
  RooArgSet *params      = (RooArgSet *)simPdf->getParameters(observables);
  RooArgSet* savedParams = (RooArgSet *)params->snapshot() ;
  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data", 
                            reco_vars,
                            Index(sample), 
                            Import(map)); 
  RooDataSet* combData = 0;

  // Results' containers
  double fitTime, imprTime, minTime;
  double co1, co4, co5;
  double boundDistFit, boundDist;
  bool boundCheck, convCheck, usedPenalty;

  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // TTree with the MINOS output
  vector<double> vResult  (pars.getSize());
  vector<double> vConfInterLow  (pars.getSize());
  vector<double> vConfInterHigh (pars.getSize());
  vector<double> vFreeResult  (pars.getSize());
  vector<double> vFreeConfInterLow  (pars.getSize());
  vector<double> vFreeConfInterHigh (pars.getSize());
  if (save>0) fout->cd();
  TTree* fitResultsTree = new TTree("fitResultsTree","fitResultsTree");
  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
    RooRealVar* par = (RooRealVar*)pars.at(iPar);
    fitResultsTree->Branch(Form("%s_low",par->GetName()),&vConfInterLow[iPar]);
    fitResultsTree->Branch(Form("%s_high",par->GetName()),&vConfInterHigh[iPar]);
    fitResultsTree->Branch(Form("%s_best",par->GetName()),&vResult[iPar]);
    fitResultsTree->Branch(Form("%s_free_low",par->GetName()),&vFreeConfInterLow[iPar]);
    fitResultsTree->Branch(Form("%s_free_high",par->GetName()),&vFreeConfInterHigh[iPar]);
    fitResultsTree->Branch(Form("%s_free_best",par->GetName()),&vFreeResult[iPar]);
  }
  fitResultsTree->Branch("fitTime",&fitTime);
  fitResultsTree->Branch("imprTime",&imprTime);
  fitResultsTree->Branch("minTime",&minTime);
  fitResultsTree->Branch("co1",&co1);
  fitResultsTree->Branch("co4",&co4);
  fitResultsTree->Branch("co5",&co5);
  fitResultsTree->Branch("boundDist",&boundDist);
  fitResultsTree->Branch("boundDistFit",&boundDistFit);
  fitResultsTree->Branch("boundCheck",&boundCheck);
  fitResultsTree->Branch("convCheck",&convCheck);
  fitResultsTree->Branch("usedPenalty",&usedPenalty);

  // Timer for fitting time
  TStopwatch subTime;

  // counters to monitor results' status
  int cnt[9];
  for (int iCnt=0; iCnt<9; ++iCnt) cnt[iCnt] = 0;

  Fitter* fitter = 0;
  vector<Fitter*> vFitter (0);

  for (uint is = firstSample; is <= lastSample; is++) {

    string samplename = Form("_subs%d", is);
    samplename =  "data%i" + samplename;
    string the_cut = Form(("sample==sample::"+samplename).c_str(),years[0]);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }
 
    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    if (nSample>0) cout<<"Fitting data subsample "<<is+1<<" with "<<combData->numEntries()<<" entries"<<endl;
    else cout<<"Fitting data sample with "<<combData->numEntries()<<" entries"<<endl;

     // set penalty term power parameter
    int combEntries = combData->numEntries();
    penTerm->setPower(power/combEntries);

    // to start the fit, parameters are restored to the center of the parameter space
    *params = *savedParams ;

    // run the fit
    fitter = new Fitter (Form("fitter%i",is),Form("fitter%i",is),pars,combData,simPdf,simPdf_penalty,boundary,bound_dist,penTerm,&c_vars);
    vFitter.push_back(fitter);

    fitter->setUnbl(unblind);

    bool runPostFitSteps = true;
    if (nSample==0 && (q2Bin==4 || q2Bin==6)) {
      // define if run improvAng and minosAng
      runPostFitSteps = false;
      fitter->setNCPU(8);
      // for bin 4, do not run the penalty even if needed
      if (q2Bin==4) fitter->runSimpleFit = true;
    }

    if (q2Bin==3 || q2Bin==5){
      fitter->DisableConstOpt();
      std::cout<<Form("DisableConstOpt for q2Bin=%d",q2Bin)<<std::endl;
    }  

    subTime.Start(true);
    int status = fitter->fit();
    subTime.Stop();

    fitTime=subTime.CpuTime();
    cout<<(fitter->runSimpleFit?"Fit time: ":"Fit+boundDist time: ")<<fitTime<<endl;

    co1=0;
    co4=0;
    co5=0;
    boundDist=0;
    boundDistFit=0;
    minTime=0;

    convCheck = false;
    boundCheck = false;

    if (status<2) {
      
      convCheck = true;
      boundCheck = boundary->getValV() == 0;

      if (unblind>3) fitter->result()->Print("v");

      if (fitter->runSimpleFit) boundDistFit = boundDist = -1;
      else boundDistFit = boundDist = fitter->boundDist;

      usedPenalty = fitter->usedPenalty;
	
      if (usedPenalty) {
	// save coefficient values
	co1 = fitter->coeff1;
	co4 = fitter->coeff4;
	co5 = fitter->coeff5;

	TStopwatch improvTime;
	improvTime.Start(true);
	if (runPostFitSteps) fitter->improveAng(1,300000);
	improvTime.Stop();
	imprTime = improvTime.CpuTime();
	cout<<"Improv time: "<<imprTime<<" s"<<endl;

	boundDist = fitter->boundDist;

      }

      // run MINOS error
      TStopwatch minosTime;
      minosTime.Start(true);
      
      if (runPostFitSteps) fitter->MinosAng();
      
      minosTime.Stop();
      minTime = minosTime.CpuTime();
      cout<<"MINOS errors computed in "<<minTime<<" s"<<endl;
      
      // cout<<"Error difference [custMINOS - fit], lower and higher:"<<endl;
      // for (int iPar = 0; iPar < pars.getSize(); ++iPar)
      // 	cout<<vResult[iPar]-vConfInterLow[iPar]+vFitErrLow[iPar]<<"   \t"
      // 	    <<vConfInterHigh[iPar]-vResult[iPar]-vFitErrHigh[iPar]<<endl;


      // save results in tree
      for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
      	vResult[iPar] = fitter->vResult[iPar];
      	vFreeResult[iPar] = fitter->vFreeFitResult[iPar];
      	if (runPostFitSteps) {
      	  vConfInterLow[iPar] = fitter->vConfInterLow[iPar] - vResult[iPar];
      	  vConfInterHigh[iPar] = fitter->vConfInterHigh[iPar] - vResult[iPar];
      	} else {
      	  vConfInterLow[iPar] = fitter->vFitErrLow[iPar];
      	  vConfInterHigh[iPar] = fitter->vFitErrHigh[iPar];
	}
	vFreeConfInterLow[iPar] = fitter->vFreeFitErrLow[iPar];
	vFreeConfInterHigh[iPar] = fitter->vFreeFitErrHigh[iPar];
      }
      fitResultsTree->Fill();

      if (plot && !multiSample && unblind>2) {

	string plotString = shortString + "_" + all_years + stat + XGBstr + Form("_unbl%i",unblind);
	addInfo = evalMistagSys ? "_mistagSyst_exp" : "";
	string plotname = "plotSimFit4d_d/simFitResult_data_fullAngularMass_Swave_" + plotString + addInfo + ".pdf";
	fitter->plotSimFitProjections(plotname.c_str(),{samplename,sigpdfname,bkgpdfname},years,true);

      }

      if (save>1 && (q2Bin!=4 || years.size()<3 || nSample>0) && unblind>1) {
	wsp_out = new RooWorkspace("wsp_out","wsp_out");
	if (unblind<4)
	  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
	    ((RooRealVar*)pars.at(iPar))->setVal(0);
	    ((RooRealVar*)pars.at(iPar))->setError(0);
	  }
	wsp_out->import(*combData);
	wsp_out->import(*simPdf);
      }

    }

    // fill fit-status-dependent counters
    ++cnt[8];
    int iCnt = 0;
    if (!convCheck) iCnt += 4;
    if (!boundCheck) iCnt += 2;
    if (usedPenalty) iCnt += 1;
    ++cnt[iCnt];

    // print fit status and time
    if (!boundCheck)
      if (convCheck) cout<<"Converged in unphysical region";
      else cout<<"Not converged";
    else
      if (convCheck)
	if (usedPenalty) cout<<"Converged with penalty term with coeff: "<<fitter->coeff1<<" "<<fitter->coeff4<<" "<<fitter->coeff5;
	else cout<<"Converged without penalty";
      else cout<<"This should never be printed";
    cout<<" ("<<fitTime<<"s)"<<endl;

  }  


  if (multiSample) {
    cout<<"Fitted subsamples: "<<cnt[8]<<" of which good: "<<cnt[0]+cnt[1]<<" ("<<cnt[1]<<" with the use of the penalty term)"<<endl;
    cout<<"Bad fits: "<<cnt[3]<<" converging outside physical region, "<<cnt[5]+cnt[7]<<" not converged ("<<cnt[5]<<" in ph region)"<<endl;
  }

  if (save>0 && unblind>1) {
    fout->cd();
    if (unblind>3) fitResultsTree->Write();
    if (wsp_out) wsp_out->Write();
    fout->Close();
  }

}


void simfit_data_fullAngularMass_SwaveBin1(int q2Bin, int parity, bool multiSample, uint nSample, uint q2stat, int fitOption, int XGBv, int unblind, bool localFiles, bool plot, int save, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_data_fullAngularMass_SwaveBin(q2Bin, parity, multiSample, nSample, q2stat, fitOption, XGBv, unblind, localFiles, plot, save, years);
  else
    simfit_data_fullAngularMass_SwaveBin(q2Bin, parity, multiSample, nSample, q2stat, fitOption, XGBv, unblind, localFiles, plot, save, years);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively
  // unblind options: [0] no fit in signal q2 bins
  //                  [1] fit reports only time, convergence and penalty
  //                  [2] adds root file with RooWorkspace and blinded POI
  //                  [3] adds fit projections plot
  //                  [4] adds RooFitResult and MINOS errors to log, unblind POI and add TTree in root file

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  bool multiSample = false;
  uint nSample = 0;
  uint q2stat = 0;
  if ( argc > 3 && atoi(argv[3]) > 0 ) multiSample = true;
  if ( argc > 4 ) nSample = atoi(argv[4]);
  if ( argc > 5 ) q2stat = atoi(argv[5]);

  if (nSample==0) multiSample = false;

  int fitOption=0;
  // possible configurations and meanings of alternative configurations:
  // bin 4: fitOpt == 1 -> bkg mass pdf bernstein polynomial instead of expo
  // bin 6: fitOpt == 1 -> include Z component, all pars fixed to MC
  //        fitOpt == 2 -> include Z component, some pars are floating 
  //        fitOpt == 3 -> include Z component, Z mass model from the psi2S MC
  if ( argc > 6 ) fitOption = atoi(argv[6]);

  int XGBv = 0;
  if ( argc > 7 ) XGBv = atoi(argv[7]);

  int unblind=0;
  if ( argc > 8 ) unblind = atoi(argv[8]);

  bool localFiles = false;
  if ( argc > 9 && atoi(argv[9]) > 0 ) localFiles = true;

  bool plot = true;
  int save = true;

  if ( argc > 10 && atoi(argv[10]) == 0 ) plot = false;
  if ( argc > 11 ) save = atoi(argv[11]);

  std::vector<int> years;
  if ( argc > 12 && atoi(argv[12]) != 0 ) years.push_back(atoi(argv[12]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 13 && atoi(argv[13]) != 0 ) years.push_back(atoi(argv[13]));
  if ( argc > 14 && atoi(argv[14]) != 0 ) years.push_back(atoi(argv[14]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "q2stat      " << q2stat       << endl;
  cout <<  "fitOption   " << fitOption    << endl;
  cout <<  "XGB version " << XGBv         << endl;
  cout <<  "local files " << localFiles   << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;

  cout <<  "UNBLINDIG STEP " << unblind << endl;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2stat  <  0 || q2stat  >= nBins ) return 1;

  if ( XGBv < 0 ) return 1;

  // Protectrion against accidental unblinding
  if ( unblind < 1 && q2Bin != 4 && q2Bin != 6 ) {
    cout<<"The analysis is blind!"<<endl;
    return 1;
  }

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_data_fullAngularMass_SwaveBin1(q2Bin, parity, multiSample, nSample, q2stat, fitOption, XGBv, unblind, localFiles, plot, save, years);
  else
    simfit_data_fullAngularMass_SwaveBin1(q2Bin, parity, multiSample, nSample, q2stat, fitOption, XGBv, unblind, localFiles, plot, save, years);

  return 0;

}
