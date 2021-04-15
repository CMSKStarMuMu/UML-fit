#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>

// #include <RooRealVar.h>
#include <RooAbsPdf.h>
// #include <RooWorkspace.h>
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
// #include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include "RooExponential.h"

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
std::map<int,float> scale_to_data;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [4*nBins];

double power = 1.0;

void simfit_data_fullAngularMassBin(int q2Bin, int parity, bool multiSample, uint nSample, bool localFiles, bool plot, bool save, std::vector<int> years)
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
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;
  
  std::vector<TFile*> fin_eff;
  std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_sb;
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

  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6);
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

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("recoDATADataset_b%i_%i.root", q2Bin, years[iy]);
    if (!localFiles) filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/", years[iy]) + filename_data;

    // import data (or MC as data proxy)
    retrieveWorkspace( filename_data, wsp, Form("ws_b%ip0", q2Bin ));

    // import KDE efficiency histograms and partial integral histograms
    string filename = Form((parity==0 ? "KDEeff_b%i_ev_%i.root" : "KDEeff_b%i_od_%i.root"),q2Bin,years[iy]);
    if (!localFiles) filename = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/",years[iy]) + filename;
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
    data.push_back( createDatasetInData( wsp[iy],  q2Bin,  observables,  shortString  )); 
    // data.push_back( createDataset( nSample,  firstSample,  lastSample, wsp[iy],  
    //                                q2Bin,  parity,  years[iy], 
    //                                reco_vars, observables,  shortString  )); 

    // Mass Component
    // import mass PDF from fits to the MC
    string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_Jpsi_newbdt.root",years[iy]);
    // string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",years[iy]);
    if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w"))  return;

    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    RooRealVar* mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
    RooRealVar* sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
    RooRealVar* alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* alpha_rt2     = new RooRealVar (Form("#alpha_{RT2}^{%i}",years[iy] ), "alphart2"    , wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal() ,    -10,   10 );
    RooRealVar* n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.01,  100.);
    RooRealVar* n_rt2         = new RooRealVar (Form("n_{RT2}^{%i}",years[iy])      , "nrt2"        , wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin))->getVal()      ,      0.01,  100.);

    RooAbsPdf* dcb_rt;
    RooRealVar* sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,   0 , 0,   0.12, "GeV");
    RooRealVar* f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   0 , 0.,  1.);
    if (q2Bin >= 4){
      sigma_rt2-> setVal(wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal() );
      f1rt     -> setVal(wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal() );
      dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
    } 
    else{
        alpha_rt2->setRange(0,10);
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt  );
    }
    
    /// create constrained PDF for RT mass
    RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
    constr_rt_list.add(*dcb_rt);
    RooProdPdf * c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list );
    c_vars.add(c_vars_rt);
   
    /// create WT component
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));

    RooRealVar* mean_wt     = new RooRealVar (Form("mean_{WT}^{%i}",years[iy])      , "masswt"     ,  wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal()    ,      5,    6, "GeV");
    RooRealVar* sigma_wt    = new RooRealVar (Form("#sigma_{WT1}^{%i}",years[iy])   , "sigmawt"    ,  wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal() ,      0,    1, "GeV");
    RooRealVar* alpha_wt1   = new RooRealVar (Form("#alpha_{WT1}^{%i}",years[iy] )  , "alphawt1"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* alpha_wt2   = new RooRealVar (Form("#alpha_{WT2}^{%i}",years[iy] )  , "alphawt2"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);
    RooRealVar* n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);

    RooAbsPdf* dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );

    /// create constrained PDF for WT mass
    RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
    constr_wt_list.add(*dcb_wt);
    RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
    c_vars.add(c_vars_wt);

    cout << "deltap built --> constraint not added yet (to be done)" << endl;
    //// creating constraints for the difference between the two peaks
//     RooFormulaVar* deltaPeaks = new RooFormulaVar(Form("deltaPeaks^{%i}", years[iy]), "@0 - @1", RooArgList(*mean_rt, *mean_wt))  ;
//     c_deltaPeaks.push_back(     new RooGaussian(Form("c_deltaPeaks^{%i}", years[iy]), "c_deltaPeaks", *deltaPeaks, 
//                                                 RooConst( deltaPeaks->getVal() ), 
//                                                 RooConst( 0.0005 ) 
//                                                ) );
//     c_vars.add(*deltaPeaks);       c_pdfs.add(*c_deltaPeaks[iy]);

    RooRealVar* mFrac = new RooRealVar(Form("f_{M}^{%i}",years[iy]),"mistag fraction",1, 0.5, 1.5);
    /// create constraint on mFrac (mFrac = 1, constraint derived from stat scaling)
    double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
    double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
    double fraction = nwt_mc / (nrt_mc + nwt_mc);
    double frac_sigma = fM_sigmas[years[iy]][q2Bin]/fraction;
    RooGaussian* c_fm = new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                        RooConst(1.) , 
                                        RooConst(frac_sigma)
                                        );
//     cout << fraction << " +/- " << fM_sigmas[years[iy]][q2Bin] << " ---> << "1 +/- " << frac_sigma << endl;                                    
    c_vars.add(*mFrac); 

    // Angular Component
    RooAbsReal* ang_rt = new ShapeSigAng( ("PDF_sig_ang_rt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_rt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
         		                 *effC[iy], intCVec[iy],
         		                 true
         		                 );
    
    RooAbsReal* ang_wt = new ShapeSigAng( ("PDF_sig_ang_wt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_wt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
         		                 *effW[iy], intWVec[iy],
         		                 false 
         		                 );
    
    PdfSigAngMass* pdf_sig_ang_mass = nullptr;
    PdfSigAngMass* pdf_sig_ang_mass_penalty = nullptr;
    if (q2Bin < 5)  {
        pdf_sig_ang_mass = new PdfSigAngMass( ("PDF_sig_ang_mass_"+shortString+"_"+year).c_str(),
					      ("PDF_sig_ang_mass_"+year).c_str(),
					      *ctK,*ctL,*phi,*mass,
					      *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
					      *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
					      *mFrac, *c_fm,
					      *ang_rt, *ang_wt,
					      *c_dcb_rt, *c_dcb_wt
					      );
    
        pdf_sig_ang_mass_penalty = new PdfSigAngMass( ( "PDF_sig_ang_mass_penalty_"+shortString+"_"+year).c_str(),
						      ( "PDF_sig_ang_mass_penalty_"+year).c_str(),
						      *ctK,*ctL,*phi,*mass,
						      *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
						      *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
						      *mFrac, *c_fm,
						      *penTerm,
						      *ang_rt, *ang_wt,
						      *c_dcb_rt, *c_dcb_wt
						      );
    }      		                                           
    else {
        pdf_sig_ang_mass = new PdfSigAngMass( ("PDF_sig_ang_mass_"+shortString+"_"+year).c_str(),
					      ("PDF_sig_ang_mass_"+year).c_str(),
					      *ctK,*ctL,*phi,*mass,
					      *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
					      *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
					      *mFrac, *c_fm,
					      *ang_rt, *ang_wt,
					      *c_dcb_rt, *c_dcb_wt
					      );
    
        pdf_sig_ang_mass_penalty = new PdfSigAngMass( ("PDF_sig_ang_mass_penalty_"+shortString+"_"+year).c_str(),
						      ("PDF_sig_ang_mass_penalty_"+year).c_str(),
						      *ctK,*ctL,*phi,*mass,
						      *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
						      *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
						      *mFrac, *c_fm,
						      *penTerm,
						      *ang_rt, *ang_wt,
						      *c_dcb_rt, *c_dcb_wt
						      );
    } 
    

    // Read angular pdf for sidebands from external file 
    string filename_sb = Form("savesb_%i_b%i.root", years[iy], q2Bin );
    // if (!localFiles) filename_sb = "/eos/cms/store/user/fiorendi/p5prime/sidebands/" + filename_sb;
    retrieveWorkspace( filename_sb, wsp_sb, "wsb");

    RooBernsteinSideband* bkg_ang_pdf = (RooBernsteinSideband*) wsp_sb[iy]->pdf(Form("bkg_mass_sb_bin%i_%i", q2Bin, years[iy]));
    // RooBernsteinSideband* bkg_ang_pdf = (RooBernsteinSideband*) wsp_sb[iy]->pdf(Form("BernSideBand_bin%i_%i", q2Bin, years[iy]));
  RooArgSet* bkg_ang_params = (RooArgSet*)bkg_ang_pdf->getParameters(observables);
  auto iter = bkg_ang_params->createIterator();
  RooRealVar* ivar =  (RooRealVar*)iter->Next();
  while (ivar) {
    ivar->setConstant(true);
    ivar = (RooRealVar*) iter->Next();
  }

    // read mass pdf for background
    RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , -1., -10., 0.);
    // RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);
    RooExponential* bkg_exp = new RooExponential(Form("bkg_exp_%i",years[iy]),  Form("bkg_exp_%i",years[iy]) ,  *slope,   *mass  );

    // create 4D pdf  for background and import to workspace
    RooProdPdf* bkg_pdf = new RooProdPdf(Form("bkg_pdf_%i",years[iy]), Form("bkg_pdf_%i",years[iy]), RooArgList(*bkg_ang_pdf,*bkg_exp)); 


    // sum signal and bkg pdf 
    RooRealVar *fsig = new RooRealVar( ("fsig_"+shortString+"_"+year).c_str(), ("fsig_"+shortString+"_"+year).c_str(),0,1 );
    RooAddPdf* full_pdf = new RooAddPdf( ("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
                                          RooArgList(*pdf_sig_ang_mass, *bkg_pdf),
                                          RooArgList(*fsig)
                                       );
    PDF_sig_ang_mass_bkg.push_back(full_pdf);
    
    RooAddPdf* full_pdf_penalty = new RooAddPdf( ("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
                                                 ("PDF_sig_ang_fullAngularMass_bkg_penalty_"+shortString+"_"+year).c_str(),
                                                  RooArgList(*pdf_sig_ang_mass_penalty, *bkg_pdf),
                                                  RooArgList(*fsig)
                                       );
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

  if (nSample>0)   stat = stat + Form("-%i",firstSample);
  if (multiSample) stat = stat + Form("-%i",lastSample);
  TFile* fout = 0;
  if (save) fout = new TFile(("simFitResults4d/simFitResult_data_fullAngularMass" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"RECREATE");
  
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
  if (save) fout->cd();
  TTree* fitResultsTree = new TTree("fitResultsTree","fitResultsTree");
  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
    RooRealVar* par = (RooRealVar*)pars.at(iPar);
    fitResultsTree->Branch(Form("%s_low",par->GetName()),&vConfInterLow[iPar]);
    fitResultsTree->Branch(Form("%s_high",par->GetName()),&vConfInterHigh[iPar]);
    fitResultsTree->Branch(Form("%s_best",par->GetName()),&vResult[iPar]);
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

    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
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

    subTime.Start(true);
    int status = fitter->fit();
    subTime.Stop();

    fitTime=subTime.CpuTime();
    cout<<"Fit+boundDist time: "<<fitTime<<endl;

    co1=0;
    co4=0;
    co5=0;
    boundDist=0;
    boundDistFit=0;
    minTime=0;

    convCheck = false;
    boundCheck = false;

    if (status==0) {
      
      convCheck = true;
      boundCheck = boundary->getValV() == 0;

      fitter->result()->Print("v");

      boundDistFit = boundDist = fitter->boundDist;
      usedPenalty = fitter->usedPenalty;
	
      if (fitter->usedPenalty) {
	// save coefficient values
	co1 = fitter->coeff1;
	co4 = fitter->coeff4;
	co5 = fitter->coeff5;

	TStopwatch improvTime;
	improvTime.Start(true);
	// fitter->improveAng();
	improvTime.Stop();
	imprTime = improvTime.CpuTime();
	cout<<"Improv time: "<<imprTime<<" s"<<endl;

	boundDist = fitter->boundDist;

      }

      if (nSample>0) {
	// run MINOS error
	TStopwatch minosTime;
	minosTime.Start(true);

	// fitter->MinosAng();

	minosTime.Stop();
	minTime = minosTime.CpuTime();
	cout<<"MINOS errors computed in "<<minTime<<" s"<<endl;

	// cout<<"Error difference [custMINOS - fit], lower and higher:"<<endl;
	// for (int iPar = 0; iPar < pars.getSize(); ++iPar)
	// 	cout<<vResult[iPar]-vConfInterLow[iPar]+vFitErrLow[iPar]<<"   \t"
	// 	    <<vConfInterHigh[iPar]-vResult[iPar]-vFitErrHigh[iPar]<<endl;

      }

      // save results in tree
      for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
	vResult[iPar] = fitter->vResult[iPar];
	if (nSample>0) {
	  vConfInterLow[iPar] = fitter->vConfInterLow[iPar];
	  vConfInterHigh[iPar] = fitter->vConfInterHigh[iPar];
	} else {
	  vConfInterLow[iPar] = fitter->vFitErrLow[iPar];
	  vConfInterHigh[iPar] = fitter->vFitErrHigh[iPar];
	}
      }
      fitResultsTree->Fill();

    }

    // fill fit-status-dependent counters
    ++cnt[8];
    int iCnt = 0;
    if (!convCheck) iCnt += 4;
    if (!boundCheck) iCnt += 2;
    if (fitter->usedPenalty) iCnt += 1;
    ++cnt[iCnt];

    // print fit status and time
    if (!boundCheck)
      if (convCheck) cout<<"Converged in unphysical region";
      else cout<<"Not converged";
    else
      if (convCheck)
	if (fitter->usedPenalty) cout<<"Converged with penalty term with coeff: "<<fitter->coeff1<<" "<<fitter->coeff4<<" "<<fitter->coeff5;
	else cout<<"Converged without penalty";
      else cout<<"This should never be printed";
    cout<<" ("<<fitTime<<"s)"<<endl;

  }  


  if (multiSample) {
    cout<<"Fitted subsamples: "<<cnt[8]<<" of which good: "<<cnt[0]+cnt[1]<<" ("<<cnt[1]<<" with the use of the penalty term)"<<endl;
    cout<<"Bad fits: "<<cnt[3]<<" converging outside physical region, "<<cnt[5]+cnt[7]<<" not converged ("<<cnt[5]<<" in ph region)"<<endl;
  }

  if (save) {
    fout->cd();
    fitResultsTree->Write();
    fout->Close();
  }

  if (!plot || multiSample) return;

  // For plotting the effective penalty term is used
  Penalty* penTerm_eff = new Penalty(*penTerm,"penTerm_eff");
  penTerm_eff->setPower(power);
  RooFormulaVar* penLog = new RooFormulaVar("penLog","penLog","-1.0 * log(penTerm_eff)",RooArgList(*penTerm_eff));

  double xZoom = 200.0;
  if (nSample>0) xZoom = 2.0;

  cnll  = new TCanvas (("cnll_"+shortString).c_str(),("cnll_"+shortString).c_str(),1800,1800);
  cZoom = new TCanvas (("cZoom_"+shortString).c_str(),("cZoom_"+shortString).c_str(),1800,1800);
  cPen = new TCanvas (("cPen_"+shortString).c_str(),("cPen_"+shortString).c_str(),1800,1800);
  cnll->Divide(3,3);
  cZoom->Divide(3,3);
  cPen->Divide(3,3);

  RooPlot* frame [8];
  RooPlot* fZoom [8];
  RooPlot* fPenTerm [8];

  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)pars.at(iPar);

    frame[iPar] = par->frame(Name(Form("f1%s",par->GetName())),Title(Form("-log(L) scan vs %s",par->GetTitle()))) ;
    fZoom[iPar] = par->frame(Name(Form("f2%s",par->GetName())),Title(Form("zoom on -log(L) scan vs %s",par->GetTitle())),
			     Range(TMath::Max(par->getMin(),par->getValV()+xZoom*par->getErrorLo()),
				   TMath::Min(par->getMax(),par->getValV()+xZoom*par->getErrorHi()) )) ;

    fitter->nll->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;
    fitter->nll->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;

    if (iPar>0) {

      double hMax = frame[iPar]->GetMaximum();

      boundary->plotOn(frame[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      boundary->plotOn(fZoom[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));

      if (fitter->usedPenalty) {

	fitter->nll_penalty->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll_penalty->getVal()+10),LineColor(kBlue),LineWidth(2));
	fitter->nll_penalty->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll_penalty->getVal()+10),LineColor(kBlue),LineWidth(2));

	penLog->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));
	penLog->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));

      }  

      frame[iPar]->SetMaximum(hMax);

      fPenTerm[iPar] = par->frame(Name(Form("f3%s",par->GetName())),Title(Form("Penalty term vs %s",par->GetTitle()))) ;
      penTerm_eff->plotOn(fPenTerm[iPar],LineColor(4),LineWidth(2)) ;
      double hMaxP = fPenTerm[iPar]->GetMaximum();
      boundary->plotOn(fPenTerm[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMaxP,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      fPenTerm[iPar]->SetMaximum(hMaxP);
      cPen->cd(iPar+1);
      fPenTerm[iPar]->Draw();

    }

    fZoom[iPar]->SetMaximum(0.5*xZoom*xZoom);

    cnll->cd(iPar+1);
    frame[iPar]->Draw();

    cZoom->cd(iPar+1);
    fZoom[iPar]->Draw();


  }

  string plotString = shortString + "_" + all_years;
  if (nSample>0) plotString = plotString + Form("_s%i",nSample);

  cnll->Update();
  cnll->SaveAs( ("plotSimFit4d_d/recoNLL_scan_" + plotString + ".pdf").c_str() );

  cZoom->Update();
  cZoom->SaveAs( ("plotSimFit4d_d/recoNLL_scan_" + plotString + "_zoom.pdf").c_str() );

  cPen->Update();
  cPen->SaveAs( ("plotSimFit4d_d/recoPenTerm_" + plotString + ".pdf").c_str() );

   
  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections 
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),3000,1400);
  c[confIndex]->Divide(4, years.size());
  
  cout<<"plotting 4d canvas"<<endl;
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    std::vector<RooPlot*> frames;
    frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
    frames.push_back( prepareFrame( ctK ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( ctL ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( phi ->frame(Title((longString+year).c_str())) ));
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

    cout<<"canvas ready"<<endl;
    for (unsigned int fr = 0; fr < frames.size(); fr++){
        cout<<"fr " << fr<<endl;
        combData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), Cut(("sample==sample::data"+year+Form("_subs%d",firstSample)).c_str()), Name(("plData"+year).c_str()));
        
        simPdf->plotOn(frames[fr], Slice(sample, ("data"+year+Form("_subs%d",firstSample)).c_str()), 
                                     ProjWData(RooArgSet(sample), *combData), 
                                     LineWidth(1), 
                                     Name(("plPDF"+year).c_str()), 
                                     NumCPU(4));
        if (fr == 0) { 
          leg->AddEntry(frames[fr]->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
        }
        c[confIndex]->cd(iy*4+fr+1);
        gPad->SetLeftMargin(0.19); 
        frames[fr]->Draw();
        leg->Draw("same");
    }
  }
  c[confIndex]->SaveAs( ("plotSimFit4d_d/simFitResult_data_fullAngularMass_" + plotString +  ".pdf").c_str() );

}



void simfit_data_fullAngularMassBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool localFiles, bool plot, bool save, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_data_fullAngularMassBin(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years);
  else
    simfit_data_fullAngularMassBin(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  bool multiSample = false;
  uint nSample = 0;
  if ( argc > 3 && atoi(argv[3]) > 0 ) multiSample = true;
  if ( argc > 4 ) nSample = atoi(argv[4]);

  if (nSample==0) multiSample = false;

  bool localFiles = false;
  if ( argc > 5 && atoi(argv[5]) > 0 ) localFiles = true;

  bool plot = true;
  bool save = true;

  if ( argc > 6 && atoi(argv[6]) == 0 ) plot = false;
  if ( argc > 7 && atoi(argv[7]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc > 8 && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 9  && atoi(argv[9])  != 0 ) years.push_back(atoi(argv[9]));
  if ( argc > 10 && atoi(argv[10]) != 0 ) years.push_back(atoi(argv[10]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "local files " << localFiles   << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;


  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
  scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_data_fullAngularMassBin1(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years);
  else
    simfit_data_fullAngularMassBin1(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years);

  return 0;

}
