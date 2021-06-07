#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>

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
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooConstVar.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"

#include "PdfSigMass.h"
#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
std::map<int,float> scale_to_data;

TCanvas* c [4*nBins];

void simfit_recoMC_fullMassBin(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  string all_years = "";
  string year = ""; 
  string isample = ""; 
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;

  std::vector<RooWorkspace*> wsp, wsp_mcmass;
  std::vector<std::vector<RooDataSet*>> data;
  std::vector<RooAbsPdf*> PDF_sig_mass(0);
  std::vector<RooGaussian*> c_deltaPeaks, c_fm;
  RooArgSet c_vars_rt, c_pdfs_rt;
  RooArgSet c_vars_wt, c_pdfs_wt;
  RooArgSet c_vars; 
  RooWorkspace * ws_pars = new RooWorkspace("ws_pars");

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6);
  RooRealVar* rand = new RooRealVar("rand","rand", 0,1);
  RooArgSet reco_vars ( *mass, *rand);
  RooArgSet observables (*mass);


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

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/recoMCDataset_b%i_%i.root", years[iy], q2Bin, years[iy]); 

    // import data (or MC as data proxy)
    if (!retrieveWorkspace( filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity )))  return;

    // create roodataset (in case data-like option is selected, only import the correct % of data)
    data.push_back( createDataset( nSample,  firstSample,  lastSample, wsp[iy],  
                                   q2Bin,  parity,  years[iy], 
                                   reco_vars, observables,  shortString  )); 

    // Mass Component
    // import mass PDF from fits to the MC
    string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",years[iy]);
    if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w"))  return;

    RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",0.13, 0, 1);
    //   mFrac->setConstant();
 
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    RooRealVar* mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
    RooRealVar* sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
    RooRealVar* alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.,  100.);

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
      if (q2Bin < 7) 
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
      else
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
       
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
    RooRealVar* n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0., 100.);
    RooRealVar* n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0., 100.);

    RooAbsPdf* dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );

    /// create constrained PDF for RT mass
    RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
    constr_wt_list.add(*dcb_wt);
    RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
    c_vars.add(c_vars_wt);

    if (q2Bin < 5)  
        PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                               ("PDF_sig_mass_"+year).c_str(),
                                               *mass,
                                               *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                                               *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
      		                               *mFrac,
      		                               *c_dcb_rt,
      		                               *c_dcb_wt  
      		                              ));
    else if (q2Bin==7)
        PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                               ("PDF_sig_mass_"+year).c_str(),
                                               *mass,
                                               *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                               *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
      		                               *mFrac,
      		                               *c_dcb_rt ,
      		                               *c_dcb_wt
      		                               ));
    else  
        PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                               ("PDF_sig_mass_"+year).c_str(),
                                               *mass,
                                               *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                                               *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
      		                               *mFrac,
      		                               *c_dcb_rt ,
      		                               *c_dcb_wt
      		                               ));

    
    /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
    double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
    double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
    double fraction = nwt_mc / (nrt_mc + nwt_mc);
    c_fm.push_back(new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                    RooConst(fraction) , 
                                    RooConst(fM_sigmas[years[iy]][q2Bin])
                                    ) );
    cout << fraction << "   " << fM_sigmas[years[iy]][q2Bin] << endl;                                    
    c_vars.add(*mFrac);       

    RooProdPdf* final_PDF = new RooProdPdf(("final_PDF_"+year).c_str(), 
                                           ("final_PDF_"+year).c_str(),
                                           RooArgList(*PDF_sig_mass[iy],*c_fm[iy]));
   
    
    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if (multiSample) for (uint is = firstSample; is <= lastSample; is++) {
	if ( !data[iy][is] || data[iy][is]->IsZombie() ) {
	  cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
	  return;
	}
	map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
	simPdf->addPdf(*final_PDF, ("data"+year+Form("_subs%d",is)).c_str());
    }
    else {
      if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
	cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
	return;
      }
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
      simPdf->addPdf(*final_PDF, ("data"+year+Form("_subs%d",firstSample)).c_str());
    }
  }

   
  TFile* fout = new TFile(("simFitResults4d/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"UPDATE");

  // save initial par values into a workspace 
  ws_pars->import(*simPdf);
  RooArgSet *params = (RooArgSet *)simPdf->getParameters(*mass);
  // The kTRUE flag imports the values of the objects in (*params) into the workspace
  // If not set, the present values of the workspace parameters objects are stored
  ws_pars->saveSnapshot("initial_pars", *params, kTRUE);
 
  
  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data", 
                            reco_vars,
                            Index(sample), 
                            Import(map)); 
  RooDataSet* combData = 0;
  RooAbsReal* nll = 0;
 
  for (uint is = firstSample; is <= lastSample; is++) {

    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }
    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    if (nSample>0) cout<<"Fitting subsample "<<is+1<<" with "<<combData->numEntries()<<" entries"<<endl;
    else cout<<"Fitting full MC sample with "<<combData->numEntries()<<" entries"<<endl;

    ws_pars->loadSnapshot("initial_pars");

    TStopwatch subTime;
    nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                            RooFit::Extended(kFALSE),
                                            RooFit::Constrain(c_vars),
                                            RooFit::NumCPU(1)
                                            );
         
    RooMinimizer m(*nll) ;   
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    m.setPrintLevel(-1);
    m.setPrintEvalErrors(-1);
    m.setMinimizerType("Minuit2");
//     m.setVerbose(kTRUE); 

    subTime.Start(true);
    m.setStrategy(0);
    m.migrad() ;
    m.hesse() ;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;
    subTime.Stop();
    cout << "fitting done  " << subTime.CpuTime() << endl;
    
    RooFitResult* fitResult = m.save(("result_" + shortString + Form("subs%d",is)).c_str()) ; 
    fitResult->Print("v");

    // Save fit results in file
    if (save) {
      fout->cd();
      fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
    }
  }  

  fout->Close();
  if (!plot || multiSample) return;

  string plotString = shortString + "_" + all_years;
  if (nSample>0) plotString = plotString + Form("_s%i",nSample);

  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections 
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),3000,1400);
  c[confIndex]->Divide(years.size());
  cout<<"plotting 4d canvas: "<< years.size() << endl;

  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    std::vector<RooPlot*> frames;
    frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
 
    for (unsigned int fr = 0; fr < frames.size(); fr++){
        combData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), 
                         Cut(("sample==sample::data"+year+Form("_subs%d",firstSample)).c_str()), 
                         Name(("plData"+year).c_str()));
        simPdf  ->plotOn(frames[fr], LineWidth(1),
                         Slice(sample, ("data"+year+Form("_subs%i",firstSample)).c_str()), 
                         ProjWData(RooArgSet(sample), *combData),  
                         Name(("plPDF"+year).c_str()));
        if (fr == 0) { 
          leg->AddEntry(frames[fr]->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+year ).c_str()),("Mass component of the pdf "+year).c_str(),"l");
        }
        c[confIndex]->cd(iy+1);
        gPad->SetLeftMargin(0.19); 
        frames[fr]->Draw();
        leg->Draw("same");
        break;
    }
  }
  c[confIndex]->SaveAs( ("plotSimFit4d_d/simFitResult_recoMC_fullMass_" + plotString +  ".pdf").c_str() );

}



void simfit_recoMC_fullMassBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, years);
  else
    simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, years);
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

  bool plot = true;
  bool save = true;

  if ( argc > 5 && atoi(argv[5]) == 0 ) plot = false;
  if ( argc > 6 && atoi(argv[6]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc > 7 && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 8 && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  if ( argc > 9 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
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
      simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, years);
  else
    simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, years);

  return 0;

}
