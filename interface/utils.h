#include <TSystem.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooGenericPdf.h>
#include <map>

using namespace std;
using namespace RooFit;

std::map<int, float> scale_to_data = {
  {2016, 0.006 *2 /1.45 },
  {2017, 0.005 *2 /1.40 },
  {2018, 0.007 *2 /1.30 }
// *2 since we are using only odd/even events, second factor is "data-driven":
// https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
};

// redo jpsi and psi2s
std::map<int,std::vector<float>> frt_sigmas = {
  {2016, {0.021, 0.015, 0.016, 0.011, 0.009, 0.009, 0.008, 0.011}},
  {2017, {0.018, 0.013, 0.014, 0.010, 0.008, 0.008, 0.007, 0.011}},
  {2018, {0.013, 0.010, 0.011, 0.007, 0.006, 0.006, 0.006, 0.007}},
};

std::map<int,std::vector<float>> fM_sigmas = {
  {2016, {0.023, 0.015, 0.017, 0.013, 0.0005 , 0.010, 0.0018, 0.013}},
  {2017, {0.018, 0.014, 0.015, 0.010, 0.0004 , 0.008, 0.0016, 0.011}},
  {2018, {0.015, 0.010, 0.011, 0.008, 0.00027, 0.006, 0.0011, 0.008}},
};

std::map<int,std::vector<float>> nbkg_years = {
  {2016, {162, 535, 462,  810, 0.005, 1342, 0.006, 467}},
  {2017, {185, 496, 441,  711, 0.004, 1363, 0.005, 379}},
  {2018, {288, 842, 734, 1270, 0.002, 2954, 0.003, 779}},
};

std::map<int,std::vector<float>> nsig_years = {
  {2016, {205, 454, 391,  689, 0.005, 1174, 0.006,  704}},
  {2017, {307, 581, 495, 1013, 0.004, 1524, 0.005,  835}},
  {2018, {500, 981, 821, 1608, 0.002, 3018, 0.003, 1836}},
};

std::map<int,std::vector<float>> constr_width_fac = {
  {2016, {2.10, 4.02, 1.32, 3.08, 1.47, 2.88, 1.22, 3.37, 1.00, 1.00, 2.33, 1.83, 0.00, 0.00, 2.16, 2.18}},
  {2017, {1.76, 3.61, 1.18, 2.32, 1.36, 2.90, 1.06, 2.50, 1.00, 1.00, 2.19, 1.76, 0.00, 0.00, 1.99, 2.01}},
  {2018, {1.87, 4.17, 1.16, 2.72, 1.27, 2.99, 1.08, 2.13, 1.00, 1.00, 1.78, 1.70, 0.00, 0.00, 1.84, 2.62}}
};


std::map<int,std::vector<float>> deltamv = {
  {2016, {0.20,0.10,0.08,0.11}},
  {2017, {0.19,0.09,0.07,0.08}},
  {2018, {0.19,0.09,0.07,0.08}},
};

float binBorders [10] = { 1.1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};



RooPlot* prepareFrame(RooPlot* frame){
    frame->GetYaxis()->SetTitleOffset(1.8);
    frame->SetMaximum(frame->GetMaximum()*1.15);
    frame->SetMinimum(0);
    return frame;
}


RooGaussian* constrainVar(RooRealVar* var,
			  string inVarName,
			  RooWorkspace *w,
			  int year,
			  bool addToList,
			  RooArgSet &c_vars,
			  RooArgSet &c_pdfs,
			  int scaleWidth = -1
                          ){

  float widthFac = 1.;
  if (scaleWidth>=0) widthFac = constr_width_fac[year][scaleWidth];

    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( widthFac*w->var(inVarName.c_str())->getError())
                                                 ); 
//     gauss_constr ->Print();
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
    return gauss_constr;                 
}

void constrainVar2(RooRealVar* var, 
		   string inVarName,
		   RooWorkspace *w,
		   int year,
		   bool addToList,
		   RooArgSet &c_vars,
		   RooArgSet &c_pdfs,
		   int scaleWidth = -1
		   ){

  float widthFac = 1.;
  if (scaleWidth>=0) widthFac = constr_width_fac[year][scaleWidth];

    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( widthFac*w->var(inVarName.c_str())->getError())
                                                 ); 
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
}


bool retrieveWorkspace(string filename, std::vector<RooWorkspace*> &ws, std::string ws_name){

    TFile* f =  TFile::Open( filename.c_str() ) ;
    if ( !f || !f->IsOpen() ) {
      cout << "File not found: " <<  gSystem->GetFromPipe(Form("readlink -m %s",filename.c_str())) << endl;
      return false;
    }else{
      cout << "Opening: " << gSystem->GetFromPipe(Form("readlink -m %s",filename.c_str())) << endl;
    }
    RooWorkspace* open_w = (RooWorkspace*)f->Get(ws_name.c_str());
    if ( !open_w || open_w->IsZombie() ) {
      cout<<"Workspace "<< ws_name <<  "not found in file: " << filename << endl;
      return false;
    }
    ws.push_back( open_w );
    f->Close();
    return true;
}


std::vector<RooDataSet*> createDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year, //std::map<int,float> scale_to_data,
                                       RooArgSet vars, std::string shortString  ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;

	RooDataSet* isample = new RooDataSet(*dataCT,("data_"+shortString + Form("_subs%i", is)).c_str());
        isample->append(*dataWT);
        datasample.push_back (isample);
      }
    }
    else{

      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;

      RooDataSet* isample = new RooDataSet( *dataCT, ("data_"+shortString + "_subs0").c_str() );
      isample->append(*dataWT);
      /* isample->append(*((RooDataSet*)dataCT->reduce(EventRange(0,0.1*dataCT->numEntries())))); */
      /* isample->append(*((RooDataSet*)dataWT->reduce(EventRange(0,0.1*dataWT->numEntries())))); */
      datasample.push_back (isample);
    }
    return datasample;
}

std::vector<RooDataSet*> createDatasetInData(RooWorkspace *ws, int q2Bin, RooArgSet vars, std::string shortString ){

  RooDataSet* data;
  std::vector<RooDataSet*> datasample;

  RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
				       ("data_"+shortString + "_subs0").c_str(), 
				       RooArgSet(vars));
  data = (RooDataSet*)ws->data(Form("data_b%i",q2Bin)) ;
  isample->append(*data);
  /* isample->append(*((RooDataSet*)data->reduce(EventRange(0,0.1*data->numEntries())))); */
  datasample.push_back (isample);

  return datasample;
}

std::vector<RooDataSet*> createDatasetInData(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws,
					     int q2Bin, int year, RooArgSet vars, std::string shortString, int q2st )
{

  std::vector<RooDataSet*> datasample;

  RooDataSet* data = (RooDataSet*)ws->data(Form("data_b%i",q2Bin)) ;

  if (nSample>0){

    int nSigCR[3] = {627365, 787149, 1694952};
    int subSampStat = (int)(data->sumEntries()*nsig_years[year][q2st]/nSigCR[year-2016]);

    for (uint is = firstSample; is <= lastSample; is++) {

      RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i", is)).c_str(),
					   ("data_"+shortString + Form("_subs%i", is)).c_str(),
					   RooArgSet(vars));

      isample->append(*((RooDataSet*)data->reduce(EventRange(is*subSampStat+1,
							     (is+1)*subSampStat))));
      datasample.push_back (isample);
    }
  }
  else{
    RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(),
					 ("data_"+shortString + "_subs0").c_str(),
					 RooArgSet(vars));
    isample->append(*data);
    datasample.push_back (isample);

  }


  return datasample;
}

RooGenericPdf* createPdfCuts(int q2Bin, int year, RooRealVar *var, RooRealVar *slo){
      if((q2Bin!=3&&q2Bin!=5&&q2Bin!=7) || (year<2016||year>2018)){
       std::cout<<Form("ERROR!!! ==> PdfCut not defined for q2Bin=%i, year=%i",q2Bin,year)<<std::endl;
       exit(1);
      }
      RooGenericPdf* PdfCut  = 0;
      string expression;
      double PDGB0Mass       = 5.2797;
      double PDGJpsiMass     = 3.0969;
      double PDGPsiPrimeMass = 3.6861;
      double deltam[4];
      deltam[0]=deltamv[year][0];
      deltam[1]=deltamv[year][1];
      deltam[2]=deltamv[year][2];
      deltam[3]=deltamv[year][3];
      double resmass[4] = {PDGJpsiMass,PDGJpsiMass,PDGPsiPrimeMass,PDGPsiPrimeMass};
      double binlow[4]  = {sqrt(binBorders[3]),sqrt(binBorders[5]),sqrt(binBorders[5]),sqrt(binBorders[7])};
      double binhigh[4] = {sqrt(binBorders[4]),sqrt(binBorders[6]),sqrt(binBorders[6]),sqrt(binBorders[8])};
      double m1[4];
      double m2[4];
      double m3[4];
      double m4[4];
      string formula[4] = {"","","",""};

      for (int i=0; i<4; ++i) {
        m1[i] = PDGB0Mass-resmass[i]+binhigh[i]+deltam[i];
        m2[i] = PDGB0Mass-resmass[i]+binlow[i]+deltam[i];
        m3[i] = PDGB0Mass-resmass[i]+binhigh[i]-deltam[i];
        m4[i] = PDGB0Mass-resmass[i]+binlow[i]-deltam[i];
        formula[i] = formula[i] + Form(" + (@0<%.5f)*(@0>%.5f)*(@0-%.5f)/%.5f",m1[i],m2[i],m2[i],binhigh[i]-binlow[i]);
        formula[i] = formula[i] + Form(" + (@0<%.5f)*(@0>%.5f)*(%.5f-@0)/%.5f",m3[i],m4[i],m3[i],binhigh[i]-binlow[i]);
      }
//      
     if(q2Bin==3)  expression = Form("(@0<%.5f)+(@0>%.5f)",m4[0],m1[0])+formula[0];
     if(q2Bin==5)  expression = Form("(@0<%.5f)+(@0>%.5f)*(@0<%.5f)+(@0>%.5f)",m4[1],m1[1],m4[2],m1[2])+formula[1]+formula[2];
     if(q2Bin==7)  expression = Form("(@0<%.5f)+(@0>%.5f)",m4[3],m1[3])+formula[3];
     expression="("+expression+")*exp(@0*@1)";
     std::cout<<Form("setting delatm for Era=%d => %f, %f, %f, %f",year,deltam[0],deltam[1],deltam[2],deltam[3])<<std::endl;
     std::cout<<Form("PdfCut expression defined for q2Bin=%i, year=%i as: %s",q2Bin,year,expression.c_str())<<std::endl;
     
//     
     PdfCut = new RooGenericPdf(Form("pdf%i_%i",q2Bin,year),Form("pdf%i_%i",q2Bin,year),expression.c_str(),RooArgList(*var,*slo));
     PdfCut->Print();
     return PdfCut;
}
