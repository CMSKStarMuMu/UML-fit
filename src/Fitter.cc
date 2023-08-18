#include "Riostream.h" 

#include "Fitter.h" 

ClassImp(Fitter) 

Fitter::Fitter(const char *_name, const char *_title,
	       RooArgList _angPars,
	       RooDataSet* _combData,
	       RooAbsPdf* _simPdf,
	       RooAbsPdf* _simPdf_penalty,
	       BoundCheck* _boundary,
	       BoundDist* _bound_dist,
	       Penalty* _penTerm,
               RooArgSet* _constrVars
	       ) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist),
  penTerm(_penTerm),
  constrVars(_constrVars)
{

  SetDefConf();
  usedPenalty = false;
  runSimpleFit = false;

}

Fitter::Fitter(const char *_name, const char *_title,
	       RooArgList _angPars,
	       RooDataSet* _combData,
	       RooAbsPdf* _simPdf,
	       RooAbsPdf* _simPdf_penalty,
	       BoundCheck* _boundary,
	       BoundDist* _bound_dist,
      	       RooArgSet* _constrVars 	       
	       ) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist),
  constrVars(_constrVars)
{

  penTerm = 0;

  SetDefConf();
  usedPenalty = false;
  runSimpleFit = false;

}

Fitter::Fitter(const Fitter& other, const char* name) :  
  name(other.name),
  title(other.title),
  angPars(other.angPars),
  combData(other.combData),
  simPdf(other.simPdf),
  simPdf_penalty(other.simPdf_penalty),
  boundary(other.boundary),
  bound_dist(other.bound_dist),
  penTerm(other.penTerm),
  constrVars(other.constrVars)
{


  fac1 = other.fac1;
  fac4 = other.fac4;
  base1 = other.base1;
  base4 = other.base4;
  max1 = other.max1;
  max4 = other.max4;
  base1_corr = other.base1_corr;
  base4_corr = other.base4_corr;
  min_base = other.min_base;

  maxCoeff = other.maxCoeff;
  coeff1 = other.coeff1;
  coeff4 = other.coeff4;
  coeff5 = other.coeff5;

  usedPenalty = other.usedPenalty;
  runSimpleFit = other.runSimpleFit;

}

void Fitter::SetDefConf()
{
  fac1 = 0.1;
  fac4 = 0.01;
  base1 = 0.05;
  base4 = 0.05;
  max1 = 0.;
  max4 = 0.;

  min_base = 1.05;

  base1_corr = base1*sqrt(combData->numEntries());
  base4_corr = base4*sqrt(combData->numEntries());
  if (base1_corr<min_base) base1_corr = min_base;
  if (base4_corr<min_base) base4_corr = min_base;

  maxCoeff = 1e8;
  coeff1 = 0;
  coeff4 = 0;
  coeff5 = 0;

  boundDist = -1.;
  boundDistTime = 0.;

  minParError = 0.01;
  widthScale = 0.1;

  vFreeFitResult = std::vector<Double_t>(angPars.getSize(),0);
  vFreeFitErrLow  = std::vector<Double_t>(angPars.getSize(),0);
  vFreeFitErrHigh = std::vector<Double_t>(angPars.getSize(),0);
  vFitResult = std::vector<Double_t>(angPars.getSize(),0);
  vFitErrLow  = std::vector<Double_t>(angPars.getSize(),0);
  vFitErrHigh = std::vector<Double_t>(angPars.getSize(),0);
  vImprovResult = std::vector<Double_t>(angPars.getSize(),0);

  vResult = std::vector<Double_t>(angPars.getSize(),0);
  vConfInterLow  = std::vector<Double_t>(angPars.getSize(),0);
  vConfInterHigh = std::vector<Double_t>(angPars.getSize(),0);

  nCPU = 1;
  nCPU_Pen=1;

  unbl = 0;

}


Int_t Fitter::fit()
{

    coeff1 = 0;
    coeff4 = 0;
    coeff5 = 0;

    // set up free fit
    nll = simPdf->createNLL(*combData,
                            RooFit::Extended(kFALSE),
                            RooFit::NumCPU(nCPU)
                            );
    if ( constrVars != nullptr && constrVars->getSize() > 0  ){
      nll = simPdf->createNLL(*combData,
                              RooFit::Extended(kFALSE),
                              RooFit::NumCPU(nCPU),
                              RooFit::Constrain(*constrVars)
                              );
    }                              
         
    RooMinimizer m(*nll) ;
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    // m.setVerbose(kTRUE);
    m.setPrintLevel(-1);
    m.setPrintEvalErrors(-1);
    //  Minuit2.setEps(1e-16) ;
    m.setMinimizerType("Minuit2");

    // free fit
    m.setStrategy(0);
    // m.setEvalErrorWall(false);
    m.migrad() ;
    m.hesse() ;
    // std::cout << std::endl;
    // std::cout << "######################### now strategy 2 #########################"<< std::endl;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;
    // m.minos() ;
    
    result_free = m.save("result") ;

    result_penalty = 0;
    usedPenalty = false;

    std::cout<<"============ FREE FIT ==========="<<std::endl;
    if (unbl>3) result_free->Print("v");

    // if free fit is good return its result
    if ( result_free->status()==0 && result_free->covQual()==3 && boundary->getValV() == 0 ) {
      if ( !runSimpleFit ) computeBoundaryDistance(); //allow to skip the boundary-distance computation (for full-stat MC and control region fits)
      fillResultContainers(false,true);
      return 0;
    }

    if ( runSimpleFit ) { //allow to skip the penalised fit (for full-stat MC and control region fits)
      if ( result_free->status()==0 && result_free->covQual()==3 ) {
	fillResultContainers(false,true);
	return 1;
      }

      return 2;

    }
    
    std::cout<<"============ START FIT WITH PENALTY ==========="<<std::endl;
    usedPenalty = true;

    // optional: if a partial boundary is satisfied
    // do not apply the corresponding penalty term
    // bool inCTL4  = true;
    // bool inCTL15 = true;
    // if ( !boundary->isInCTL4()  ) inCTL4  = false;
    // if ( !boundary->isInCTL15() ) inCTL15 = false;

    for (int totCoeff=0; fac1*pow(base1_corr,totCoeff)<=maxCoeff; ++totCoeff) {

      for (int iCoeff1=totCoeff; iCoeff1>=0; --iCoeff1) {

	// set penalty coefficients
	coeff1 = fac1 * pow(base1_corr,iCoeff1);
	if (max1>0 && coeff1>max1) continue;

	coeff4 = fac4 * pow(base4_corr,totCoeff-iCoeff1);
	if (max4>0 && coeff4>max4) continue;

	coeff5 = pow(coeff1,1.5) / 316.2;

	// optional: if a partial boundary is satisfied
	// do not apply the corresponding penalty term
	// if ( inCTL15 ) {
	//   if ( iCoeff1>0 ) continue;
	//   coeff1 = 0;
	//   coeff5 = 0;
	// }
	// if ( inCTL4 ) {
	//   if ( totCoeff-iCoeff1>0 ) continue;
	//   coeff4 = 0;
	// }

	penTerm->setCoefficient(1,coeff1);
	penTerm->setCoefficient(4,coeff4);
	penTerm->setCoefficient(5,coeff5);

	// set up the penalised fit
	nll_penalty = simPdf_penalty->createNLL(*combData,
						RooFit::Extended(kFALSE),
						RooFit::NumCPU(nCPU_Pen)
						);

	RooMinimizer m_penalty (*nll_penalty) ;
	m_penalty.optimizeConst(kTRUE);
	m_penalty.setOffsetting(kTRUE);
	// m_penalty.setVerbose(kTRUE);
	m_penalty.setMinimizerType("Minuit2");
	// m_penalty.setProfile(kTRUE);
	m_penalty.setPrintLevel(-1);
	m_penalty.setPrintEvalErrors(-1);
	m_penalty.setStrategy(2);
    
	// penalised fit
	m_penalty.migrad() ;
	m_penalty.hesse() ;
	result_penalty = m_penalty.save("result");
	    
	// std::cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<std::endl;
	if (unbl>3) result_penalty->Print("v");

	// if a good fit is found return good status
	if ( result_penalty->status()==0 && result_penalty->covQual()==3 ) {
	  if ( boundary->getValV()==0 ) {
	    if (unbl>0) std::cout<<"P "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<std::endl;
	    computeBoundaryDistance();
	    fillResultContainers();
	    return 0;
	  } else
	    if (unbl>0) std::cout<<"O "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<std::endl;
	} else
	  if (unbl>0) std::cout<<"N "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<std::endl;

      }

    }
    
    // if no good fit is found return bad status
    return 2;

}

Int_t Fitter::improveAng(int seed, int nGen)
{

  Double_t preBoundDist = boundDist;
  
  // Improve global fit result
  TRandom3 randGen (seed);

  std::vector<double> vTestPar(angPars.getSize());
  std::vector<double> vImprovPar(angPars.getSize());
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    vImprovPar[iPar] = ((RooRealVar*)angPars.at(iPar))->getValV();

  double NLL_before = nll->getValV();
  double improvNLL = NLL_before;
  double testNLL = 0;
  int iImprove = 0;
  std::cout << Form("NLL_before = %.5f",NLL_before) << std::endl;

  do {

    for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)angPars.at(iPar);
      do vTestPar[iPar] = randGen.Gaus(vImprovPar[iPar],TMath::Max(boundDist,0.002));
      while (vTestPar[iPar]>par->getMax() || vTestPar[iPar]<par->getMin());
      par->setVal(vTestPar[iPar]);
    }

    if (boundary->getValV()>0) continue;

    testNLL = nll->getValV();
    if (improvNLL>testNLL) {
      improvNLL = testNLL;
      for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
	vImprovPar[iPar] = vTestPar[iPar];
    }

    ++iImprove;

  } while (iImprove<nGen);

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    ((RooRealVar*)angPars.at(iPar))->setVal(vImprovPar[iPar]);
  
  computeBoundaryDistance();
  std::cout << Form("NLL_after = %.5f",improvNLL) << std::endl;
  std::cout<<"Improved fit result: deltaNLL = "<<NLL_before-improvNLL<<" bound dist: "<<preBoundDist<<" -> "<<boundDist<<std::endl;

  fillResultContainers(true);

  return 0;

}


Int_t Fitter::MinosAng(int seed, int nGenMINOS)
{

  // NLL of the best-fit result
  double NLL_min = nll->getValV();

  TRandom3 randGenMinos (seed);
  double probedNLL;

  // Loop over the parameters
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)angPars.at(iPar);

    // get and print the best-fit result
    double p_best = vResult[iPar];
    if (unbl>3) std::cout<<par->GetName()<<" best: "<<p_best<<std::endl;

    // vectors for TGraph plots
    // std::vector<double> vPval (0);
    // std::vector<double> vdNLL (0);
    // if (is==0) {
    //   vPval.push_back(p_best);
    //   vdNLL.push_back(0);
    // }

    // firstly low, then high error
    for (int isErrHigh=0; isErrHigh<2; ++isErrHigh) {

      std::vector<double> vLastHit(0);
      for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1)
	vLastHit.push_back(vResult[iPar1]);

      TH1D* parRandomPool = 0;
      int nHistBins = 0;
      if (isErrHigh>0) {
	nHistBins = (int)((par->getMax()-p_best)/0.0005);
	// nHistBins = (int)((par->getMax()-p_best+0.0005)/0.001);
	if (nHistBins<1) nHistBins=1;
	parRandomPool = new TH1D(Form("hRandPoolH%i_%s",iPar,name.Data()),
				 Form("hRandPoolH%i_%s",iPar,name.Data()),
				 nHistBins,par->getMax()-0.0005*nHistBins,par->getMax());
	// nHistBins,par->getMax()+0.0005-0.001*nHistBins,par->getMax()+0.0005);
      } else {
	nHistBins = (int)((p_best-par->getMin())/0.0005);
	// nHistBins = (int)((p_best-par->getMin()+0.0005)/0.001);
	if (nHistBins<1) nHistBins=1;
	parRandomPool = new TH1D(Form("hRandPoolL%i_%s",iPar,name.Data()),
				 Form("hRandPoolL%i_%s",iPar,name.Data()),
				 nHistBins,par->getMin(),par->getMin()+0.0005*nHistBins);
	// nHistBins,par->getMin()-0.0005,par->getMin()-0.0005+0.001*nHistBins);
      }
      double sigma = fabs(isErrHigh>0?vFitErrHigh[iPar]:vFitErrLow[iPar]);
      if (sigma<minParError) sigma = minParError;
      for (int iBin=1; iBin<=nHistBins; ++iBin) {
	double x = (parRandomPool->GetBinCenter(iBin)-p_best) / sigma;
	parRandomPool->SetBinContent(iBin,fabs(x)*exp(-0.5*x*x));
      }

      double p_in = p_best;
      double p_test = 0;

      int iPnt=0;

      do {

	do p_test = parRandomPool->GetRandom();
	while (p_test>par->getMax() || p_test<par->getMin());
	par->setVal(p_test);
	for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
	  if (iPar1==iPar) continue;
	  RooRealVar* par1 = (RooRealVar*)angPars.at(iPar1);
	  double par1val = 0;
	  do par1val = randGenMinos.Gaus(vLastHit[iPar1],widthScale*TMath::Max(0.5*(vFitErrHigh[iPar1]-vFitErrLow[iPar1]),minParError));
	  while (par1val>par1->getMax() || par1val<par1->getMin());
	  par1->setVal(par1val);
	}
	// check if the point is physical
	if (boundary->getValV()>0) continue;
	// get and test the local likelihood
	probedNLL = nll->getValV();
	if (probedNLL<=NLL_min+0.5) {
	  p_in = p_test;
	  if ( isErrHigh > 0 ) { if ( p_in > par->getMax()-parRandomPool->GetBinWidth(1) ) break; }
	  else if ( p_in < par->getMin()+parRandomPool->GetBinWidth(1) ) break;
	  for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
	    RooRealVar* par1 = (RooRealVar*)angPars.at(iPar1);
	    vLastHit[iPar1] = par1->getValV();
	  }
	  if (isErrHigh>0)
	    for (int iBin=1; iBin<=parRandomPool->FindBin(p_test); ++iBin)
	      parRandomPool->SetBinContent(iBin,0);
	  else 
	    for (int iBin=parRandomPool->FindBin(p_test); iBin<=nHistBins; ++iBin)
	      parRandomPool->SetBinContent(iBin,0);
	} else {
	  parRandomPool->Fill(p_test,0.02/(probedNLL-NLL_min-0.5));
	}
	  
	// fill the plotting vectors
	// if (is==0 && probedNLL-NLL_min<4.5) {
	//   vPval.push_back(p_test);
	//   vdNLL.push_back(probedNLL-NLL_min);
	// }

	++iPnt;
	// apply conditions
      } while ( iPnt < nGenMINOS );

      if (isErrHigh>0) {
	vConfInterHigh[iPar] = p_in;
	if (unbl>3) std::cout<<par->GetName()<<" high: "<<p_in<<std::endl;
      } else {
	vConfInterLow[iPar] = p_in;
	if (unbl>3) std::cout<<par->GetName()<<" low:  "<<p_in<<std::endl;
      }

    }

    // if (is==0) {
    //   // produce deltaNLL vs parameter graph, with the probed points
    //   TCanvas* canNLL = new TCanvas(Form("canNLL_%s",par->GetName()),"canNLL",1000,1000);
    //   TGraph* grNLL = new TGraph(vPval.size(),&vPval[0],&vdNLL[0]);
    //   grNLL->SetName(Form("grNLL_%s",par->GetName()));
    //   grNLL->SetTitle(Form("deltaNLL scan for %s",par->GetTitle()));
    //   grNLL->GetXaxis()->SetTitle(par->GetName());
    //   grNLL->GetYaxis()->SetTitle("deltaNLL");
    //   grNLL->SetMarkerStyle(7);
    //   // grNLL->SetMarkerStyle(20);
    //   // grNLL->SetMarkerSize(2);
    //   grNLL->SetMarkerColor(9);
    //   canNLL->cd();
    //   grNLL->Draw("AP");
	  
    //   TLine errLow (vConfInterLow[iPar],grNLL->GetYaxis()->GetXmin(),vConfInterLow[iPar],grNLL->GetYaxis()->GetXmax());
    //   TLine errHigh(vConfInterHigh[iPar],grNLL->GetYaxis()->GetXmin(),vConfInterHigh[iPar],grNLL->GetYaxis()->GetXmax());
    //   TLine DeltaNLL0p5 (grNLL->GetXaxis()->GetXmin(),0.5,grNLL->GetXaxis()->GetXmax(),0.5);
    //   errLow .SetLineColor(46);
    //   errHigh.SetLineColor(46);
    //   DeltaNLL0p5.SetLineColor(13);
    //   DeltaNLL0p5.SetLineStyle(9);
    //   errLow .Draw();
    //   errHigh.Draw();
    //   DeltaNLL0p5.Draw();

    //   canNLL->SaveAs(Form("plotSimFit_d/profiledNLL-%s_toy%i_%s_%s.pdf",par->GetName(),seed,shortString.c_str(),all_years.c_str()));
    // }

  }

  return 0;

}


void Fitter::fillResultContainers(bool fromImprov, bool isFree)
{

  // fill results' containers

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)angPars.at(iPar);
    vResult[iPar] = par->getValV();

    if (!fromImprov) {
      if (isFree) {
	vFreeFitResult[iPar] = vResult[iPar];
	vFreeFitErrLow[iPar] = par->getErrorLo();
	vFreeFitErrHigh[iPar] = par->getErrorHi();
	vFitResult[iPar] = vResult[iPar];
	vFitErrLow[iPar] = vFreeFitErrLow[iPar];
	vFitErrHigh[iPar] = vFreeFitErrHigh[iPar];
      } else {
	vFitResult[iPar] = vResult[iPar];
	vFitErrLow[iPar] = par->getErrorLo();
	vFitErrHigh[iPar] = par->getErrorHi();
      }
    } else vImprovResult[iPar] = vResult[iPar];

  }

}


Double_t Fitter::computeBoundaryDistance()
{

  TStopwatch distTime;
  distTime.Start(true);
  std::cout<<"== Start to compute  boundary distance =="<<std::endl;

  // Compute distance from boundary
  boundDist = bound_dist->getValV();

  distTime.Stop();
  // if the result was not cached (producing time=~0) update the time
  if (distTime.CpuTime()>0.01)
    boundDistTime = distTime.CpuTime();

  std::cout<<"Distance from boundary: "<<boundDist<<" (computed in "<<boundDistTime<<" s)"<<std::endl;
  
  return boundDist;

}

void Fitter::plotSimFitProjections(const char* filename, std::vector<std::string> catnames, std::vector<int> years, bool is4D)
{

  std::cout<<Form("== Start filling plots in : %s ==",filename)<<std::endl;
  auto canv = new TCanvas ("canv","canv",is4D?2000:1500,500*years.size());
  canv->Divide(is4D?4:3, years.size());
  std::vector<std::string> catnamesresolv(catnames.size());
  
  for (unsigned int iy = 0; iy < years.size(); iy++) {

    for (uint i=0; i<catnames.size(); ++i)
      catnamesresolv[i] = Form(catnames[i].c_str(),years[iy]);

    auto singleYearPdf = ((RooSimultaneous*)simPdf)->getPdf(catnamesresolv[0].c_str());
    auto singleYearData = combData->reduce(("sample==sample::"+catnamesresolv[0]).c_str());

    std::string frametitle = "%s projection of ";
    frametitle = frametitle + Form("%i dataset",years[iy]);

    plotProjections(singleYearPdf, singleYearData, catnamesresolv, frametitle, is4D, canv, iy*(is4D?4:3));

  }

  canv->SaveAs(filename);
  delete canv;

}

void Fitter::plotProjections(RooAbsPdf* singleYearPdf, RooAbsData* singleYearData, std::vector<std::string> catnames, std::string frametitle, bool is4D, TCanvas* canv, int ipad)
{

    TLegend* leg = new TLegend (0.60,0.75,0.9,0.9);

    std::vector<std::string> varname(0);
    if (is4D) varname = {"mass", "ctK", "ctL", "phi"};
    else varname = {"ctK", "ctL", "phi"};

    std::vector<RooPlot*> frames;
    TLatex Tl;
    Tl.SetTextSize(0.04);

    for (unsigned int fr = 0; fr < varname.size(); fr++){

      auto var = (RooRealVar*)singleYearPdf->getObservables(singleYearData)->find(varname[fr].c_str());
      frames.push_back( var->frame(RooFit::Title(Form(frametitle.c_str(),var->GetTitle()))) );
    
      singleYearData->plotOn(frames[fr],
			     RooFit::MarkerColor(kRed+1),
			     RooFit::LineColor(kRed+1),
			     RooFit::Binning(40),
			     RooFit::Name(Form("plData%i",ipad)) );
        
      singleYearPdf->plotOn(frames[fr],
			    RooFit::LineWidth(1),
			    RooFit::Name(Form("plPDF%i",ipad)),
			    RooFit::NumCPU(4));

      if (catnames.size()>1) {
	singleYearPdf->plotOn(frames[fr],
			      RooFit::LineWidth(1),
			      RooFit::Name(Form("plPDFbkg%i",ipad)),
			      RooFit::NumCPU(4),
			      RooFit::LineColor(8),
			      RooFit::Components( catnames[2].c_str() ));

	singleYearPdf->plotOn(frames[fr],
			      RooFit::LineWidth(1),
			      RooFit::Name(Form("plPDFsig%i",ipad)),
			      RooFit::NumCPU(4),
			      RooFit::LineColor(880),
			      RooFit::Components( catnames[1].c_str() ));
      }

      if (fr == 0) { 
	leg->AddEntry(frames[fr]->findObject(Form("plData%i",ipad)),
		      "Data",	"lep");
	if (catnames.size()>1) {
	  leg->AddEntry(frames[fr]->findObject(Form("plPDF%i",ipad)),
			"Total PDF","l");
	  leg->AddEntry(frames[fr]->findObject(Form("plPDFsig%i",ipad)),
			"Signal","l");
	  leg->AddEntry(frames[fr]->findObject(Form("plPDFbkg%i",ipad)),
			"Background","l");
	} else
	  leg->AddEntry(frames[fr]->findObject(Form("plPDF%i",ipad)),
			"Signal PDF","l");
      }

      double chiSq = frames[fr]->chiSquare(Form("plPDF%i",ipad),Form("plData%i",ipad));

      canv->cd(ipad+fr+1);
      frames[fr]->GetYaxis()->SetTitleOffset(2);
      gPad->SetLeftMargin(0.16); 
      gPad->SetRightMargin(0.04); 
      frames[fr]->Draw();
      Tl.DrawLatexNDC(0.2,0.86,Form("#chi^{2}/N_{bins} = %.2f",chiSq));
      if (fr == 0) leg->Draw("same");
    }
}
