/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "PdfSigAng.h" 
#include <math.h>
#include "TMath.h"

ClassImp(PdfSigAng) 

PdfSigAng::PdfSigAng(const char *name, const char *title, 
		     RooAbsReal& _ctK,
		     RooAbsReal& _ctL,
		     RooAbsReal& _phi,
		     RooAbsReal& _Fl,
		     RooAbsReal& _P1,
		     RooAbsReal& _P2,
		     RooAbsReal& _P3,
		     RooAbsReal& _P4p,
		     RooAbsReal& _P5p,
		     RooAbsReal& _P6p,
		     RooAbsReal& _P8p,
		     RooAbsReal& _mFrac,
                     RooAbsReal& _rtAngTerm,
                     RooAbsReal& _wtAngTerm
		     ) :
  RooAbsPdf(name,title), 
  ctK("ctK","ctK",this,_ctK),
  ctL("ctL","ctL",this,_ctL),
  phi("phi","phi",this,_phi),
  Fl("Fl","Fl",this,_Fl),
  P1("P1","P1",this,_P1),
  P2("P2","P2",this,_P2),
  P3("P3","P3",this,_P3),
  P4p("P4p","P4p",this,_P4p),
  P5p("P5p","P5p",this,_P5p),
  P6p("P6p","P6p",this,_P6p),
  P8p("P8p","P8p",this,_P8p),
  mFrac("mFrac","mFrac",this,_mFrac),
  rtAngTerm("rtAngTerm","rtAngTerm",this,_rtAngTerm),
  wtAngTerm("wtAngTerm","wtAngTerm",this,_wtAngTerm)
{

  isPenalised = false;

}

PdfSigAng::PdfSigAng(const char *name, const char *title, 
		     RooAbsReal& _ctK,
		     RooAbsReal& _ctL,
		     RooAbsReal& _phi,
		     RooAbsReal& _Fl,
		     RooAbsReal& _P1,
		     RooAbsReal& _P2,
		     RooAbsReal& _P3,
		     RooAbsReal& _P4p,
		     RooAbsReal& _P5p,
		     RooAbsReal& _P6p,
		     RooAbsReal& _P8p,
		     RooAbsReal& _mFrac,
                     RooAbsReal& _rtAngTerm,
                     RooAbsReal& _wtAngTerm,
		     RooAbsReal& _PenTerm) :
  RooAbsPdf(name,title), 
  ctK("ctK","ctK",this,_ctK),
  ctL("ctL","ctL",this,_ctL),
  phi("phi","phi",this,_phi),
  Fl("Fl","Fl",this,_Fl),
  P1("P1","P1",this,_P1),
  P2("P2","P2",this,_P2),
  P3("P3","P3",this,_P3),
  P4p("P4p","P4p",this,_P4p),
  P5p("P5p","P5p",this,_P5p),
  P6p("P6p","P6p",this,_P6p),
  P8p("P8p","P8p",this,_P8p),
  mFrac("mFrac","mFrac",this,_mFrac),
  rtAngTerm("rtAngTerm","rtAngTerm",this,_rtAngTerm),
  wtAngTerm("wtAngTerm","wtAngTerm",this,_wtAngTerm),
  PenTerm("PenTerm","PenTerm",this,_PenTerm)
{

  isPenalised = true;

}


PdfSigAng::PdfSigAng(const PdfSigAng& other, const char* name) :  
  RooAbsPdf(other,name), 
  ctK("ctK",this,other.ctK),
  ctL("ctL",this,other.ctL),
  phi("phi",this,other.phi),
  Fl("Fl",this,other.Fl),
  P1("P1",this,other.P1),
  P2("P2",this,other.P2),
  P3("P3",this,other.P3),
  P4p("P4p",this,other.P4p),
  P5p("P5p",this,other.P5p),
  P6p("P6p",this,other.P6p),
  P8p("P8p",this,other.P8p),
  mFrac("mFrac",this,other.mFrac),
  rtAngTerm("rtAngTerm", this, other.rtAngTerm),
  wtAngTerm("wtAngTerm", this, other.wtAngTerm)
{

  if (other.isPenalised) {
    PenTerm = RooRealProxy("PenTerm",this,other.PenTerm);
    isPenalised = true;
  } else isPenalised = false;

}



Double_t PdfSigAng::evaluate() const 
{
  double penalty = 1;
  if (isPenalised) penalty = penTermVal()->getVal();

  double decCT_times_eff = ((RooAbsReal&)(rtAngTerm.arg())).getVal();
  double decWT_times_eff = ((RooAbsReal&)(wtAngTerm.arg())).getVal();

  double ret = (decCT_times_eff + mFrac * decWT_times_eff ) * penalty;
  std::cout<<"eval: "<<decCT_times_eff<<"\t"<<mFrac<<"\t"<<decWT_times_eff<<"\t"<<penalty<<std::endl;
  return ret;
}

namespace {
  Bool_t fullRangeCosT(const RooRealProxy& x ,const char* range)
  {
    // set accepted integration range for cosTheta variables
    return range == 0 || strlen(range) == 0
      ? std::fabs(x.min() + 1.) < 1.e-5 && std::fabs(x.max() - 1.) < 1.e-5
      : std::fabs(x.min(range) + 1.) < 1.e-5 && std::fabs(x.max(range) - 1.) < 1.e-5;
  }
  Bool_t fullRangePhi(const RooRealProxy& x ,const char* range)
  {
    // set accepted integration range for phi variable
    return range == 0 || strlen(range) == 0
      ? std::fabs(x.min() + TMath::Pi()) < 1.e-3 && std::fabs(x.max() - TMath::Pi()) < 1.e-3
      : std::fabs(x.min(range) + TMath::Pi()) < 1.e-3 && std::fabs(x.max(range) - TMath::Pi()) < 1.e-3;
  }
}

Int_t PdfSigAng::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // use pre-computed integrals for the integration over all the three variables
  // after checking that integrals are in place
//   if ( intCPart.size()<1 || intCPart[0]==0 || intWPart.size()<1 || intWPart[0]==0 )
//     return 0 ;
  if ( matchArgs(allVars,analVars,ctK,ctL,phi) ){
    if ( fullRangeCosT(ctK,rangeName) && fullRangeCosT(ctL,rangeName) && fullRangePhi(phi,rangeName) ){
      // std::cout << "code 1"<<  std::endl;
      return 1 ;
    }  
  }      
  if ( matchArgs(allVars,analVars,ctL,phi) ){
    if ( fullRangeCosT(ctL,rangeName) && fullRangePhi(phi,rangeName) ){
      // std::cout << "code 2"<<  std::endl;
      return 2 ;      
    }  
  }
  if ( matchArgs(allVars,analVars,ctK,phi) ){
    if ( fullRangeCosT(ctK,rangeName) && fullRangePhi(phi,rangeName) ){
      // std::cout << "code 3"<<  std::endl;
      return 3 ;      
    }  
  }
  if ( matchArgs(allVars,analVars,ctK,ctL) ){
    if ( fullRangeCosT(ctK,rangeName) && fullRangeCosT(ctL,rangeName) ){
      // std::cout << "code 4"<<  std::endl;
      return 4 ;      
    }  
  }

  // the lack of analytical integral for the subsets of angular variables does not slow down the fit
  // since only the complete integration is used there
  // if one wants to speed up also the PDF projection for plotting, the other analytical integrals can be computed
  // but it seems a huge effort to me
  return 0 ;

}

Double_t PdfSigAng::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code>0 && code<5) ;

    double theIntegral;

    RooAbsReal & ctKarg = (RooAbsReal&)ctK.arg();
    RooAbsReal & ctLarg = (RooAbsReal&)ctL.arg();
    RooAbsReal & phiarg = (RooAbsReal&)phi.arg();

    RooAbsReal & rtAng = (RooAbsReal&)rtAngTerm.arg();
    RooAbsReal & wtAng = (RooAbsReal&)wtAngTerm.arg();
    
    double rtAngIntegral, wtAngIntegral;

    if (code ==1){
  
      rtAngIntegral = ((RooAbsReal* )rtAng.createIntegral(RooArgSet(ctKarg,ctLarg,phiarg)))->getVal();
      wtAngIntegral = ((RooAbsReal* )wtAng.createIntegral(RooArgSet(ctKarg,ctLarg,phiarg)))->getVal();
  
      if (rtAngIntegral<=0) {
          if (rtAngIntegral<0) std::cout<<"ERROR! Negative ct pdf integral, fake value returned"<<std::endl;
          else std::cout<<"ERROR! Null ct pdf integral, fake value returned"<<std::endl;
          return 1e-55;
      }
      if (wtAngIntegral<=0) {
          if (wtAngIntegral<0) std::cout<<"ERROR! Negative wt pdf integral, fake value returned"<<std::endl;
          else std::cout<<"ERROR! Null wt pdf integral, fake value returned"<<std::endl;
          return 1e-55;
      }

      std::cout<<"3D int: "<<rtAngIntegral<<"\t"<<wtAngIntegral<<std::endl;
      // std::cout <<  "PdfSigAng:analyticalIntegral1:\t" << rtAngIntegral  << "\t" << wtAngIntegral << std::endl;
    }   
    else if (code ==2){
      //matchArgs(allVars,analVars,ctL,phi)
      rtAngIntegral = ((RooAbsReal* )rtAng.createIntegral( RooArgSet(ctLarg,phiarg) ))->getVal();
      wtAngIntegral = ((RooAbsReal* )wtAng.createIntegral( RooArgSet(ctLarg,phiarg) ))->getVal();
      // std::cout <<  "PdfSigAng:analyticalIntegral2:\t" << rtAngIntegral  << "\t" << wtAngIntegral << std::endl;
    }
    else if (code ==3){
      //matchArgs(allVars,analVars,ctK,phi)
      rtAngIntegral = ((RooAbsReal* )rtAng.createIntegral( RooArgSet(ctKarg,phiarg) ))->getVal();
      wtAngIntegral = ((RooAbsReal* )wtAng.createIntegral( RooArgSet(ctKarg,phiarg) ))->getVal();
      // std::cout <<  "PdfSigAng:analyticalIntegral3\t" << rtAngIntegral  << "\t" << wtAngIntegral << std::endl;
    }
    else if (code ==4){
      //matchArgs(allVars,analVars,ctK,ctL)
      rtAngIntegral = ((RooAbsReal* )rtAng.createIntegral( RooArgSet(ctKarg,ctLarg) ))->getVal();
      wtAngIntegral = ((RooAbsReal* )wtAng.createIntegral( RooArgSet(ctKarg,ctLarg) ))->getVal();
      // std::cout <<  "PdfSigAng:analyticalIntegral4\t" << rtAngIntegral  << "\t" << wtAngIntegral << std::endl;
    }

//       theIntegral =        rtAngIntegral / ((RooAbsReal* )rtAng.createIntegral(RooArgSet(ctKarg,ctLarg,phiarg)))->getVal() + \
//                      mFrac*wtAngIntegral / ((RooAbsReal* )wtAng.createIntegral(RooArgSet(ctKarg,ctLarg,phiarg)))->getVal();
    theIntegral =  (rtAngIntegral + mFrac*wtAngIntegral);
    return theIntegral ;
}
