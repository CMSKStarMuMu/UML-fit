/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "PdfSigAng_Pen.h" 

ClassImp(PdfSigAng_Pen) 

PdfSigAng_Pen::PdfSigAng_Pen(const char *name, const char *title, 
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
			     RooAbsReal& _EffC,
			     RooAbsReal& _EffW,
			     std::vector<double> _intCPart,
			     std::vector<double> _intWPart,
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
  EffC("EffC","corr-tag efficiency",this,_EffC),
  EffW("EffW","wrong-tag efficiency",this,_EffW),
  intCPart(_intCPart),
  intWPart(_intWPart),
  PenTerm("PenTerm","PenTerm",this,_PenTerm)
{
}


PdfSigAng_Pen::PdfSigAng_Pen(const PdfSigAng_Pen& other, const char* name) :
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
  EffC("EffC",this,other.EffC),
  EffW("EffW",this,other.EffW),
  intCPart(other.intCPart),
  intWPart(other.intWPart),
  PenTerm("PenTerm",this,other.PenTerm)
{
}



Double_t PdfSigAng_Pen::evaluate() const 
{

  double decCT = ( 0.75 * (1-Fl) * (1-ctK*ctK) +
		   Fl * ctK*ctK +
		   ( 0.25 * (1-Fl) * (1-ctK*ctK) - Fl * ctK*ctK ) * ( 2 * ctL*ctL -1 ) +
		   0.5 * P1 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * cos(2*phi) +
		   2 * cos(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( P4p * ctL * sqrt(1-ctL*ctL) + P5p * sqrt(1-ctL*ctL) ) +
		   2 * sin(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( P8p * ctL * sqrt(1-ctL*ctL) - P6p * sqrt(1-ctL*ctL) ) +
		   2 * P2 * (1-Fl) * (1-ctK*ctK) * ctL -
		   P3 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * sin(2*phi) );

  double decWT = ( 0.75 * (1-Fl) * (1-ctK*ctK) +
		   Fl * ctK*ctK +
		   ( 0.25 * (1-Fl) * (1-ctK*ctK) - Fl * ctK*ctK ) * ( 2 * ctL*ctL -1 ) +
		   0.5 * P1 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * cos(2*phi) -
		   2 * cos(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( -1. * P4p * ctL * sqrt(1-ctL*ctL) + P5p * sqrt(1-ctL*ctL) ) +
		   2 * sin(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( -1. * P8p * ctL * sqrt(1-ctL*ctL) - P6p * sqrt(1-ctL*ctL) ) -
		   2 * P2 * (1-Fl) * (1-ctK*ctK) * ctL +
		   P3 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * sin(2*phi) );

  double effCValue = effCVal()->getVal();
  if (effCValue<0)  std::cout<<"ERROR! NEGATIVE CT EFFICIENCY SPOTTED AT ("<<ctK<<","<<ctL<<","<<phi<<"): "<<effCValue<<std::endl;
  if (effCValue==0) std::cout<<"ERROR! ZERO CT EFFICIENCY SPOTTED AT ("    <<ctK<<","<<ctL<<","<<phi<<"): "<<effCValue<<std::endl;

  double effWValue = effWVal()->getVal();
  if (effWValue<0)  std::cout<<"ERROR! NEGATIVE WT EFFICIENCY SPOTTED AT ("<<ctK<<","<<ctL<<","<<phi<<"): "<<effWValue<<std::endl;
  if (effWValue==0) std::cout<<"ERROR! ZERO WT EFFICIENCY SPOTTED AT ("    <<ctK<<","<<ctL<<","<<phi<<"): "<<effWValue<<std::endl;

  double penalty = penTermVal()->getVal();

  double ret = (9./(32 * 3.14159265) * (effCValue * decCT + mFrac * effWValue * decWT) * penalty);

  if ( ret < 1e-250 ) return 1e-250;

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

Int_t PdfSigAng_Pen::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // use pre-computed integrals for the integration over all the three variables
  // after checking that integrals are in place
  if ( intCPart.size()<1 || intCPart[0]==0 || intWPart.size()<1 || intWPart[0]==0 )
    return 0 ;
  if ( matchArgs(allVars,analVars,ctK,ctL,phi) )
    if ( fullRangeCosT(ctK,rangeName) && fullRangeCosT(ctL,rangeName) && fullRangePhi(phi,rangeName) )
      return 1 ;

  // the lack of analytical integral for the subsets of angular variables does not slow down the fit
  // since only the complete integration is used there
  // if one wants to speed up also the PDF projection for plotting, the other analytical integrals can be computed
  // but it seems a huge effort to me
  return 0 ;

}

Double_t PdfSigAng_Pen::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code>0 && code<2) ;

  // use the pre-computed integrals from histogram
  Double_t retCT =  9./(32*3.14159265) * (
					  0.75*(1-Fl)              * intCPart[0]
					  + Fl                     * intCPart[1]
					  + 0.25*(1-Fl)            * intCPart[2]
					  - Fl                     * intCPart[3]
					  + 0.5*P1*(1-Fl)          * intCPart[4]
					  + 0.5*sqrt(Fl-Fl*Fl)*P4p * intCPart[5]
					  + sqrt(Fl-Fl*Fl)*P5p     * intCPart[6]
					  - sqrt(Fl-Fl*Fl)*P6p     * intCPart[7]
					  + 0.5*sqrt(Fl-Fl*Fl)*P8p * intCPart[8]
					  + 2*(1-Fl)*P2            * intCPart[9]
					  - P3*(1-Fl)              * intCPart[10]
					  );
  
  Double_t retWT =  9./(32*3.14159265) * (
					  0.75*(1-Fl)              * intWPart[0]
					  + Fl                     * intWPart[1]
					  + 0.25*(1-Fl)            * intWPart[2]
					  - Fl                     * intWPart[3]
					  + 0.5*P1*(1-Fl)          * intWPart[4]
					  + 0.5*sqrt(Fl-Fl*Fl)*P4p * intWPart[5]
					  - sqrt(Fl-Fl*Fl)*P5p     * intWPart[6]
					  - sqrt(Fl-Fl*Fl)*P6p     * intWPart[7]
					  - 0.5*sqrt(Fl-Fl*Fl)*P8p * intWPart[8]
					  - 2*(1-Fl)*P2            * intWPart[9]
					  + P3*(1-Fl)              * intWPart[10]
					  );


  if (retCT<=0) {
    if (retCT<0) std::cout<<"ERROR! Negative ct pdf integral, fake value returned"<<std::endl;
    else std::cout<<"ERROR! Null ct pdf integral, fake value returned"<<std::endl;
    return 1e-55;
  }
  if (retWT<=0) {
    if (retWT<0) std::cout<<"ERROR! Negative wt pdf integral, fake value returned"<<std::endl;
    else std::cout<<"ERROR! Null wt pdf integral, fake value returned"<<std::endl;
    return 1e-55;
  }
 
  return (retCT+mFrac*retWT) ;

}
