/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef PDFSIGMASS
#define PDFSIGMASS

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooCategoryProxy.h"
#include "RooSetProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooFit.h"
#include "Riostream.h"
#include "RooObjCacheManager.h"
#include "PdfCBShape.h"

 
class PdfSigMass : public RooAbsPdf {
 protected:


  RooRealProxy m ;
  
  RooRealProxy mean_rt     ;
  RooRealProxy sigma_rt1   ;
  RooRealProxy sigma_rt2   ;
  RooRealProxy alpha_rt1   ;
  RooRealProxy alpha_rt2   ;
  RooRealProxy n_rt1       ;
  RooRealProxy n_rt2       ;
  RooRealProxy f1rt        ;

  RooRealProxy mean_wt     ;
  RooRealProxy sigma_wt1   ;
  RooRealProxy alpha_wt1   ;
  RooRealProxy alpha_wt2   ;
  RooRealProxy n_wt1       ;
  RooRealProxy n_wt2       ;
  
  RooRealProxy mFrac ;

  RooRealProxy rtMassTerm;
  RooRealProxy wtMassTerm;

  Double_t evaluate() const ;

 public:
  PdfSigMass() {} ; 
  PdfSigMass(const char *name, const char *title,
	    RooAbsReal& _m,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac    ,
	    RooAbsReal& _rtMassTerm,
	    RooAbsReal& _wtMassTerm
	    );

  PdfSigMass(const char *name, const char *title,
	    RooAbsReal& _m,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _sigma_rt2,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
	    RooAbsReal& _f1rt     ,  
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac ,
	    RooAbsReal& _rtMassTerm,
	    RooAbsReal& _wtMassTerm
  	    );

  PdfSigMass(const PdfSigMass& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new PdfSigMass(*this,newname); }
  inline virtual ~PdfSigMass() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  ClassDef(PdfSigMass,1) // PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events
    };
 
#endif
