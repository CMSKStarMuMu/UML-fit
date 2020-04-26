/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "BoundCheck.h" 

ClassImp(BoundCheck) 

BoundCheck::BoundCheck(const char *name, const char *title, 
		       RooAbsReal& _P1,
		       RooAbsReal& _P2,
		       RooAbsReal& _P3,
		       RooAbsReal& _P4p,
		       RooAbsReal& _P5p,
		       RooAbsReal& _P6p,
		       RooAbsReal& _P8p,
		       bool _useCTL4,
		       double _useCTL15,
		       bool _verbose) :
  RooAbsReal(name,title), 
  P1("P1","P1",this,_P1),
  P2("P2","P2",this,_P2),
  P3("P3","P3",this,_P3),
  P4p("P4p","P4p",this,_P4p),
  P5p("P5p","P5p",this,_P5p),
  P6p("P6p","P6p",this,_P6p),
  P8p("P8p","P8p",this,_P8p)
{

  useCTL4 = _useCTL4;
  useCTL15 = _useCTL15;

  verbose = _verbose;

}


BoundCheck::BoundCheck(const BoundCheck& other, const char* name) :  
  RooAbsReal(other,name), 
  P1("P1",this,other.P1),
  P2("P2",this,other.P2),
  P3("P3",this,other.P3),
  P4p("P4p",this,other.P4p),
  P5p("P5p",this,other.P5p),
  P6p("P6p",this,other.P6p),
  P8p("P8p",this,other.P8p)
{

  useCTL4 = other.useCTL4;
  useCTL15 = other.useCTL15;

  verbose = other.verbose;

}



Double_t BoundCheck::evaluate() const 
{

  double ctL4phi1 = P4p*P4p + P5p*P5p + P6p*P6p + P8p*P8p - 2 + 2*fabs( 2*P2 - P4p*P5p +P6p*P8p );

  if (verbose && ctL4phi1>0) std::cout<<"[OUT] 4 "<<ctL4phi1<<std::endl;

  if (useCTL4 && ctL4phi1>0) return 1;

  if (useCTL15<0) return 0;

  double a0 = 1 - P1*P1 - P6p*P6p*(1+P1) - P8p*P8p*(1-P1) - 4*P2*P2 - 4*P2*P6p*P8p; 
  double a4 = 1 - P1*P1 - P4p*P4p*(1+P1) - P5p*P5p*(1-P1) - 4*P2*P2 + 4*P2*P4p*P5p; 

  double a1 = 4*P3*P8p*P8p - 4*P3*P6p*P6p - 8*P1*P3 + 2*P5p*P6p*(1+P1) - 2*P4p*P8p*(1-P1) - 4*P2*P4p*P6p + 4*P2*P5p*P8p;
  double a3 = 4*P3*P4p*P4p - 4*P3*P5p*P5p + 8*P1*P3 + 2*P5p*P6p*(1-P1) - 2*P4p*P8p*(1+P1) - 4*P2*P4p*P6p + 4*P2*P5p*P8p;

  double a2 = 2 + 2*P1*P1 - 8*P2*P2 - 16*P3*P3 - (P4p*P4p+P6p*P6p)*(1-P1) - (P5p*P5p+P8p*P8p)*(1+P1) + 4*P2*P4p*P5p - 4*P2*P6p*P8p + 8*P3*P4p*P8p + 8*P3*P5p*P6p;

  double b0 = P8p*P8p - 1 - P1 + 2*P2 + P6p*P8p; 
  double b2 = P4p*P4p - 1 + P1 + 2*P2 - P4p*P5p; 
  double b1 = P4p*P8p - 2*P3 + 0.5 * ( P4p*P6p - P5p*P8p );
  
  double c0 = P8p*P8p - 1 - P1 - 2*P2 - P6p*P8p; 
  double c2 = P4p*P4p - 1 + P1 - 2*P2 + P4p*P5p; 
  double c1 = P4p*P8p - 2*P3 - 0.5 * ( P4p*P6p - P5p*P8p );
  
  int nSteps = 500;
  double phi, sin2, sincos, cos2;
  double ctL1, ctL5p, ctL5m;
  for (int step = 0; step<nSteps; ++step) {
    phi = 3.14159 * step / nSteps;
    if (useCTL15>0) phi=useCTL15;

    sin2 = sin(phi)*sin(phi);
    sincos = sin(phi)*cos(phi);
    cos2 = cos(phi)*cos(phi);

    ctL5p = b0*sin2 + b1*sincos + b2*cos2;
    if (ctL5p >= 0) continue;

    ctL5m = c0*sin2 + c1*sincos + c2*cos2;
    if (ctL5m >= 0) continue;

    ctL1 = a0*sin2*sin2 + a1*sin2*sincos + a2*sin2*cos2 + a3*sincos*cos2 + a4*cos2*cos2;
    if (ctL1 >= 0) continue;

    if (verbose) std::cout<<"[OUT] 1 "<<ctL1<<" 5p "<<ctL5p<<" 5m "<<ctL5m<<" phi "<<phi<<std::endl;

    return 1;
    if (useCTL15>0) break;

  }
  
  return 0;

}