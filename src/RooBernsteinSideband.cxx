/***************************************************************************** 
 * Project: RooBernsteinSideband                                             * 
 *                                                                           * 
 * P.Dini fecit, Anno Domini MMXVIII                                         *
 *                                                                           * 
 * "Est modus in rebus."                                                     * 
 *                                                                           * 
 * Class to describe 3D angular Sideband in B0->K*mumu Analysis              * 
 *                                                                           * 
 * Ipse mutavit, Anno Domini MMXXIV                                          * 
 *                                                                           * 
 * "Et de hoc satis"                                                         * 
 *                                                                           * 
 * Now with a trivial circularity condition added                            * 
 *                                                                           * 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooBernsteinSideband.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h"
#include "RooRealProxy.h" 
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include <math.h> 
#include "TMath.h"  
#include <algorithm>  

ClassImp(RooBernsteinSideband); 
   
  
   RooBernsteinSideband::RooBernsteinSideband(const char *name, const char *title, 
                        RooAbsReal& x, RooAbsReal& y, RooAbsReal& z, 
			const RooArgList& coefList,
			int maxDegree1, int maxDegree2, int maxDegree3
			) :
   RooAbsPdf(name,title), 
   _x("_x","_x",this,x),
   _y("_y","_y",this,y),
   _z("_z","_z",this,z), 
   _coefList("_coefList","_coefList",this),
   _maxDegree1(maxDegree1),
   _maxDegree2(maxDegree2),
   _maxDegree3(maxDegree3)
{ 
    TIterator* coefIter = coefList.createIterator() ;
     RooAbsArg* coef ;
     while((coef = (RooAbsArg*)coefIter->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef)) {
    	 std::cout << "RooBernsteinSideband::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList.add(*coef) ;
//   	 std::cout << "RooBernsteinSideband::ctor(" << GetName() << ") coSidebandcient " << coef->GetName()<<std::endl ;
      }
      int icc=0;
      for(int i = 0; i <= maxDegree1 ; ++i) {
     	for(int j = 0; j <= maxDegree2 ; ++j) {
     	 for(int k = 0; k <  maxDegree3 ; ++k) {
 
     	 printf("RooBernsteinSideband: Par(%d)=%f cosL=%d cosK=%d phi=%d\n",icc,((RooAbsReal&) _coefList[icc]).getVal(), i,j,k);
     	 icc++;
 
     	 }
     	}
      }
     delete coefIter;
     _maxDegreeV = std::max(maxDegree3,std::max(maxDegree1,maxDegree2));
     printf("RooBernstein: Max(Numbers of degree) = %d\n",_maxDegreeV);
     if ( _maxDegreeV>25) {
       printf("RooBernsteinSideband error: Max(Numbers of degree) > 25 = %d\n",_maxDegreeV);
       assert(0) ;
     }
 }   

RooBernsteinSideband::RooBernsteinSideband(const RooBernsteinSideband& other, const char* name) :  
   RooAbsPdf(other,name), 
   _x("_x",this,other._x),
   _y("_y",this,other._y),
   _z("_z",this,other._z), 
   _coefList("_coefList",this,other._coefList),
   _maxDegree1(other._maxDegree1),
   _maxDegree2(other._maxDegree2),
   _maxDegree3(other._maxDegree3),
   _maxDegreeV(other._maxDegreeV)
{ 
} 
//
Double_t RooBernsteinSideband::evaluate() const 
{ 
// //   
        
//       double __z=fabs(_z);
	
       double xmin = _x.min();
       double xdif = _x.max()-xmin;
       double x    =(_x-xmin)/xdif;
       double ymin = _y.min();
       double ydif = _y.max()-ymin;
       double y    =(_y-ymin)/ydif;
       double zmin = _z.min();
       double zdif = _z.max()-zmin;
       double z    =(_z-zmin)/zdif;
       
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       
       
       sx[0]=1.0;
       sy[0]=1.0;
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegreeV ; ++i){
//         sx[i]= sx[i-1]*(1.-x);
//         sy[i]= sy[i-1]*(1.-y);
//         sz[i]= sz[i-1]*(1.-z);
        if(i<=_maxDegree1) sx[i]= sx[i-1]*(1.-x);
        if(i<=_maxDegree2) sy[i]= sy[i-1]*(1.-y);
        if(i<=_maxDegree3) sz[i]= sz[i-1]*(1.-z);
       }
       
//        fptype tx = 1.0;
//        fptype ty = 1.0;
//        fptype tz = 1.0;
//        for(  int i = 0; i <= _maxDegreeV ; ++i){
//         int ii = _maxDegreeV-i;
//         if(ii<=_maxDegree1) {sx[ii]*= tx*device_coeffbinomial(_maxDegree1,ii);tx*=x;}
//         if(ii<=_maxDegree2) {sy[ii]*= ty*device_coeffbinomial(_maxDegree2,ii);ty*=y;}
//         if(ii<=_maxDegree3) {sz[ii]*= tz*device_coeffbinomial(_maxDegree3,ii);tz*=z;}
// 	
//        }
   
// 
       double bernknvalx = 0.;
       double bernknvaly = 0.;
       double bernknvalz = 0.;
       
//        double sx_ini = pow(1.0-x,_maxDegree1);
//        double sy_ini = pow(1.0-y,_maxDegree2);
//        double sz_ini = pow(1.0-z,_maxDegree3);
//       long double sx_div = 1.0/(1.0-x);
//       long double sy_div = 1.0/(1.0-y);
//       long double sz_div = 1.0/(1.0-z);
       
       int ipar =0;
       int ipa0 =0;
       double func =0.0;
       double intg_1 =0.0;
       
       double tx = 1.;
//       long double sx = sx_ini;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         //if(((RooAbsReal&) _coefList[ipar]).getVal()!=0.000000000) 
	 bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
       
         double ty = 1.;
//         long double sy =sy_ini;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	 //if(((RooAbsReal&) _coefList[ipar]).getVal()!=0.000000000)  
	 bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
	 
          double tz = 1.;
//          long double sz = sz_ini;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
//	   ipar++;
//	   if((((RooAbsReal&) _coefList[ipar-1]).getVal())==0.0) continue;
// 	   func += ((RooAbsReal&) _coefList[ipar-1]).getVal()*sx[_maxDegree1-i]*sy[_maxDegree2-j]*sz[_maxDegree3-k];
//           intg_1 +=((RooAbsReal&) _coefList[ipar-1]).getVal();
//	   if(fabs(((RooAbsReal&) _coefList[ipar]).getVal())>0.0){
	    if (k==0) ipa0=ipar;
 	    bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
 	    if(k==_maxDegree3){
 	     func   += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	     intg_1 += ((RooAbsReal&) _coefList[ipa0]).getVal();
	    }else{
 	     func   += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	     intg_1 += ((RooAbsReal&) _coefList[ipar]).getVal();
	     ipar++;
	    } 
// 	   } 
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
	   tz *= z;
	  }
	   ty *= y;
         }
	   tx *= x;
       }
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//      if(func<1.E-30)  std::cout<<"evaluate = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;
//
// evaluate non deve essere normalizzato!
//
//      intg_1 = (1.0+_maxDegree1)*(1.0+_maxDegree2)*(1.0+_maxDegree3)/intg_1;
//      intg_1=intg_1/(xdif*ydif*zdif);
//      func = func*intg_1;
 //     if(func<1.E-30)  std::cout<<"WARNING! eval = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;;
//      if(func<1.E-30)  return 1.E-30;
      return  func;
} 

fptype RooBernsteinSideband::device_coeffbinomial  (fptype enne, fptype kappa) const {
        fptype factor=1.;
        for(fptype i = 1; i <=kappa; ++i) {
          factor *= (enne+1-i)/i; 
        }	 
        if (factor<=0 ){
	 printf("Error in RooBernsteinSideband coeffbinomial=> factor = %f enne=%f kappa=%f\n",factor,enne,kappa);
         return 0;
	} 
        return factor;
} 

fptype  RooBernsteinSideband::device_bernsteinkn_func  (fptype x, fptype enne, fptype kappa) const{
   return TMath::Binomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);
//   return RooBernsteinSideband::device_coeffbinomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);
}

 fptype RooBernsteinSideband::evaluateInt(fptype xBinw,fptype yBinw,fptype zBinw) const {

// evaluate (integral of the function over the bin) / (bin volume) 
 
    int maxDegree1      = _maxDegree1;
    int maxDegree2      = _maxDegree2;
    int maxDegree3      = _maxDegree3;
    fptype x    = _x; 
    fptype y    = _y; 
    fptype z    = _z; 

//     fptype xBinw = 2.0/25.;
//     fptype yBinw = 2.0/25.;
//     fptype zBinw = 2.0*TMath::Pi()/25.;


    fptype xmin = _x.min();
    fptype xdif = _x.max() - xmin ;

    fptype xLeft  = ((x-xBinw/2.)-xmin)/xdif;
    fptype xRight = ((x+xBinw/2.)-xmin)/xdif;
//    x=(x-xmin)/xdif;
    fptype ymin = _y.min();
    fptype ydif = _y.max()-ymin;
    fptype yLeft  = ((y-yBinw/2.)-ymin)/ydif;
    fptype yRight = ((y+yBinw/2.)-ymin)/ydif;
//    y=(y-ymin)/ydif;

    fptype zmin = _z.min();
    fptype zdif = _z.max()-zmin;
    fptype zLeft  = ((z-zBinw/2.)-zmin)/zdif;
    fptype zRight = ((z+zBinw/2.)-zmin)/zdif;
//    z=(z-zmin)/zdif;
    
       int ipar = 0 ;
       int ipa0 =0;
       fptype ret   =0;
       fptype intg_1 =0.;
       for(int i = 0; i <= maxDegree1 ; ++i) {
         for(int j = 0; j <= maxDegree2 ; ++j) {
//	  std::cout<<"func = par["<<ipar<<"]*x^"<<kk<<"*y^"<<jj<<std::endl;
          for(int k = 0; k <= maxDegree3 ; ++k) {
//
	   if (k==0) ipa0=ipar;
           fptype bernknintgbinx = device_SidebandBernsteinkn_intgBin(xLeft,xRight,maxDegree1,i);
           fptype bernknintgbiny = device_SidebandBernsteinkn_intgBin(yLeft,yRight,maxDegree2,j);
           fptype bernknintgbinz = device_SidebandBernsteinkn_intgBin(zLeft,zRight,maxDegree3,k);
 	   if(k==maxDegree3){
            ret   +=((RooAbsReal&) _coefList[ipa0]).getVal()*bernknintgbinx*bernknintgbiny*bernknintgbinz;
            intg_1+=((RooAbsReal&) _coefList[ipa0]).getVal();
 	   }else{
            ret   +=((RooAbsReal&) _coefList[ipar]).getVal()*bernknintgbinx*bernknintgbiny*bernknintgbinz;
            intg_1+=((RooAbsReal&) _coefList[ipar]).getVal();
	    ipar++;
	   } 
	  }
         }
       }
     //okkio!!! if(ret<1.E-30) ret = 1.E-30;
    intg_1 = (1.0+maxDegree1)*(1.0+maxDegree2)*(1.0+maxDegree3)/intg_1;
//    ret=ret/(xBinw*yBinw*zBinw);
//P.S. Il termmine xdif*ydif*zdif si elide sopra e sotto
      ret = ret*intg_1;
//    std::cout<<"ret = "<<ret<<std::endl;
   return ret;

 }
fptype  RooBernsteinSideband::device_SidebandBernsteinkn_intgBin( fptype xLeft, fptype xRight, fptype enne, fptype kappa) const{

      if(xLeft==0.0&&xRight==1.0) return 1.0/(enne+1.0);
 
      fptype integbernkn = 0.0;
      fptype ifactni = 0.0;
      fptype ifactik = 0.0;
      
      
      fptype powxL = pow(xLeft ,kappa+1)   ;
      fptype powxR = pow(xRight,kappa+1) ;
       for(fptype i = kappa; i <=enne ; ++i) {
// n!/(i!(n-i)!)

//        ifactni =  device_coeffbinomial(enne,i);
        ifactni =  TMath::Binomial(enne,i);
// i!/(k!(i-k)!)
 
        ifactik =  TMath::Binomial(i,kappa);
//        ifactik =  device_coeffbinomial(i,kappa);
        integbernkn += ifactni*ifactik*pow(-1.0,i-kappa)*(powxR-powxL)/(i+1);
	powxL*=xLeft ;
	powxR*=xRight;
//        integbernkn += ifactni*ifactik*pow(-1.0,i-kappa)*(pow(xRight, i+1)-pow(xLeft,i+1))/(i+1);
       }

       if (integbernkn<=0.0 ){
          printf(" Error in SidebandBernsteinkn_intgbin xLeft=%f xRight=%f kappa=%f enne=%f integral = %5.15f\n",xLeft,xRight,kappa,enne,integbernkn);
        integbernkn=1.E-30;
       }
       return integbernkn;
}
 
Int_t RooBernsteinSideband::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if       (matchArgs(allVars, analVars, _x,_y,_z)){
//    std::cout<<"getAnalyticalIntegral==1: Analytic integral over x y z"<<std::endl;
   return 1;
  }else if (matchArgs(allVars, analVars, _y,_z)){
//    std::cout<<"getAnalyticalIntegral==2: Analytic integral over y z"<<std::endl;
   return 2;
  }else if (matchArgs(allVars, analVars, _x,_z)){
//    std::cout<<"getAnalyticalIntegral==3: Analytic integral over x z"<<std::endl;
   return 3;
  }else if (matchArgs(allVars, analVars, _x,_y)){
//    std::cout<<"getAnalyticalIntegral==4: Analytic integral over x y"<<std::endl;
   return 4;
  }else if (matchArgs(allVars, analVars, _x   )){
//    std::cout<<"getAnalyticalIntegral==5: Analytic integral over x  "<<std::endl;
   return 5;
  }else if (matchArgs(allVars, analVars, _y   )){
//    std::cout<<"getAnalyticalIntegral==6: Analytic integral over y  "<<std::endl;
   return 6;
  }else if (matchArgs(allVars, analVars, _z   )){
//    std::cout<<"getAnalyticalIntegral==7: Analytic integral over z  "<<std::endl;
   return 7;
  }else{  
//    std::cout << "Error in RooBernsteinSideband::analyticalIntegral: set getAnalyticalIntegral==0" << std::endl;
//    std::cout << "==> analVars = " << analVars << std::endl ;
   return 0;
 } 
}
Double_t RooBernsteinSideband::analyticalIntegral(Int_t code, const char* rangeName) const
{
    if(code==1){ // X Y Z INTEGRATED

//        printf("integrale xyz\n");
//       fptype xBinw=_x.max(rangeName)-_x.min(rangeName);
//       fptype yBinw=_y.max(rangeName)-_y.min(rangeName);
//       fptype zBinw=_z.max(rangeName)-_z.min(rangeName);

      fptype xmin = _x.min();
      fptype xdif = _x.max() - xmin ;

      fptype xLeft  = (_x.min(rangeName)-xmin)/xdif;
      fptype xRight = (_x.max(rangeName)-xmin)/xdif;
//      x=(x-xmin)/xdif;
      fptype ymin = _y.min();
      fptype ydif = _y.max()-ymin;
      fptype yLeft  = (_y.min(rangeName)-ymin)/ydif;
      fptype yRight = (_y.max(rangeName)-ymin)/ydif;
//      y=(y-ymin)/ydif;
      fptype zmin = _z.min();
      fptype zdif = _z.max()-zmin;
      fptype zLeft  = (_z.min(rangeName)-zmin)/zdif;
      fptype zRight = (_z.max(rangeName)-zmin)/zdif;
//      z=(z-zmin)/zdif;
 
      fptype bernknintgbinx[9] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      fptype bernknintgbiny[9] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      fptype bernknintgbinz[9] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
 
      fptype powdx[9]	       =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      fptype powdy[9]	       =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      fptype powdz[9]	       =  {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      
      int enne = _maxDegreeV;
      
      if ( enne>8) {
       printf("RooBernsteinSideband::analyticalIntegral: Max(Numbers of degree) > 8 = %d\n",enne);
       return 0.0;
      }

      fptype x_Left  = xLeft ;
      fptype y_Left  = yLeft ;
      fptype z_Left  = zLeft ;
      fptype x_Right = xRight;
      fptype y_Right = yRight;
      fptype z_Right = zRight;
      for(int k = 0; k <=enne;++k) {
      	powdx[k] = x_Right - x_Left;
      	powdy[k] = y_Right - y_Left;
        powdz[k] = z_Right - z_Left;
	x_Left	*= xLeft ;
	y_Left	*= yLeft ;
	z_Left	*= zLeft ;
	x_Right *= xRight;
	y_Right *= yRight;
	z_Right *= zRight;
      }
      for(int kappa = 0; kappa <=enne; ++kappa) {
      int powmenuno = 1;
//
       for(int i = kappa; i <=enne ; ++i) {
         fptype comb = device_coeffbinomial(i,kappa)*powmenuno/(i+1);
	 powmenuno *=-1;
         if (i<=_maxDegree1) {bernknintgbinx[kappa] += device_coeffbinomial(_maxDegree1,i)*comb*powdx[i];};
         if (i<=_maxDegree2) {bernknintgbiny[kappa] += device_coeffbinomial(_maxDegree2,i)*comb*powdy[i];};
         if (i<=_maxDegree3) {bernknintgbinz[kappa] += device_coeffbinomial(_maxDegree3,i)*comb*powdz[i];};
       }
      }
         int ipar = 0 ;
         int ipa0 =0;
         fptype ret   =0;
         for(int i = 0; i <= _maxDegree1 ; ++i) {
           for(int j = 0; j <= _maxDegree2 ; ++j) {
//          std::cout<<"func = par["<<ipar<<"]*x^"<<kk<<"*y^"<<jj<<std::endl;
            for(int k = 0; k <= _maxDegree3 ; ++k) {

	     if (k==0) ipa0=ipar;
//              fptype bernknintgbinx = device_SidebandBernsteinkn_intgBin(xLeft,xRight,_maxDegree1,i);
//              fptype bernknintgbiny = device_SidebandBernsteinkn_intgBin(yLeft,yRight,_maxDegree2,j);
//              fptype bernknintgbinz = device_SidebandBernsteinkn_intgBin(zLeft,zRight,_maxDegree3,k);
//           fptype bernknintgbinx = 1./(maxDegree1+1.);
//           fptype bernknintgbiny = 1./(maxDegree2+1.);
//           fptype bernknintgbinz = 1./(maxDegree3+1.);
//             ret   +=((RooAbsReal&) _coefList[ipar]).getVal()*bernknintgbinx*bernknintgbiny*bernknintgbinz;
	     if(k==_maxDegree3){
              ret   +=((RooAbsReal&) _coefList[ipa0]).getVal()*bernknintgbinx[i]*bernknintgbiny[j]*bernknintgbinz[k];
 	     }else{
              ret   +=((RooAbsReal&) _coefList[ipar]).getVal()*bernknintgbinx[i]*bernknintgbiny[j]*bernknintgbinz[k];
              ipar++;
	     } 
            }
 
 
           }
         }
      ret=ret*(xdif*ydif*zdif);
       //okkio!!! if(ret<1.E-30) ret = 1.E-30;
//              std::cout<<"int xyz = "<<ret<<std::endl;
      return ret;
//     
   } 
   if(code>1){
//    
    fptype xmin = _x.min();
    fptype xdif = _x.max()-xmin;
    fptype x=(_x-xmin)/xdif;
    fptype ymin = _y.min();
    fptype ydif = _y.max()-ymin;
    fptype y=(_y-ymin)/ydif;
    fptype zmin = _z.min();
    fptype zdif = _z.max()-zmin;
    fptype z=(_z-zmin)/zdif;
//    
    fptype xLeft  = (_x.min(rangeName)-xmin)/xdif;
    fptype xRight = (_x.max(rangeName)-xmin)/xdif;
    fptype yLeft  = (_y.min(rangeName)-ymin)/ydif;
    fptype yRight = (_y.max(rangeName)-ymin)/ydif;
    fptype zLeft  = (_z.min(rangeName)-zmin)/zdif;
    fptype zRight = (_z.max(rangeName)-zmin)/zdif;
  
    
//    std::cout<<"CODE = "<<code<<std::endl;
//
    if (code==2) { // Y Z INTEGRATED
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sx[0]=1.0;
       for( int i = 1; i <= _maxDegree1 ; ++i){
        sx[i]= sx[i-1]*(1.-x);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       fptype tx = 1.;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         for(int j = 0; j <= _maxDegree2 ; ++j) {
          for(int k = 0; k <= _maxDegree3 ; ++k) {
//           fptype bernknvalx =  device_bernsteinkn_func(x,_maxDegree1,i);
	   if (k==0) ipa0=ipar;
           fptype bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
 	   fptype bernknvaly =  device_SidebandBernsteinkn_intgBin(yLeft,yRight,_maxDegree2,j);
	   fptype bernknvalz =  device_SidebandBernsteinkn_intgBin(zLeft,zRight,_maxDegree3,k);
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   }
	  }
         }
	 tx*=x;
       }
//       std::cout<<"int yz = "<<func/xdif<<std::endl;
       return  func*ydif*zdif;
    }
    if (code==3) { // X Z INTEGRATED
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sy[0]=1.0;
       for( int i = 1; i <= _maxDegree2 ; ++i){
        sy[i]= sy[i-1]*(1.-y);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         fptype ty = 1.;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   if (k==0) ipa0=ipar;
 	   fptype bernknvalx =  device_SidebandBernsteinkn_intgBin(xLeft,xRight,_maxDegree1,i);
//           fptype bernknvaly =  device_bernsteinkn_func(y,_maxDegree2,j);
           fptype bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
	   fptype bernknvalz =  device_SidebandBernsteinkn_intgBin(zLeft,zRight,_maxDegree3,k);
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   } 
	  }
	  ty*=y;
         }
       }
       return  func*xdif*zdif;
    }
    if (code==4) {// X Y INTEGRATED
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegree3 ; ++i){
        sz[i]= sz[i-1]*(1.-z);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	  fptype tz=1.;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   if (k==0) ipa0=ipar;
 	   fptype bernknvalx =  device_SidebandBernsteinkn_intgBin(xLeft,xRight,_maxDegree1,i);
 	   fptype bernknvaly =  device_SidebandBernsteinkn_intgBin(yLeft,yRight,_maxDegree2,j);
//           fptype bernknvalz =  device_bernsteinkn_func(z,_maxDegree3,k);
           fptype bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   } 
	   tz*=z;
	  }
         }
       }
       return  func*xdif*ydif;
    }
    if (code==5) {// X INTEGRATED
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sy[0]=1.0;
       for( int i = 1; i <= _maxDegree2 ; ++i){
        sy[i]= sy[i-1]*(1.-y);
       }
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegree3 ; ++i){
        sz[i]= sz[i-1]*(1.-z);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
	 fptype ty=1.;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	  fptype tz=1.;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   if (k==0) ipa0=ipar;
 	   fptype bernknvalx =  device_SidebandBernsteinkn_intgBin(xLeft,xRight,_maxDegree1,i);
           fptype bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
           fptype bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   } 
	   tz*=z;
	  }
	  ty*=y;
         }
       }
       return  func*xdif;
    }
    if (code==6) {// Y INTEGRATED
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sx[0]=1.0;
       for( int i = 1; i <= _maxDegree1 ; ++i){
        sx[i]= sx[i-1]*(1.-x);
       }
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegree3 ; ++i){
        sz[i]= sz[i-1]*(1.-z);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       fptype tx=1.;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	  fptype tz=1.;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   if (k==0) ipa0=ipar;
           fptype bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
 	   fptype bernknvaly =  device_SidebandBernsteinkn_intgBin(yLeft,yRight,_maxDegree2,j);
           fptype bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   } 
	   tz*=z;
	  }
         }
	 tx*=x;
       }
       return  func*ydif;
    }
    if (code==7) {// Z INTEGRATED
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       sx[0]=1.0;
       for( int i = 1; i <= _maxDegree1 ; ++i){
        sx[i]= sx[i-1]*(1.-x);
       }
       sy[0]=1.0;
       for( int i = 1; i <= _maxDegree2 ; ++i){
        sy[i]= sy[i-1]*(1.-y);
       }
       int ipar =0;
       int ipa0 =0;
       fptype func =0.0;
       fptype tx=1.;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
	 fptype ty=1.;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   if (k==0) ipa0=ipar;
           fptype bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
           fptype bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
 	   fptype bernknvalz =  device_SidebandBernsteinkn_intgBin(zLeft,zRight,_maxDegree3,k);
	   if(k==_maxDegree3){
	    func += ((RooAbsReal&) _coefList[ipa0]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	   }else{
	    func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) coefList[ipar]).getVal()<<std::endl;
	    ipar++;
	   }
	  }
	  ty*=y;
         }
	 tx*=x;
       }
       return  func*zdif;
    }
   } 
    return 0.;

}
