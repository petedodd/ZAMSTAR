//#pragma once			/* guard */
#ifndef __RESULTS_INCLUDED__
#define __RESULTS_INCLUDED__

#include "parameters.h"
//#include "population.h"

//DECLARATION RESULTS CLASS
class results{

 public:
  bool light;			/* record lots? */
  int noruns;//the number of runs
  int notime;//number of times
  double Dm, ARIm, Im, Dsqm, ARIsqm, Isqm, Itmp, aritmp;
  //double dt, starttime;
  vector< double > time;
  vector< vector< double > >  ibmN, ibmS, ibmL, ibmrL, ibmLP,ibmLN,ibmD,ibmD2,ibmT,ibmT2, ibmHIV, ibmHIV2, ibm1549, ibmu12;
  vector< vector< double > >  ibmFI,ibmFIP, ibmFD,  ibmDS, ibmIH, ibmDH, ibmSP;
  vector< vector< double > >  ibmFOI, ibmARI, ibmARIA, TBdeaths, HIVdeaths, TBHIVdeaths, ibmFHH, ibmPHH, ibmPHHH, ibmPhh, ibmPHhh;;
  vector< vector< double >  > propThiv, tbnotes, ibmHHPI, ibmHHipt, ibmHHart, ibmpECF;
  vector< vector< double >  >  ibmvECF, ibmARTp, ibmARTe, ibmPRec, ibmSHIVcount;
  vector< vector< double >  > Inc, ibmstate;//incidence, matrix for popn SS
  vector< double > EP1, EP2, EP3;	    /* each run's measure */
  vector< vector< double >  > mortz, cd4hist;	  /* deaths rates of types */
  vector< vector< double >  > IRR, IRRa, strIn, strIp, strIa, strIpna, strIf, strIpf, strInf;    /* IRR tbhiv and art */
  vector< vector< double >  > IPTprev, IPTprevP, IPTprevA, IPTprevN, IPTch;
  /* for dss */
  vector< vector<double> > P;
  unsigned int Tmax, Tmin, Mmax, Mmin, Wmax, Wmin;
  int hh0tot, hhtot;

  //function
   //  int record( ofstream& lgf, population& pop, int runno, int i, parameters& parz );
  int analyze( ofstream& lgf, parameters& parz );  
  int write( ofstream& lgf, parameters& parz, time_t& writetime );
  int writedss( ofstream& lgf, int i, parameters& parz );
  //int demosnapshot( ofstream& lgf, population& pop, int runno, int i, parameters& parz );

  //constructor initializer: noruns to record, real noruns
 results( int n = 1, int n2 = 1, int m = 1000, bool L = false):
  light(L), noruns(n), notime(m), Dm(0), ARIm(0), Im(0), Dsqm(0), ARIsqm(0), Isqm(0),Itmp(0), aritmp(0),
    time( m, 0 ),ibmN( n, vector< double >(m,0) ),ibmS( n, vector< double >(m,0) ), ibmL( n, vector< double >(m,0) ),ibmrL( n, vector< double >(m,0) ), ibmLP( n, vector< double >(m,0) ),ibmLN( n, vector< double >(m,0) ), ibmD( n, vector< double >(m,0) ),ibmD2( n, vector< double >(m,0) ), ibmT( n, vector< double >(m,0) ),ibmT2( n, vector< double >(m,0) ), ibmHIV( n, vector< double >(m,0) ), 
    ibmHIV2( n, vector< double >(m,0) ),ibm1549( n, vector< double >(m,0) ),ibmu12( n, vector< double >(m,0) ),
    ibmFI( n, vector< double >(m,0) ),ibmFIP( n, vector< double >(m,0) ), ibmFD( n, vector< double >(m,0) ),  ibmDS( n, vector< double >(m,0) ), 
    ibmIH( n, vector< double >(m,0) ), ibmDH( n, vector< double >(m,0) ), ibmSP( n, vector< double >(m,0) ),
    ibmFOI( n, vector< double >(m,0) ), ibmARI( n, vector< double >(m,0) ),ibmARIA( n, vector< double >(m,0) ), TBdeaths( n, vector< double >(m,0) ), HIVdeaths( n, vector< double >(m,0) ), TBHIVdeaths( n, vector< double >(m,0) ),
    ibmFHH( n, vector< double >(m,0) ),ibmPHH( n, vector< double >(m,0) ), ibmPHHH( n, vector< double >(m,0) ),ibmPhh( n, vector< double >(m,0) ), ibmPHhh( n, vector< double >(m,0) ), 
    propThiv( n, vector< double >(m,0) ), tbnotes( n, vector< double >(m,0) ), ibmHHPI( n, vector< double >(m,0) ),   
    ibmHHipt( n, vector< double >(m,0) ), ibmHHart( n, vector< double >(m,0) ),
    ibmpECF( n, vector< double >(m,0) ), ibmvECF( n, vector< double >(m,0) ), ibmARTp( n, vector< double >(m,0) ),ibmARTe( n, vector< double >(m,0) ), ibmPRec( n, vector< double >(m,0) ),ibmSHIVcount( n, vector< double >(m,0) ), Inc( n, vector< double >(m,0) ),
    EP1(n2,0), EP2(n2,0), EP3(n2,0), mortz( n2, vector< double >(4,0) ), cd4hist( m, vector< double >(20,0) ), IRR( n, vector< double >(m,0) ), IRRa( n, vector< double >(m,0) ), strIn( n, vector< double >(m,0) ), strIp( n, vector< double >(m,0) ), strIa( n, vector< double >(m,0) ), strIpna( n, vector< double >(m,0) ), strIf( n, vector< double >(m,0) ),strIpf( n, vector< double >(m,0) ),strInf( n, vector< double >(m,0) ), IPTprev( n, vector< double >(m,0) ), IPTprevP( n, vector< double >(m,0) ), IPTprevA( n, vector< double >(m,0) ),IPTprevN( n, vector< double >(m,0) ), IPTch( n, vector< double >(m,0) ), P( 10, vector< double > (10 , 0 ) ) {}
  
 private:
  //nowt
 
};//end results class
/* the 20 in cd4 hist is where it's set */

#endif
