//#pragma once			/* guard */
#ifndef __POPULATION_INCLUDED__
#define __POPULATION_INCLUDED__

#include "person.h"
#include "parameters.h"
#include "results.h"



/////////////////////////////////////
//DECLARATION POPULATION CLASS
class population{
 public:
  //declarations
  double TBRATEN, TBRATEP, relpop;
  int hiv, hiv1549, hiv18, hiv15;//HIV population
  unsigned int size, size0, size1549, nohh, sizeyoung;
  unsigned int noevs;
  unsigned int noindstates;
  vector< double > ghaz;//population event hazards
  vector< person > v;//the vector of persons
  vector< household > hh;//the vector of households
  vector< unsigned int > popsN, popsP;//number of HIV-/+ people in each state

  double FOI, rFOI, time, hFOI;
  double hhLL;//household log-likelihood
  //annoying counters -- see pop::scheduledevents for defns

  int noh, noa, noaf, noafp, noafn, noc, nnd, noart, noah, noahs, noaa, noap,kidi, adi,slow, norx, nob, noart15;
  int ntbd, ntbdet, ntbsc, ntbdh, ntbdeth, nod, dethhno, dethhno1549, ntbdetprev, ntbretx, catII, iptends;
  int ntbdetsp, ntbdetf, newhhcount, bdcount;//smr+ dets, fast dets, kid infections, background death count
  double av, av1,av2, av3, av4, av5, tauav, prevhh, tbavd, tbav, ARTnum, ARTdenom, shivcount, foibar, hhfoibar, tbcoprev, hivcoprev, tbavaod;
  double avAdeath;
  double propFI, propFIP, propFD, propIH, propDH, propDS, prevHHH, taurx, foik, rfoi, propECF, propViaF;
  int nophhi, nohhart, nohhipt, notecf, nodtecf, cipth;//number household experienced/ hhART /hhIPT  
  int size18, prev18, prevSP18, prevT18, prevTH18, rL;  /* adult population */
  vector< int > prev_rx, prev_rx_cp;                                /* treatment counters */
  int nbc_rx, nint_rx; //* treatment counters
  double iptnos, iptnosP, iptnosA; /* numbers on IPT */

  //vector< vector< double > > hhNnm;//matrix of household numbers by gender
  vector< vector< unsigned int > > hhNnm;//matrix of household numbers by gender
  vector< unsigned int > hhNn;//this is just the householdsize distn..
  unsigned int hhmaxn, hhmaxm, hhmax;//max hh dimensions
  vector< double > f_weight, a_weight;//vectors for detailed balance hh method
  vector< double > cd4hist;	      /* current hist */
  //functions
  int event( unsigned int who,  int etype, parameters& parz, const gsl_rng *r );//event enactor
  int event( unsigned int who, string evtype, map < string, int >& evk, parameters& parz, const gsl_rng *r );
  int build_new_households( parameters& parz, int nob, const gsl_rng *r );
  double update_hazards( ofstream& lgf, parameters& parz );//hazard and foi updater; returns foibar, takes hh foi and returns by sideeffect
  int update_ages( const gsl_rng *r, parameters& parz ); //incrementer of times
  unsigned int setup_weights(ofstream& lgf, parameters& parz, int kappa, map<int,int>& ghk, vector< double >& weights, vector< unsigned int>& where, unsigned int& noe, unsigned int& num, const gsl_rng *r );
  int hhmovement( const gsl_rng *r, ofstream& lgf, int num, parameters& parz, bool record );//move them around between households
  int counts( ofstream& lgf, parameters& parz );//counts popsP etc.
  int scheduledevents( ofstream& lgf, const gsl_rng *r, bool qlog, parameters& parz, map < string, int >& evk );//for things going on in events
  int hhcounts( ofstream& lgf, parameters& parz );//safer way to count up the household states
  vector< double > FOIvec, rFOIvec;//foi across different age cats (resistant)
  vector< double > BINvec;//population sizes across different categories
  int HHpostinit( const gsl_rng *r, ofstream& lgf, parameters& parz );//for post-processing a read population hh data
  int snapshot( ofstream& lgf, vector< vector< double > >& ibmstate );//snapshot of popn
  int initialinfect( const gsl_rng *r, ofstream& lgf, parameters& parz ); /* for initialisation */
  int reinitialize( const gsl_rng *r, ofstream& lgf, parameters& parz  ); /* hmmm... */
  int logstate( ofstream& lgf, int runno, int i );
  int record( ofstream& lgf, results& rez, int runno, int i, parameters& parz ); /* record population in results */
  int demosnapshot( ofstream& lgf, results& rez, int runno, int i, parameters& parz );
  int event_sweep( map < string, int >& swk, map < string, int >& evk, map<int,int>& ghk, const gsl_rng *r, ofstream& lgf, parameters& parz );
  int ibmrun( const gsl_rng *r, ofstream& lgf, parameters& parz, results& rez, int runno); /* the actual IBM */

  /* HIV MAC counters */
  int GLRn, GLRp, iGLR, GLRlast, N15, N15n, N15_0, prev15, h15u200, h15u200an, GLRpo350, GL15, GLnoteli, prevu200, prevo200, prevART, GLpnoART, Noteli;
  double GLRnCD41, GLRnCD42, GLRnCD43, GLRpCD41, GLRpCD42, GLRpCD43; /* GLR or not by CD4 count */
  double prev15n, prev15an3, prev15an2, prev15an1;//some TB indicator counter
  double hpop1, hpop2, hpop3;
  /* HIVMAC */

  //constructor initializer (n.b. needs to be fu with full state mods after file reading...)
 population(int s = 1000, unsigned int h = 300, unsigned int n = 2, unsigned int m = 2, unsigned int ncat = 1): //takes popsize, noevs, nostates
  size(s), nohh(h), noevs(n), noindstates(m), ghaz( n, 0 ),  v( s, person(0,23,n) ), hh( h, household() ), popsN( m, 0 ), popsP( m, 0), FOIvec( ncat, 0 ), BINvec( ncat, 0 ) {}

private:
  //nowt

};// end population class

#endif
