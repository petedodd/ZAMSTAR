//#pragma once			/* guard */
#ifndef __PERSON_INCLUDED__
#define __PERSON_INCLUDED__

#include "parameters.h"

//DECLARATION PERSON CLASS
class person{

public:
  //declarations

  double tbrate; 
  unsigned int gender;//0 = W, 1 = M.
  unsigned int hiv, Lhiv, hhhiv, dart, bbs, shiv, dshiv, hhtb;//HIV status, hh HIV status, destined art, below before art start,self-hiv status, destined to learn
  // unsigned int state;// Susc = 0, Latent = 1,E(fast) = 2,Disease = 3,Treated = 4,X (dead) = 5
  int ARTrefusenik;		/* are they going to say no to TnT */

  int hhid;//household id, -1 if nowt
  unsigned int acat, scat;//age category, sexual activity cat
  double age, tau, taumax, taumax0, taustart, sigma, Rxtime, 
    tsi, tbout, cd4, cd42, cd4st, tshh;//hiv, tbinfection, tbactive (tboutcome). timesince HH
  unsigned int dTBsc, dTBd, dTBdet, dTBact;//destined for self-cure, death, detection, activation
  //age and hiv time-since-infection and hiv-lifetime
  unsigned int noevs;// number of events in system
  vector< double > indhaz;// hazards of events for this fella
  int is1549;//are they?
  int isS;//MTB uninfected
  int isL, isrL;//MTB latent
  int isPI;//previously infected?
  int isD, isrD;//TB disease (resistant)
  int isSmr;//smr +
  int isT;//on TB Rx
  int isX;//dead
  int isART;//on ART
  int isIPT;//on IPT
  int prevTB;//previous TB?
  int artdefr;			/* ART defaulter */
  int phhi, hhi, hhart, hhipt, tecf, iptd;//has received hh intervention in the past; on ART from HH; had IPT from HH; touched by real ECF
  double aluck;//their luck with activation times(runif), their fast luck for activation
  int fast;// are they fast from last infection? could be changed to be personal char'c
  double  tonipt;			/* time of last hazard adjustment, ipt-time */
  double LE;				/* lifeexpectancy */
  /* HIVMAC */
  double artstrate;
  int GL, CD4cat;		/* -1;1,2,3,4 */


  //functions
  void update_age( const gsl_rng *r, parameters& parz );// age ages
  void update_haz( double* foiarg, parameters& parz, double year );//update hazards
  // int event ( unsigned int etype );//carries out event --NB -- don't think actually exists or used...
  double HIVprob( parameters& parz );//relative chances of being infected
  double fertility( parameters& parz );//relative fertility
  double probsmrpos( parameters& parz );//probability smr+ at that age/hiv etc
  double probtbdetect( parameters& parz, double t );//detection probability (inc HIV, correlations)
  double detectime( const gsl_rng *r, parameters& parz, double t );//TB detection time
  void setTBoutcome( const gsl_rng *r, parameters& parz, double t );//TB sc/death time
  //  double Rinv( parameters& parz );//inverse cumulative for activation
  //  double HazReset2( parameters& parz, double H0 );//needs other things
  double a( parameters& parz );			  /* use in hazreset */
  double b( parameters& parz );			  /* use in hazreset */
  double HIVrelinfness( parameters& parz );	  /* inf'ness */

  //constructor initializer
 person(unsigned int g = 0,  double a = 23, unsigned int n = 2):
  tbrate(0.0),gender(g), hiv(0), Lhiv(0), hhhiv(0),  dart(0), bbs(0), shiv(0), dshiv(0), hhtb(0), ARTrefusenik(0), 
    acat(0), scat(0), age(a), tau(0.0), taumax(99), taumax0(99), taustart(-1.0), sigma(0.0), Rxtime(0.0), 
    tsi(0.0), tbout(99), cd4(1000), tshh(-1.0), dTBact(0),
    noevs(n), indhaz(noevs,0),
    is1549(0), isS(1), isL(0),isrL(0), isPI(0),
    isD(0),isrD(0), isSmr(0), isT(0), isX(0), isART(0), isIPT(0), prevTB(0), artdefr(0),
    phhi(0), hhi(-1), hhart(0), hhipt(0), tecf(0), iptd(0), aluck(-1),fast(0), tonipt(-1), artstrate(0), GL(-1), CD4cat(-1) {}

private:
  //nowt
};// end person class
/////////////////////////////////////



//DECLARATION HOUSEHOLD CLASS
class household{

public:
  //declarations
  list< unsigned int > members;//who's there
  unsigned int nmen, nwomen;//gender numbers
  double hfoi, rhfoi;
  unsigned int hiv, tb;//hiv status of hh, undxd tb
  //functions
  int add( unsigned int who, unsigned int sex );//someone to the household
  int remove( unsigned int who, unsigned int sex );//take someone out.

  //constructor
 household():
  members(), nmen(0), nwomen(0) {}

private:
  //nowt

};


#endif
