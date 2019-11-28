//#pragma once			/* guard */
#ifndef __PARAMETERS_INCLUDED__
#define __PARAMETERS_INCLUDED__

#include "utilities.h"

//DECLARATION PARAMETERS CLASS
class parameters{

 public:
  int noruns, jobid;
  double starttime, stoptime, dt, time;		     /* copy of the time */
  double fertage_alpha, fertage_beta, mf_deathratio;//fertility by age, and death
  double beta;//infection strength
  double lamn, nun;//selfcurerate, tbdeathrate -- not sure used any more...
  double r0, r1, deltan, deltap;// reactivation,  actual rate for CD
  double CDRn, RxL, Rxp, Rxpd, Iprotn, IprotnS, Pprotn, PprotnS, deltanSmr;//prob E not L, CDR, duration Rx, prob Rx success, protection if L, protection after IPT, detection smr+
  double  fSmr,  tbsurvn_k,tbsurvn_L,tbmortn, hivdurf;//, tbsurvp_k, tbsurvp_L;//TB self-cure probs and survival timescales
  //tbmortn is prob tb death norx, tbsurvx are timescales //probSCp,
  double smrfac;	       /* reduction in smr+ty hiv- */
  double tbd_betan, tbd_betap;//detection timings - fraction speedups
  double H2Lhr;		      /* hr for HIV 2 LTBI - more clustering */
  double tbd_betan2, tbd_betap2;//detection timings w/interventions
  int S0, L0, D0, T0;//initial states 
  int notime, popsize0;//to be computed
  double thiv_alpha, thiv_beta;// HIV survival parms
  double mhiv_alpha, mhiv_beta, fhiv_alpha, fhiv_beta, fm_hivratio, cd40;//M/F HIV and by age parms
  /* cd4 at start time */
  double h_t0, h_peak, h_alpha, h_beta, h_theta, h_peaktime, h_peakiness, decline2010;
  //HIV incidence SM model
  double  f, g, rho, CDRp, Rxpp, Rxppd, deltapSmr, shivt, shivp;
  int massart;			/* TnT round */
  double artstart, artT, artP;	/* start time, period and prob refusenik */
  double cd4bl;			/* baseline cd4 */
  double IPTdurnN, IPTdurnP, IPThrN, IPThrP, IPThrA, IPTcprobN, IPTcprobP; 
  /* IPT characteristics */
  int qlog, lightFLAG, snaps, bigSS, Rxflag, Rflag;	      /* switches */
  int DHflag;
  double sbFracH, sbFracM, sbFrac1, sbFrac2, sbFrac3, sbc1, sbc2, sbc3, sba1, sba2, sba3, sbA, lam, riA, riE, durA, durE, ARTeps, hiv0; /* dynamic HIV */
						      /* switches */
  double ARTf, cd42a, cd42b, cd42c, expET, artTk, artTl, dCD4et, dCD4ep, artDR, artCU;
  /* ART modelling parms */
  double dCD4st, cd41;		/* not read in - set in computerest */
  double Rxtime, Rxp2, Rxpp2;	/* changes in Rx */
  double Rfrac0;
  /* IPT ramp-ups */
  double IPTcovP2, IPTcovPt1, IPTcovPt2, IPTcovP1;
  /* HIVMAC stuff */
  double artsig, pscale, pn, pmaxval, kq, gg, gi, glrst, nglstr, hseed,  SQGLprop;
  int hivmac, SQ;
  int iGLRc;			/* oddly, used as counter, HIVMAC */
  int artp15x, artn15x, hivn15x, u15x, ihiv15, ihiv1549;
  int hivpgp15x, hivpgn15x, iARTgp, iARTgn, iHIVgp, iRx;
  double GLfracn,GLfracpn,GLfract, GLfracnoart, GLfracnoteli;
  double prevART, hivu200, X15tbn, X15tban, X15tbap;
  double rxt0, rxt1, rxstart, rxrate, mdrprop0, mdrprop1, mdrt0, mdrt1, mdrretxOR;    

  vector< unsigned int> intsex, inthh;//initial population pars
  vector< double > intage;//initial population pars
  vector< double > mr_dat, mr_yrz, br_dat, br_yrz, dr_dat, dr_yrz, h1549_yrz, h1549_dat, cdr_yrz, cdr_dat, rpop_yrz, rpop_dat;//data on birth/death/migration/hiv
  vector< double > mr_spl, br_spl, dr_spl, h1549_spl, cdr_spl, rpop_spl;//data for spline interpolations
  //vector< vector< double > > hhsizes;//matrix of hh sizes
  //vector< double > hhtotsizes;//vector of same but only total hh sizes
  //vector< double > f_weight, a_weight;//vectors for detailed balance hh method
  vector< double > agebins;//vector of age bin tops (last one irrelevant)
  vector< vector< double > > agedata;//vector of age-mixing data
  vector< vector< vector< double > > > rxdata0, rxdata1;//treatment outcome data - first is files
  vector< vector< double > > UNLTb;//UN life-table
  vector< vector< int > > popdata;//initial population data
  unsigned int nbins;//number of age categories
  int nohh0;//initial no hh -- won't be read but calculated
  double nophh0, hhmoverate;//initial no/hh, household moverate
  int hhFLAG, AGEflag;//hhflag, flag for agestructured infection
  double betaH, hhhivRR;//household infection parm, hh hiv risk ratio
  int tECFflag, ECFflag, IPTflag, ARTflag, HHIflag, HIPflag, HARflag;//intervention flags
  double tECFst, tECFet,ECFst, ECFet, IPTst, IPTet, ARTst, ARTrt, tbart, HHIst, HHIet;//start/end-times and tbart is spectrum of disease
  double hhdOR, hhdF, deltaSmrI;//past HH detection OR, factor reduction in time-t-detect, differential intervention effect by smr status
  double CDRn2, CDRp2, IPTcov, IPTT, IPTart, IPTshiv, IPTcovP, IPTpprot, ARTcovmax, HHIcov, TSTu16sens, cd4h, cd4g, cd4g2;//other int parz
  int IPTmultiple, art750;	   /* multple courses allowed?, art cd42=750 override */
  /* cd4s typical, guideline and household */
  double deltan2, deltap2, tecfdOR, tecfdF;//detection rates for ECF; true ECF OR and factor
  double maxdOR, maxdF, tECFhaz; //the maximum of effects when both apply, tECF application rate
  double EVA,EVB,EVC,EVD, fscale, smrnmort;//Emilia's parz
  string wkdir, birthFN, deathFN, migrateFN, amFN, synpop0, LTbFN, abFN, cdrFN, rpopFN;//filenames for external data
  double ART2TB, ART2TBst;
  int ART2HHhx;		   /* extra interventions */
  double ARThxT;
  double PACFlag, PACFst, PACFet, PACFT, PACFeffsn, PACFeffsp, PACFeffart, PACFhivOR; /* ACF parms */
  double HIVflag, HIVst, HIVet, HIVT, HIVeff; /* HIV screening parms */
  int hivincF;				      /* using HIV inc file? */
  string hivincFN, rxFN0, rxFN1;			      /* location if yes */
  double partrefusenik;/*probability of become a permanent refusenik*/
  double TSTspec, TSTsens;


  //functions
  int computerest( void );
  int read( ofstream& lgf, string& wdir );//reading parameter files
  int readpopdata( const gsl_rng *r, ofstream& lgf );//read the data only
  double br ( double yr );//birth rate
  double dr ( double yr );//death rate
  double rpop( void );		/* rel population */
  double mr ( double yr );//migration  rate
  double h1549 ( double yr );//hiv incidence 15-49
  double fastrisk( double x );//age-dept fast risk
  double A( double agee, double agor );//age-mixing function
  double ARTcov( double t );//ART coverage at a given point in time
  double extratime( double lifeleft, double cd4 ); /* life expectancy starting art */
  double Hexp( double a, double b, double x );
  double invHexp( double a, double b, double x );
  double lifeexpectancy( const gsl_rng *r, double E0, double age );
  double lifeexpectancy( const gsl_rng *r, double age ); /* wrapper */
  double ARTcd4asymptote( double cd4A );				 /* cd42 */
  double Rxpfn(); double Rxppfn();	
  /* HIVMAC */
  double GL_g( void );double GL_p( void );double GL_q( void );
  double GL_newp( double frac );  double GL_newq( double frac );
  /* new changes for interventions */
  double IPTcovPf( void );
  double rxINTcov( void );
  double rxMDRprop( int retx );
  double rxprob( int reg, int rxst, int n );
 private:
  //nowt
 
};//end parameters class
/////////////////////////////////////



#endif
