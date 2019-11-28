//#######################
#include "parameters.h"
#include "results.h"
#include "population.h"

ofstream lgf;


int main(int argc, char **argv){
 
  //TAKES A DIRECTORY PREFIX 
  // N.B. jobid now used to set rng seed
  cout << argc << endl;
  string wdir;			// the working directory
  int jobid;		        // integers to uniquify output
  switch( argc ){
  case 2:
    wdir.assign( argv[1] );	// assign working directory
    jobid = 1;	// job id
    break;
  case 3:
    wdir.assign( argv[1] );	// assign working directory
    jobid = atoi( argv[2] );	// job id
    break;
  default:
    cout << "**error**! need 1 (working full-path-directory/) argument!" << endl;
    cout << "...and optionally one job id ...bailing!" << endl;
    exit(-1);
  } // end switch for aruments
  cout <<"..."<<jobid<< "...working in: " << wdir << endl;

  //OPENING LOGFILE
  ostringstream oss;
  time_t starttime;
  time( &starttime );
  oss << wdir << jobid << "logfile.txt";
  lgf.open( oss.str().c_str() ); // unique output 
//   lgf.open( (wdir + "logfile.txt").c_str() );
  if (!lgf){cout << "failed to open logfile!" << endl; exit (-1);}
  lgf << "logfile opened: timestamp = " << starttime << endl;
  lgf << "working directory = " << wdir << endl;

  //PARALLEL FLAGS
#ifdef PARA
  lgf << "PARA flag turned on" << endl;
  cout << "PARA flag turned on" << endl;
#endif 

 //GSL RNG INITIALIZATION
  lgf << "initializing GSL rngs...";
  const gsl_rng_type *T;
  gsl_rng *r;
  T = gsl_rng_taus;//tausworthy //gsl_rng_default;//type of rng
  //T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);//
//   //time-seeded
//   time_t seconds;
//   seconds = time(NULL);
//   gsl_rng_set(r,seconds);//the seed
  lgf << "done!" << endl;

  //READ PARAMETER FILES and DECLARATIONS
  lgf << "reading parms..." << endl;
  parameters parz;// the parameters - parameter object
  parz.read( lgf, wdir );// read in parameters
  parz.jobid = jobid;  // apply cmdline args
  parz.computerest();//compute other dependent parz
  parz.readpopdata( r, lgf );//read in population data and initialize

  population pop( 0, 1, 4, 4, parz.nbins);//declare the population: size, hhs, noevs, nostates, abins

  bool light; int nr;
  if( parz.lightFLAG ){
    light = true; nr = 1;	// initialize light
  } else {
    light = false; nr = parz.noruns; // initialize heavy
  }
  results rez( nr, parz.noruns, parz.notime, light );//initialize

  // LOOP OVER RUNS
  //gets wiped at the start of each loop
  //record is ok as it uses run no. 
  int runno = 0;
  for ( runno = 0; runno < parz.noruns; ++runno ){//run loop
    cout << parz.wkdir <<" runno= "<< runno << endl;
    unsigned long int seed = (unsigned long int)(10000*parz.jobid + runno);
    gsl_rng_set(r,seed);//the seed
    lgf << "starting with seed="<< seed << "..."<<endl;

    //setup
    pop.reinitialize( r, lgf, parz );// initialize from read in data
    pop.initialinfect( r, lgf, parz );//setting up the initial infections

    //snapshot
    bool snapshot = parz.snaps;	     // controls both
    if ( parz.hhFLAG && runno == 0 && snapshot ){//for first run
      pop.demosnapshot( lgf, rez, 0, 1, parz);
    }

    //runs
    pop.ibmrun( r, lgf, parz, rez, runno );//do model run
    
    //snapshot
    if ( parz.hhFLAG && runno == 0 && snapshot ){//for first run
      pop.demosnapshot( lgf, rez, 0, 2, parz);
    }
  }//end run loop
  
  //PROCESS AND RECORD RESULTS
  rez.analyze( lgf, parz );
  time_t writetime;
  time( &writetime );//time for unique output mode
  rez.write( lgf, parz, writetime );//writing simple results

  //TIDYING
  gsl_rng_free(r);
  lgf << "freeing gsl rng..." << endl;
  lgf << "exiting..." << endl;
  lgf.close();
  return(0);
}
