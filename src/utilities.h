//#pragma once			/* guard */
#ifndef __UTILITIES_INCLUDED__
#define __UTILITIES_INCLUDED__
#define auxFloor(x) ((float)(long)(x))
#define auxCeil(x) ((float)(long)((x)+1))
/* fwd declarations */
/* includes */
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <list>
//#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// declarations std lib
using namespace std;

//FUNCTION DECLARATIONS
void donowt(void);
unsigned int howmany_events(const gsl_rng *r, vector< double >& haz, double dt);
unsigned int howmany_events(const gsl_rng *r, double haz, double dt);
void gsl_ran_multinomial_vec (const gsl_rng * r,  unsigned int K,  unsigned int N, vector< double >& p, vector< unsigned int >& n);
void tridiagsolver( vector< double >& a, vector< double >& b, vector< double >& c, vector< double >& x, vector< double >& d );
void splinetable( vector< double >& X, vector< double >& Y, vector< double >& Z );
double splinefun(double x, vector< double >& X, vector< double >& Y, vector< double >& Z );
int bisect( double x, vector< double >& V, int lo, int hi);
unsigned int bisect_find (double x, vector< double >& V);
unsigned int spline_bisect_find( double x, vector< double >& V);
void find_events(const gsl_rng *r, list< unsigned int >& ES, unsigned int noeli, unsigned int noev, vector< double >& W );
int find_one_event(const gsl_rng *r, vector< double >& W );
int find_one_event(const gsl_rng *r, vector< int >& Wi );
int find_one_event(const gsl_rng *r, vector< unsigned int >& Wu );
void shuffle( const gsl_rng *r, vector< int >& w);
void MAsmooth( vector< double >& V, vector< double >& W, int L );
void MAsmooth( vector< double >& V,  int L );
int collapse_result_averager( vector< vector< double > >& ibmres, vector< vector< double > >& ibmres2 );
int collapse_result_averager( vector< vector< double > >& ibmres );


#endif
