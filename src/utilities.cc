//#######################
#include "utilities.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//JUST FOR TESTING COMPILATION
void donowt(void){
  int i = 0;
  i++;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FOR DETERMINING THE TOTAL NUMBER OF EVENTS
unsigned int howmany_events(const gsl_rng *r, vector< double >& haz, double dt){
  double ht = 0;
  for(unsigned int i = 0; i < haz.size(); ++i){
    ht += haz[i];
  }
  double mu = dt*ht;
  return gsl_ran_poisson(r, mu);
}

// polymorphic version
unsigned int howmany_events(const gsl_rng *r, double haz, double dt){
  double mu = dt*haz;
  return gsl_ran_poisson(r, mu);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MULTINOMIAL GENERATOR FOR VECTORS
void gsl_ran_multinomial_vec (const gsl_rng * r,  unsigned int K,
                      unsigned int N,  vector< double >& p, vector< unsigned int >& n){

  if( !((K==p.size())&&(K==n.size()) ) ){cout<<"gsl_ran_mult... wanring!: sizes="<<K<<","<<p.size()<<","<<n.size()<<endl; }

  unsigned int k;
  double norm = 0.0;
  double sum_p = 0.0;

  unsigned int sum_n = 0;
  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (k = 0; k < K; k++){
    norm += p.at(k);
  }

  for (k = 0; k < K; k++){
    if (p[k] > 0.0){
      n.at(k) = gsl_ran_binomial (r, p.at(k) / (norm - sum_p), N - sum_n);
    } else {
      n.at(k) = 0;
    }

    sum_p += p.at(k);
    sum_n += n.at(k);
  }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//SPLINE-FITTING
//solves tridiagonal systems with Thomas algorithm
void tridiagsolver( vector< double >& a, vector< double >& b, vector< double >& c, vector< double >& x, vector< double >& d ){
  //see http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  /* Modify the coefficients. */
  unsigned int n = d.size();
  if (!(b.size() == x.size() && x.size() == d.size() && a.size() == c.size() && b.size() == c.size() )  ){
    cout << "tri: wrong input sizes!" << endl;
    return;
  }
  c.at(0) /= b.at(0);	/* Division by zero risk. */
  d.at(0) /= b.at(0);	/* Division by zero would imply a singular matrix. */
  for (unsigned int i = 1; i < n; i++){
    double id = 1 / (b.at(i) - c.at(i-1) * a.at(i));  /* Division by zero risk. */
    c.at(i) *= id;	                         /* Last value calculated is redundant. */
    d.at(i) = (d.at(i) - d.at(i-1) * a.at(i)) * id;
  }
  /* Now back substitute. */
  x.at(n-1) = d.at(n-1);
  for (unsigned int i = n - 2; i > 0; i--){
    x.at(i) = d.at(i) - c.at(i) * x.at(i+1);
  }
  x.at(0) = d.at(0) - c.at(0) * x.at(1);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void splinetable( vector< double >& X, vector< double >& Y, vector< double >& Z ){
  //doing calculation for natural splines
  //system passed to tri-solver is 1 less in size than no. pts
  //changes  Z to C from which the spline coeffs are computed
  if ( !(Y.size() == X.size() && X.size() == (Z.size()+1) ) ){
    cout << "splt:size mismatch in input!" << endl;
    return;
  } else if ( !Y.size()>3 ){ cout << "splt: not enough data for interpolation!" << endl; return;}

  vector< double > h( (Y.size()-1), 0 ), a( (Y.size()-1), 0 ), b( (Y.size()-1), 0 ), c( (Y.size()-1), 0 ), z( (Y.size()-1), 0 );
  for ( unsigned int ii = 0; ii < h.size(); ++ii){
    h.at(ii) = X.at(ii+1) - X.at(ii);
  }

  a.at(0) = 0;
  b.at(0) = 1;
  c.at(0) = 0;
  for ( unsigned int ii = 1; ii < (h.size()-1); ++ii){
    a.at(ii) = h.at(ii-1);
    b.at(ii) = 2 * (h.at(ii-1) + h.at(ii));
    c.at(ii) = h.at(ii);
  }
  a.back() = 0;
  b.back() = 1;
  c.back() = h.back();

  for ( unsigned int ii = 1; ii < (z.size()-1); ++ii){
    z.at(ii) = (Y.at(ii+1) - Y.at(ii)) / h.at(ii) - (Y.at(ii) - Y.at(ii-1)) / h.at(ii-1);
  }

  tridiagsolver( a, b, c, Z, z);//solve the system, z is soln

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double splinefun(double x, vector< double >& X, vector< double >& Y, vector< double >& Z ){
  //determine interval
  //Z 1 shorter than X,Y
  if ( !( x <= X.back() && x >= X.front()) ){
    cout << "splf: value not in range!" << endl;
    return 0;
  } else{
    unsigned int ii = spline_bisect_find( x, X );//which spline
    double A = Y.at(ii);
    double h = ( X.at(ii+1) - X.at(ii) );
    double B(0), C(0), D(0);
    C = Z.at(ii) ;
    if( ii+1 >= Z.size() ){	// safety for end ...
      D = ( 0 - Z.at(ii) ) / (3 * h);
      B = ( Y.at(ii+1) - Y.at(ii) ) / h - h * ( 2 * Z.at(ii) + 0 ) / 3;    
    } else {
      D = ( Z.at(ii+1) - Z.at(ii) ) / (3 * h);    
      B = ( Y.at(ii+1) - Y.at(ii) ) / h - h * ( 2 * Z.at(ii) + Z.at(ii+1) ) / 3;    
    }
    double dx = x - X.at(ii);
    return ( A + B * dx + C * pow( dx, 2 ) + D * pow( dx, 3 ) ); 
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FINDING A SINGLE EVENT (bisection utilities)
//general recursion, indices refer to bins, V[i] is value at top of bin
// V must be an increasing sequence
int bisect( double x, vector< double >& V, int lo, int hi){
  if (  x > V.back() ){//NB no check below
    cerr << "bisect range error!" <<  endl;
    cerr << "x,Vf,Vb: "<< x <<","<< V.front()<<","<<V.back()<<endl;
    exit(-1);//bail
    return -1;
  } else {
    //cout << "lo=" << lo << ",hi="<<hi << endl;
    if(lo == hi) {return lo;}
    int mid = (lo + hi) / 2;
    if ( x <= V.at(mid) ){
      return bisect( x, V, lo, mid );
    } else {
      return bisect( x, V, mid+1, hi);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//wrapper for bisection in event-context
// V must be an increasing sequence
unsigned int bisect_find( double x, vector< double >& V){
  int ans = bisect( x, V, 0, (V.size()-1) );
  if( ans == -1){ 
    return 0;
    cerr << "x,Vf,Vb: "<< x <<","<< V.front()<<","<<V.back()<<endl;
    exit(-1);//bail
  } else { return (unsigned int) ans;}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//wrapper for bisection in spline context
unsigned int spline_bisect_find( double x, vector< double >& V){
  if ( x < V.front() ){
    cerr << "below bisection range!" << endl;
    return 0;
  } else {
    int ans =  ( bisect( x, V, 1, (V.size()-1) ) - 1 );
    if ( ans < 0){ return 0;}
    return (unsigned int) ans;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FINDING MULTIPLE EVENTS
//version using a set to test whether done or not
//working with this setted version
void find_events(const gsl_rng *r, list< unsigned int >& ES, unsigned int noeli, unsigned int noev, vector< double >& W ){
  double u;
  double totw;
  unsigned int k;
  set< unsigned int > done;
  vector< double > newW( noeli, 0);//slow?
  newW[0] = W[0];
  for( k = 1; k < noeli; ++k ){newW[k] = newW[k-1] + W[k];}// a shorter thing to search through, which is cumulative
  totw = newW.back();//the total weight
  //DEBUG
  if( W[0] < 0 || totw < 0 ){
    cerr << "in find_events"<<endl;
    cerr << "noeli,noev= "<<noeli<<","<< noev <<endl;
    cerr << "Wb,Wf= "<< W.back()<<","<<W.front()<<endl;
    cerr << "newW[0],totw= "<< newW[0]<<","<<totw<<endl;
    cerr << "k,W.size()= "<< k<<","<< W.size()<<endl;
    exit(-1);
  }
  //cout << totw <<"<"<<endl;
  if ( noev > noeli) {noev = noeli;}//no double hits implies this
  unsigned int ui = 0;
  int t(0), a(0);
  //parallelize??
  while( ui < noev ){
    ++t;
    u = gsl_ran_flat( r, 0, totw );
    k = bisect_find( u, newW );
    //u = gsl_rng_uniform(r);//runif
    //k = bisect_find( totw*u, newW );//find single event by iterated bisection - scaled to avoid normalisation
    if ( done.count( k ) == 0 ){//not previously done -- ?
      done.insert( k );//add  done list (possibly costly and nb duplication with event list)
      ES.push_back( k );//place in event list
      ++ui;//advance event count
      ++a;//a/t is accpetance rate
    }
    if ( t > 50 * (int) noeli){//arbitrary time-out
	cerr << "bailing from event finder (timeout)!" << endl;
	cerr << "done.size()="<<done.size()<<",ES.size()="<<ES.size()<<",ui="<<ui<<endl;
	cerr << "k=" <<k<< ",W_k=" <<W[k] << ",|d|=" <<done.count(k)<< endl;
	cerr << "totw="<< totw <<endl;
	cerr << "noeli="<<noeli<<",noevs="<<noev<< endl;
	//cerr<<"ES=";for(list< unsigned int >::iterator ei = ES.begin(); ei != ES.end(); ++ei){cerr<<*ei<<",";}cerr<<endl;
	exit (-1);//break
      }
  }//end while
  //cout << "done while! acceptance="<< 100*a/t << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FINDING A SINGLE EVENT (WEIGHTED)
//wrapper for above: NB polymorphic
//not nec efficient...but convenient
int find_one_event(const gsl_rng *r, vector< double >& W ){
  list< unsigned int > ES;
  ES.clear();
  find_events(r, ES, W.size(), 1, W );
  int answer = (int) ES.back();
  return answer;
}

int find_one_event(const gsl_rng *r, vector< int >& Wi ){
  vector< double > W ( Wi.begin(), Wi.end() );//casting
  list< unsigned int > ES;
  ES.clear();
  find_events(r, ES, W.size(), 1, W );
  int answer = (int) ES.back();
  return answer;
}

int find_one_event(const gsl_rng *r, vector< unsigned int >& Wu ){
  vector< double > W ( Wu.begin(), Wu.end() );//casting
  list< unsigned int > ES;
  ES.clear();
  find_events(r, ES, W.size(), 1, W );
  int answer = (int) ES.back();
  return answer;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FOR RANDOM PERMUTATIONS
void shuffle( const gsl_rng *r, vector< int >& w ){//from Knuth vol.2 pg 145
  unsigned int t = w.size();
  unsigned int j = t;
  unsigned int k;
  double u;
  int temp;
  while ( j > 1){
    u = gsl_rng_uniform(r);//runif
    k = auxFloor( j*u ) + 1;
    //swap w_k and w_j
    temp = w.at( j - 1 );//w_j
    w.at( j - 1 ) = w.at( k - 1 );//w_j=w_k
    w.at( k - 1 ) = temp;//w_k=old w_j 
    --j;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MOVING AVERAGES (overloaded)
void MAsmooth( vector< double >& V, vector< double >& W, int L ){//doesn't overite
  //L->2L+1 window
  W.assign( V.size(), 0 );
  for( int u = 0; u < (int)V.size(); ++u ){
     W.at(u) = V.at(u);
     int count = 1;
     for ( int i = 1; i <= L; ++i ){
       if( (u - i) >= 0 ){
	W.at(u) += V.at( u - i );
	++count;
      }
       if( (u + i) < (int)V.size() ){
	W.at(u) += V.at( u + i );
	++count;
      }
    }
    W.at(u) /= count;
  }
}

void MAsmooth( vector< double >& V,  int L ){//overwrites
  //L->2L+1 window
  vector< double > W( V.size(), 0.0);
  for( int u = 0; u < (int)V.size(); ++u ){
     W.at(u) = V.at(u);
     int count = 1;
     for ( int i = 1; i <= L; ++i ){
       if( (u - i) >= 0 ){
	W.at(u) += V.at( u - i );
	++count;
      }
       if( (u + i) < (int)V.size() ){
	W.at(u) += V.at( u + i );
	++count;
      }
    }
    W.at(u) /= count;
  }
  V = W;//overwrite V
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// overloaded version to deal with things where should be zero but-for
// another quantity
// NB should be applied before any averaging of second quantities!
int collapse_result_averager( vector< vector< double > >& ibmres, vector< vector< double > >& ibmres2 ){
  unsigned int m = ibmres[0].size(); // how many times
  unsigned int n = ibmres.size();    // how many runs
  if( ibmres2.size()!=n || ibmres2[0].size()!=m ){cout <<"**size mismatch in results averager**...bailing!"<<endl; exit(-1);}
  unsigned int j(0), i(0);
  double av(0);
  int denom(0);
  for( i = 0; i < m; ++i ){	// time loop
    av = 0; denom = 0;
    for( j = 0; j < n; ++j ){	// run loop
      if( ibmres2[j].at(i) >= 0 && ibmres[j].at(i) >= 0  ){ 
	// -1s for rare events
	av += ibmres[j].at(i); ++denom;
      }
    }
    if( denom > 0 ){
      ibmres[0].at(i) = av / denom;		// save back average to first 
    } else {ibmres[0].at(i) = -1;}
  }
  return (int)n;
}


int collapse_result_averager( vector< vector< double > >& ibmres ){
  int n = collapse_result_averager( ibmres, ibmres );
  return (int)n;
}

/////////////////////////////////////////////////////END OF UTILITY FUNS/////////////////////////
