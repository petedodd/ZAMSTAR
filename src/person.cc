//#######################
#include "person.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PERSON MEMBER FUNS

void person::update_age( const gsl_rng *r,  parameters& parz ){
  double dt = parz.dt;
  if( !isX ){//alive! 
    age += dt;
    if( parz.AGEflag ){
      if( (age > parz.agebins.at(acat) ) && ( acat < parz.nbins - 1 ) ){
	++acat;
      }

    }
    if ( hiv ){//HIV+
      tau += dt;
      if( isART ){		// on ART -taumax0 and taustart recorded
	cd4 = cd42 * ( 1 - parz.ARTf + parz.ARTf * (taumax - tau) / (taumax - taustart));
	if( cd4 < 0 ){cd4 = 0;}			     // safety
      } else {	
	if( artdefr == 0 ){		// NOT ART defaulter
	  // including ART! not on ART
	  double frac = ( taumax - tau ) / taumax;
	  //cd4 = parz.cd40 * 0.75 * frac;//linear
	  cd4 = parz.cd40 * 0.75 * pow( frac, 2 ); // quadratic
	  if( cd4 < 0 ){cd4 = 0;}			     // safety
	} else {	  //ART defaulter
	  double frac = ( taumax - tau ) / (taumax0 - taustart);
	  //cd4 = cd4st * frac ; // linear
	  cd4 = cd4st * pow( frac, 2 ) ; // quadratic
	  if( cd4 < 0 ){cd4 = 0;}			     // safety
	}
      }

      // CD4 categories
      if( hiv && !isART ){
	if( cd4 < 200 ){
	  CD4cat = 4;
	} else{
	  if( cd4 < 350 ){
	    CD4cat = 3;
	  } else{
	    if(cd4 < 500){
	      CD4cat = 2;
	    } else{ CD4cat = 1;} 
	  }
	}
      }	// CD4 cats

    }
    if( isL ){//MTB infected
      sigma += dt;
    }
    if( isT ){//on Rx
      Rxtime += dt;
    }
    if( isD ){//has active undetected TB
      tsi += dt;
    }
    if( isIPT ){//on IPT
      tonipt += dt;
    }

    // household intervention
    if( parz.HHIflag ){
      if( hhi ){// currently household
        if (tshh <=0.5){   
          tshh += dt;
        } else {
          hhi = 0;              // no longer currently on
        }
      }
    }

    // HIVMAC stuff here
    if( parz.hivmac ){

      // put extra HIVMAC stuff here
      if( age > 15 && GL == -1 ){
	// this is aging in under the assumption things happen on becoming 15
	// only happens once
	if( gsl_rng_uniform(r) < parz.GL_g() ){
	  GL = 1;
	  if( hiv ){++parz.iHIVgp;}
	  ++parz.iGLRc;
	} else { GL = 0; }
      }

      bool logistic = true;
      if(!logistic){ 		// the original parametric version
	// these bits represent the hazards at each point in time
	if( !hiv && GL == 0 ){
	  if( gsl_rng_uniform(r) < 1-exp(-parz.GL_p() * parz.dt) ){
	    GL = 1;	// conversion while HIV-
	    ++parz.iGLRc;
	  }
	}

	if( hiv && GL == 0 ){
	  if( gsl_rng_uniform(r) < 1-exp(-parz.GL_q() * parz.dt) ){
	    GL = 1;	// conversion while HIV+
	    ++parz.iHIVgp;
	    ++parz.iGLRc;
	  }
	}
      } else {			// new logistic version

	// new logistic piece
	if( GL == 1 ){		// down-converstion for logistic model
	  if( hiv ){
	    if( gsl_rng_uniform(r) < 1-exp(+parz.GL_q() * parz.dt) ){
	      GL = 0;	// conversion while HIV+
	    }
	  } else {
	    if( gsl_rng_uniform(r) < 1-exp(+parz.GL_p() * parz.dt) ){
	      GL = 0;	// conversion while HIV-
	    }
	  }
	}

	if(GL == 0){		// conversion
	  if(hiv){
	    if( gsl_rng_uniform(r) < 1-exp(-parz.GL_q() * parz.dt) ){
	      GL = 1;	// conversion while HIV+
	      ++parz.iHIVgp;
	      ++parz.iGLRc;
	    }
	  } else {
	    if( gsl_rng_uniform(r) < 1-exp(-parz.GL_p() * parz.dt) ){
	      GL = 1;	// conversion while HIV-
	      ++parz.iGLRc;
	    }
	  }
	}
      } // end else


      // re-evaluate the art start rate
      if( hiv && !isART && GL >=0 ){

	double guide = (parz.time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
	if( GL == 1 ){		// guideline ready
	  // hopefully the above will get optimized away...
	  if( cd4 < guide ){
	    artstrate = 2.0;	// this is the set rate to start
	  } else {
	    artstrate = 0.0;	// other rates for starting
	  }
	} else {		// not guideline ready
	  if( cd4 < guide ){
	    artstrate = parz.nglstr;   // this is the set rate to start
	  } else {
	    artstrate = 0.0;	// other rates for starting
	  }
	} // end else
      }	  // end HIV+ ART- GLe

      // end of new HIVMAC stuff....


    } // end HIVMAC


  } // end not dead
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::HIVprob( parameters& parz ){//relative chances of being infected
  double ans = 0;
  if ( ( hiv == 0 ) && ( !isX ) ){//not infected or dead
    ans = gsl_ran_weibull_pdf( age, parz.mhiv_beta, parz.mhiv_alpha );
    // switch( gender ){//insert age/sex biases here...todo
    // case 0://women
    //   if ( age >= 15 ){
    // 	ans = parz.fm_hivratio * gsl_ran_weibull_pdf( (age-15), parz.fhiv_beta, parz.fhiv_alpha );
    //   }
    // case 1://men
    //   if ( age >= 15 ){ 
    // 	ans = gsl_ran_weibull_pdf( (age-15), parz.mhiv_beta, parz.mhiv_alpha );
    //   }
    // }
    // here for H2Lhr
    ans *= parz.H2Lhr;
  }
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::fertility( parameters& parz ){//relative chances of giving birth
  double ans = 0;
  if ( ( gender == 0 ) && ( !isX ) ){//not male or dead
    ans = gsl_ran_weibull_pdf( (age-15), parz.fertage_beta, parz.fertage_alpha );
  }
  return ans;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void person::update_haz( double* foiarg, parameters& parz, double year ){//update hazards
  indhaz.assign( indhaz.size(), 0.0 );
  double foi = foiarg[1];
  double rfoi = foiarg[3];
  foi -= rfoi;

  if ( !isX ){//not dead
    // Mtb 
    if ( hiv == 0 ){//HIV -ve
      if ( isS || isL ){//MTB infection
	indhaz.at(1) =  foi;
	if ( isL ){ 
	  indhaz.at(1) *= ( 1 - parz.Iprotn );
	} else {
	  if ( isPI ){//IPT, or successful Rx but PI
	    indhaz.at(1) *= ( 1 - parz.IprotnS );
	  }
	}
      }
    } else {//HIV +ve
      if ( isS || isL ){//MTB infection
	indhaz.at(1) =  foi / parz.g;
	if ( isL ){ 
	  indhaz.at(1) *=  ( 1 - parz.Iprotn );
	} else {
	  if ( isPI ){//IPT but PI
	    indhaz.at(1) *=  ( 1 - parz.IprotnS );
	  }
	}
      }
    }//end else HIV

    // rMtb
    if( parz.Rflag ){ 
      if ( hiv == 0 ){//HIV -ve
	if ( isS || isL ){//rMTB infection
	  indhaz.at(3) =  rfoi;
	  if ( isL ){ 
	    indhaz.at(3) *= ( 1 - parz.Iprotn );
	  } else {
	    if ( isPI ){//IPT, or successful Rx but PI
	      indhaz.at(3) *= ( 1 - parz.IprotnS );
	    }
	  }
	}
      } else {//HIV +ve
	if ( isS || isL ){//MTB infection
	  indhaz.at(3) =  rfoi / parz.g;
	  if ( isL ){ 
	    indhaz.at(3) *=  ( 1 - parz.Iprotn );
	  } else {
	    if ( isPI ){//IPT but PI
	      indhaz.at(3) *=  ( 1 - parz.IprotnS );
	    }
	  }
	}
      }//end else HIV
    } // end flag

  }//end not dead 


  // HIV infection
  double hivfoi = foiarg[0];
  if ( !hiv && !isX && is1549 ){//not infected/dead & 15-49
    indhaz.at(0) = hivfoi * HIVprob( parz );
    if ( parz.hhFLAG ){
      if ( hhhiv > 0 ) { indhaz.at(0) *= parz.hhhivRR; }// hh enhancement
    }
  }

  // tb rates in new method
  tbrate = 0.0;
  if( !isX && !isD && !isT && isL ){ // this matters!

    if( sigma < 2.0 ){		// alternative changing hazard COMMENT!
      tbrate = parz.r1 * parz.fastrisk( age );
    } else { tbrate = parz.r0;}
    if( hiv ){
      // works for both on/no-ART as long as we keep track of cd4-counts
      tbrate *= exp((parz.cd40 - cd4) * parz.rho * 1e-2);
      if( isIPT ){		// on IPT
	if( isART ){
	  tbrate *= parz.IPThrA;
	} else {
	  tbrate *= parz.IPThrP;
	}
      }				 // no IPT
    } else {			// HIV-
      if( isIPT ){ tbrate *= parz.IPThrN; }
      if( !dTBact ){tbrate = 0;} // not destined, i.e. protected
    }
  } // end not dead or S
  indhaz.at(2) = tbrate;	// unnecessary duplication
  // debug!
  // if( fabs(year - 2000) <.05 ){
  //   cout <<tbrate <<","<<isL<<","<<hiv<<","<<dTBact<<","<<isART<<","<<isIPT;
  //   cout <<","<<cd4<<","<<age<<","<<sigma<<","<<fast<<","<<tau;
  //   cout<<","<<year<<endl;
  // }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::probsmrpos( parameters& parz ){
  double A(0.1), B(0.65), C(0.87);//yvals
  double a(10), b(20), c(90);//xvals
  //  A = 0.0;			// NB -- experiment!
  double ans = A;
  if( age < b && age >= a ){
    ans = A * ( b - age ) / ( b - a ) + B * ( age - a ) / ( b - a );
  }
  if ( age >=b && age < c ){
    ans = B * ( c - age ) / ( c - b ) + C * ( age - b ) / ( c - b );
  }
  if ( age  > c ){ ans = C; }
  //reduction in probability of being smear +ve -- this could depend on tau...
  if( Lhiv ){ 
    ans *= parz.f; 
  }
  ans *= parz.smrfac;		// overall fudge factor
  return ans;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::probtbdetect( parameters& parz, double t ){
  //the probability of detection with TB | HIV, smr status etc...
  // see also tedious redestining if hit by ECF while prevalent
  double ans(0);
  if( parz.cdrFN == "none" ){				  // not reading in
    if ( parz.ECFflag && t > parz.ECFst && t < parz.ECFet){//CDR+
      if( Lhiv ){
	ans = parz.CDRp2;
      } else {
	ans = parz.CDRn2;
      }
      if( !isSmr ){ ans = (parz.deltapSmr * ans / (1 + (parz.deltapSmr-1) * ans));}//smr-ve OR
    } else {//normal CDR
      if( Lhiv ){
	ans = parz.CDRp;
      } else {
	ans = parz.CDRn;
      }
      if( !isSmr ){ ans = (parz.deltanSmr * ans / (1 + (parz.deltanSmr-1) * ans));}//smr-ve OR
    }
  } else {			// readin
    if( t <= parz.cdr_yrz.front() ){ ans = parz.cdr_dat.front() / 100;} else {
      if( t >= parz.cdr_yrz.back() ){ ans = parz.cdr_dat.back() / 100; } else{
	ans = splinefun( t, parz.cdr_yrz, parz.cdr_dat, parz.cdr_spl ) / 100;
      }
    }
    if( !isSmr ){ ans = (parz.deltanSmr * ans / (1 + (parz.deltanSmr-1) * ans));}//smr-ve OR
    // override with ECFflag!!
    if ( parz.ECFflag && t > parz.ECFst && t < parz.ECFet){//CDR+
      if( Lhiv ){
	ans = parz.CDRp2;
      } else {
	ans = parz.CDRn2;
      }
      if( !isSmr ){ ans = (parz.deltapSmr * ans / (1 + (parz.deltapSmr-1) * ans));}//smr-ve OR
    }
  } // end readin choice



  //HH/real ECF interventions
  if( phhi && !tecf  ){//has had a past hh intervention/no ecf
    ans = (parz.hhdOR * ans / (1 + (parz.hhdOR-1) * ans));//OR for detection having had HH
  }
  if( tecf && !phhi ){//real ECF only
    ans = (parz.tecfdOR * ans / (1 + (parz.tecfdOR-1) * ans));//OR for detection touched ECF
  }
  if( tecf && phhi ){//subject to both -- decide how to combine
    ans = (parz.maxdOR * ans / (1 + (parz.maxdOR-1) * ans));//OR for detection both
  }
  // differential effect of interevention by smr-status
  if( ( tecf || phhi ) && isSmr ){ // smr+ invervention
    ans = ( parz.deltaSmrI * ans / (1 + (parz.deltaSmrI-1) * ans));
  }

  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::detectime( const gsl_rng *r, parameters& parz, double t ){
  //the distribution of detection times from TB | HIV status etc...
  //nb this is conditional on detection...(of course)
  //can change to depend on smr+/-
  double ans(0);

  // above is old version; changing to multiply by betas
  ans = tbout;
  if ( parz.ECFflag && t > parz.ECFst && t < parz.ECFet){//CDR+
    if( hiv ){
      ans *=  parz.tbd_betap2;
    } else {
      ans *= parz.tbd_betan2;
    }
  } else {//normal CDR
    if( hiv ){
      ans *=  parz.tbd_betap;
    } else {
      ans *= parz.tbd_betan;
    }
  }
  // multiplies detection time
  if( phhi && !tecf ){	// has been subjected to HH intervention
    ans *= parz.hhdF / parz.deltaSmrI;//speed-up
  }
  if( !phhi && tecf ){//touched ECF only
    ans *= parz.tecfdF / parz.deltaSmrI;//speed-up
  }
  if( phhi && tecf ){//both
    ans *= parz.maxdF / parz.deltaSmrI;//speed-up
  } 

  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void person::setTBoutcome( const gsl_rng *r, parameters& parz, double t  ){
  dTBd = 0; dTBsc = 0; dTBdet = 0;//safe
  //destined for death or sc?

  tbout = gsl_ran_weibull( r, parz.tbsurvn_L, parz.tbsurvn_k );
  if( Lhiv ){//"Like HIV" ART dependent
    //tbout = gsl_ran_exponential( r, parz.tbsurvp );
    tbout *= parz.hivdurf;
    dTBd = 1;
  } else {

    //check order and include
    //if negative these are the survival parameters as case
    double test = ( isSmr? parz.tbmortn : parz.smrnmort );
    if( gsl_rng_uniform(r) <  test ){//death
      dTBd = 1;
    } else{//self-cure
      dTBsc = 1;
    }
    //then consider detection...
  } 
  //this is proportional hazards for death self-cure 

  //destined for detection?
  if( gsl_rng_uniform(r) < probtbdetect( parz, t ) ){
    dTBdet = 1;//destined
    dTBsc = 0; 
    dTBd = 0;
    tbout = detectime( r, parz, t );
    //when? - outcome time
  }
  // DEBUG
  // int out = dTBdet*1 + dTBsc*2 + dTBd*3;
  // cout << tbout << "," << hiv <<","<< out <<","<<isSmr<<","<<Lhiv<< endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::a( parameters& parz ){ // functions on way to ramp-ups
  double ans(0);
  if( hiv ){
    if( isART ){
      // ans = parz.cd40 - (1-parz.ARTf) * cd42;
      // ans += - parz.ARTf * cd42 * taumax / (taumax - taustart) ;
      // ans *= parz.rho * 1e-2;
      ans = ( parz.cd40 - cd42 ) * parz.rho * 1e-2;
      if( isIPT ){ ans += log( parz.IPThrA ); }   // scaling from IPT
    } else {
      ans = parz.cd40 * 0.25 * parz.rho * 1e-2;
    if( isIPT ){ ans += log( parz.IPThrP ); }   // scaling from IPT
    }
  } else {				     // HIV- :is this ever used?
    if( isIPT ){ ans += log( parz.IPThrN ); }   // scaling from IPT
  }
  return ans;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::b( parameters& parz ){
  double ans(0);
  if( hiv ){
    if( isART ){
      ans = parz.ARTf * parz.rho * 1e-2 * cd42 / (taumax - taustart);
    } else {
      ans = parz.cd40 * 0.75 * parz.rho * 1e-2/ taumax;
    }
  }
  return ans;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double person::HIVrelinfness( parameters& parz ){
  double ans(0);
  if( hiv ){
    ans = 1;
    if( tau < parz.durA ){ans = parz.riA;} else {
      if( (taumax - tau) < parz.durE ){ans = parz.riA;}
    }
    if( isART ){ans *= parz.ARTeps;}
  }
  ans *= parz.lam;		// overall scaling
  return ans;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//end person member funs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
