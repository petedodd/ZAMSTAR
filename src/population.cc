//#######################
#include "population.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::build_new_households( parameters& parz, int nob, const gsl_rng *r ){
  //creating households
  newhhcount = 0;
  double alph = 1.00;//elasticity, see newhhmoves.pdf
  newhhcount = gsl_ran_poisson( r, alph * hh.size() * nob / size );

  if ( newhhcount > 0 ){
    for ( int j = 0; j < newhhcount; ++j ){
      hh.push_back( household() );//add an empty household
      nohh++;//update no hhs
      ++hhNn.at(0);
      ++hhNnm.at(0).at(0);
    }
  }
  return newhhcount;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double population::update_hazards( ofstream& lgf, parameters& parz ){//hazard and foi updater
  //  int j = 0;//checking success
  //update FOI -- part of class!
  FOI = 0;
  rFOI = 0;
  // next 2 in case dynamic HIV model
  vector< double > hcatfoi(3,0.0); // FOI from cats
  vector< double > hcatpops(3,0.0); // sa cats

  //hazard update
  unsigned int ui = 0;//, uj;//second for complex haz
  unsigned int loc = 0;//hh location
  double sumF(0), sumrF(0), sumH(0), foihh(0), rfoihh(0);
  unsigned int count = 0;
  ////////////////////////////////////////////////
  if( !parz.AGEflag ){//no age cats
    ////////////////////////////

    cout <<"Haven't checked this code in update_hazards in ages!Bailing..."<<endl;
    exit(-1);

    for( ui = 0; ui < hh.size(); ++ui ){// set hfoi=0
      hh.at(ui).hfoi = 0;
      //hh.at(ui).hiv = 0;//hiv status of hh
    }//for nowt if hh off


    for( ui = 0; ui < v.size(); ++ui ){
      if ( v.at(ui).isD && !v.at(ui).isX ){//infectious
	loc = v.at(ui).hhid;//where
	if( loc > hh.size() ){cout << loc <<" bigger than "<< hh.size() << "in pop::update_haz...:bailing!"<<endl;exit(-1);}//safety	
	if( v.at(ui).isSmr ){//smr+
	  hh.at(loc).hfoi += parz.hhFLAG * parz.betaH;
	  sumF += parz.beta;
	} else {
	  sumF += parz.beta * parz.fSmr;
	  hh.at(loc).hfoi += parz.hhFLAG * parz.betaH * parz.fSmr;
	}
      }
      if ( v.at(ui).hiv && !v.at(ui).isART  && !v.at(ui).isX ) {hh.at(loc).hiv = 1;}
      //hiv status of hh - doesn't count if on ART(?)
    }//FOI done
    FOI = sumF / size;

  } else {//////////////////////////////////////AGEflag on!

    FOIvec.assign( parz.nbins, 0.0 );//zero the age-category things
    BINvec.assign( parz.nbins, 1e-6 );// (safely for denominator...)

    for( ui = 0; ui < hh.size(); ++ui ){// set hfoi=0
      hh.at(ui).hfoi = 0;
      hh.at(ui).hiv = 0;//hiv status of hh
      hh.at(ui).tb = 0;//tb status of hh
    }//for nowt if no hhs

    for( ui = 0; ui < v.size(); ++ui ){
      if( !v.at(ui).isX ){ ++BINvec.at( v.at(ui).acat ); }//one more in cat
      if ( v.at(ui).isD && !v.at(ui).isX ){//infected
	loc = v.at(ui).hhid;//where

	if( v.at(ui).isSmr ){//smr+
	  hh.at(loc).hfoi += parz.hhFLAG * parz.betaH;
	  FOIvec.at( v.at(ui).acat )  += parz.beta;
	} else {
	  FOIvec.at( v.at(ui).acat ) += parz.beta * parz.fSmr;
	  hh.at(loc).hfoi += parz.hhFLAG * parz.betaH * parz.fSmr;
	}
	// FOI is now the infectious weight in each acat
	// and BIN the number in each acat

      }	// end infected
      if ( v.at(ui).age >=18  && !v.at(ui).isX ) {
	if ( v.at(ui).hiv  ) { loc = v.at(ui).hhid; hh.at(loc).hiv += 1; }
	if ( v.at(ui).isD  ) { loc = v.at(ui).hhid; hh.at(loc).tb += 1;  }
      }
      //this is where hhhiv status gets set


    }//FOI done

    //N<-F/N
    for( unsigned int uu = 0; uu < parz.nbins; ++uu ){
      BINvec.at( uu ) = FOIvec.at(uu) / BINvec.at(uu); // BIN now the fraction of infection
    }
    for( unsigned int uu = 0; uu < parz.nbins; ++uu ){//matrix mult
      FOIvec.at(uu) = 0;
      for( unsigned int vv = 0; vv < parz.nbins; ++vv ){
	FOIvec.at(uu) += parz.agedata.at(uu).at(vv) * BINvec.at(vv);
      }
    }
    // this formulation implies that agedata_{uv} is the number of contacts w/v 
    // by a member of group u

    // ------ Foi for resitant strain
    if( parz.Rflag ){
      rFOIvec.assign( parz.nbins, 0.0 );//zero the age-category things
      BINvec.assign( parz.nbins, 0.0 ); // copy

      for( ui = 0; ui < hh.size(); ++ui ){// set hfoi=0
	hh.at(ui).rhfoi = 0;
      }//for nowt if no hhs

      for( ui = 0; ui < v.size(); ++ui ){
	if( !v.at(ui).isX ){ ++BINvec.at( v.at(ui).acat ); }//one more in cat
	if ( v.at(ui).isrD && !v.at(ui).isX ){//infected
	  loc = v.at(ui).hhid;//where
	  if( v.at(ui).isSmr ){//smr+
	    hh.at(loc).rhfoi += parz.hhFLAG * parz.betaH;
	    rFOIvec.at( v.at(ui).acat )  += parz.beta;
	  } else {
	    rFOIvec.at( v.at(ui).acat ) += parz.beta * parz.fSmr;
	    hh.at(loc).rhfoi += parz.hhFLAG * parz.betaH * parz.fSmr;
	  }
	}	// end infected
      }//rFOI stage 1 done

    //N<-F/N
      for( unsigned int uu = 0; uu < parz.nbins; ++uu ){
	BINvec.at( uu ) = rFOIvec.at(uu) / BINvec.at(uu); // BIN now the fraction of infection
      }
      for( unsigned int uu = 0; uu < parz.nbins; ++uu ){//matrix mult
	rFOIvec.at(uu) = 0;
	for( unsigned int vv = 0; vv < parz.nbins; ++vv ){
	  rFOIvec.at(uu) += parz.agedata.at(uu).at(vv) * BINvec.at(vv);
	}
      }

    } // rFOI complete

    // dynamic HIV
    if( parz.DHflag ){
      // see fraser/cori and refs
      for( unsigned int uu = 0; uu < v.size(); ++uu ){
	if( !v.at(uu).isX && v.at(uu).age > 15 ){
	  ++hcatpops.at( v.at(uu).scat-1 );
	  double infw = v.at(uu).HIVrelinfness( parz );
	  hcatfoi.at( v.at(uu).scat-1 ) += infw;
	}
      }	// total pops and FOIs from each cat

      double JR1(0), JR2(0), JR3(0);
      //double JA1(0), JA2(0), JA3(0); 
      // see fraser/cori and refs and correct!!
      double tot = hcatfoi.at(0)+hcatfoi.at(1)+hcatfoi.at(2); 
      JR3 = parz.sbc3 * tot;
      JR2 = parz.sbc2 * tot;
      JR1 = parz.sbc1 * tot;
      tot = JR1+JR2+JR3;
      double ptot = parz.sbc1*hcatpops.at(0)+parz.sbc2*hcatpops.at(1)+parz.sbc3*hcatpops.at(2); 
      double rbit = (1-parz.sbA) * tot / ptot;
      hcatfoi.at(0) = rbit + parz.sbA*parz.sba1 * JR1 / hcatpops.at(0); 
      hcatfoi.at(1) = rbit + parz.sbA*parz.sba2 * JR2 / hcatpops.at(1); 
      hcatfoi.at(2) = rbit + parz.sbA*parz.sba3 * JR3 / hcatpops.at(2); 

      hcatfoi.at(0) *= parz.sbc1 * parz.lam; // worry about lam 2x
      hcatfoi.at(1) *= parz.sbc2 * parz.lam;
      hcatfoi.at(2) *= parz.sbc3 * parz.lam;

    } // and DHflag

  }////////////////// end else for age flag

  //DEBUG!
  // lgf << "FOI=";
  // for(unsigned int vi = 0; vi < FOIvec.size();++vi ){lgf << FOIvec.at(vi)<<",";} lgf<<endl;


  double TBRATE(0);  // TB activations
  TBRATEN = 0;  TBRATEP = 0;	// these are stratified and part of pop
  double sumK(0); int kden(0);
  sumF = 0;
  sumH = 0;
  ghaz.assign( ghaz.size(), 0.0 );

  for( ui = 0; ui < v.size(); ++ui ){
    if ( !v.at(ui).isX ){//not dead
      // TB activation rate
      TBRATE += v.at(ui).tbrate;
      if( v.at(ui).hiv ){
	TBRATEP += v.at(ui).tbrate;
      } else {
	TBRATEN += v.at(ui).tbrate;
      }
      if ( parz.hhFLAG ){
	loc = v.at(ui).hhid;//which hh
	foihh =  hh.at( loc ).hfoi;//the household foi
	rfoihh = ( parz.Rflag>0 ? hh.at( loc ).rhfoi : 0.0 ); // rest
	if ( hh.at(loc).hiv  - v.at(ui).hiv > 0 ){ v.at(ui).hhhiv = 1;} else{ v.at(ui).hhhiv = 0;}
	if ( hh.at(loc).tb - v.at(ui).isD > 0 ){ v.at(ui).hhtb = 1;} else { v.at(ui).hhtb = 0; }
	//shares hh with someone hiv+ve, other than self
      } else { foihh = 0; }
      if( parz.AGEflag ){
	FOI = FOIvec.at( v.at(ui).acat ); // DEBUG!
	rFOI = ( parz.Rflag>0 ? rFOIvec.at( v.at(ui).acat ) : 0.0 ); // rest
	//FOI = 0.05;
      }//otherwise FOI as from above

      double HIVfoi = 1;
      if( parz.DHflag ){ HIVfoi = hcatfoi.at( v.at(ui).scat-1 );}
      double foiarg[] = { HIVfoi, FOI+foihh, 1.0, rFOI+rfoihh };
      v.at(ui).update_haz( foiarg, parz, time );//individual hazard

      //global hazard -- done directly because of individual TSIs
      sumF += FOI + foihh;//v.at(ui).indhaz.at(1);
      sumrF += rFOI + rfoihh;//
      sumH += foihh;
      ghaz.at(0) += v.at(ui).indhaz.at(0); // add on this
      ghaz.at(1) += v.at(ui).indhaz.at(1); // add on this
      ghaz.at(3) += v.at(ui).indhaz.at(3); // add on this
      // NB foi different from hazard of infection as HIV+ etc carry risks
      if( v.at(ui).age < 12 ){//kid mean
	++kden;
	sumK += FOI + foihh;
      }
    }
    ++count;
  }


  if( !count == v.size() ) {cerr << "error in hazard update" << endl;}
  ghaz.at(2) = TBRATE;		    // TB activations
  if( !parz.DHflag ){		    // HIV not dynamic
    if( !parz.hivmac ){
      ghaz.at(0) = parz.h1549( time ) * size1549;      // HIV infections
    } else {					       // HIVMAC version
      ghaz.at(0) = parz.h1549( time ) * N15;      // HIV infections
    }
  } 
  sumF /= size;//this is then the mean FOI experienced...
  sumH /= size;
  hhfoibar = sumH;//return by side-effect
  foibar = sumF;
  foik = (kden > 0 ? sumK / kden : -1);
  TBRATEN /= (size - hiv + 1e-6); // stratified HIV- tb rate
  TBRATEP /= (hiv + 1e-6);	  // stratified HIV+ tb rate

  if(!parz.qlog){
    lgf << "\t\tFOIbar=" << foibar  <<endl;
    lgf << "\t\trFOIbar=" << sumrF / size  <<endl;
    lgf << "\t\thhFOIbar=" << hhfoibar  <<endl;
    lgf << "\t\tFOIk=" << foik  <<endl;
    lgf << "\t\tmean TBRATEN(-)=" << TBRATEN  <<endl;
    lgf << "\t\tmean TBRATEP(+)=" << TBRATEP  <<endl;
    lgf << "\t\tHIV incidence=" << ghaz.at(0) / size1549  <<endl;
  }//mean FOI

  return sumF;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int population::update_ages( const gsl_rng *r, parameters& parz ){ //incrementer of times
  //  inlcudes times:
  time += parz.dt;
  parz.time = time;		// a copy of the time for some uses
  // ages
  int j = 0;//testing success
  unsigned int ui;
  unsigned int count = 0;
  for( ui = 0; ui < v.size(); ++ui ){
    if ( v.at(ui).isX != 1){//alive!
      v.at(ui).update_age( r, parz );
      if ( v.at(ui).age >= 15 && v.at(ui).is1549 == 0){//just aged in
	v.at(ui).is1549 = 1;
      }
      if ( v.at(ui).age > 50 && v.at(ui).is1549 == 1){//just aged out
	v.at(ui).is1549 = 0;
      }
    }
    ++count;
  }
  if( count == v.size() ){ j = 1;} else {cerr << "error in age updates (bail)" << endl; exit(-1);}
  return j;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::hhmovement( const gsl_rng *r, ofstream& lgf, int num, parameters& parz, bool record ){

  int j = 0;//for acceptance rate, etc

  //set weights
  vector < double > hhweights( hh.size(), 0.0 );//these are the f's at hh level
  for( unsigned int ui = 0; ui < hh.size(); ++ui ){
    unsigned int nopeo = hh.at(ui).members.size();//how many in house
    if( nopeo >= f_weight.size() ){
      hhweights.at(ui) = 0;
    } else {
      hhweights.at(ui) = f_weight.at( nopeo );
    }
  }
  
  int od(0), oo(0), odw(0), odm(0), oow(0), oom(0);
  bool  grown(false), grown2(false);
  
  while( num > 0 ){
    //pick a random
    int who = gsl_rng_uniform_int( r, v.size() );//who moving
    if( !v.at(who).isX ){//they're not dead
      int loc = find_one_event( r, hhweights );//where to
      
      try{
	oom = hh.at( v.at(who).hhid ).nmen;//origin
	oow = hh.at( v.at(who).hhid ).nwomen;
	oo =  oow + oom;
	odw = hh.at(loc).nwomen;//destination
	odm = hh.at(loc).nmen;
	od = odw + odm;
	if ( oo == 0 ){throw ;}
      }
      catch (exception& ex){
	cout << ex.what() <<" in picked person..." <<endl;
	lgf <<"who="<<who<<endl;
	lgf <<"hhid "<<v.at(who).hhid<<endl;
	lgf <<"age "<<v.at(who).age <<endl;
	lgf <<oom<<","<<oow<<endl;
	cout << "oo,od = " << oo <<"," << od <<endl;
	cout << "size: " << hhNn.size() << endl;
	exit(-1);
      }

      //move in
      --num;
      try{
	if ( v.at(who).gender == 0 ){//girl
	  hh.at( v.at(who).hhid ).remove( who, 0 );//take away from old
	  hh.at( loc ).add( who, 0 );//add to new place
	} else {//guy
	  hh.at( v.at(who).hhid ).remove( who, 1 );//take away from old
	  hh.at( loc ).add( who, 1 );//add to new place	  
	}
	v.at( who ).hhid = loc;//where are they now?
	
      }
      catch(exception& ex ){
	cout << ex.what() <<" in trying to move person"<<endl;
      }

      grown = false;
      try{
	//update hhsizes
	--hhNnm.at( odw ).at( odm );//out with the old
	--hhNnm.at( oow ).at( oom );
	if ( v.at(who).gender == 0 ){//female
	  if ( (odw + 1) >= (int) hhmaxn  ){// a new size of hh
	    ++hhmaxn;//grow max
	    grown = true;
	    vector< unsigned int > tempv; 
	    //	    vector< unsigned int > tempv; 
	    for ( unsigned int ui = 0; ui < hhNnm.at(0).size(); ++ui ){ tempv.push_back(0); }
	    hhNnm.push_back( tempv );//stick on end.
	  }
	  ++hhNnm.at( odw + 1 ).at( odm );//in with the new
	  ++hhNnm.at( oow - 1 ).at( oom );
	} else {//male
	  if ( (odm + 1) >= (int) hhmaxm  ){//a new size of hh
	    ++hhmaxm;//grow max
	    grown = true;
	    for ( unsigned int ui = 0; ui < hhNnm.size(); ++ui ){ hhNnm.at(ui).push_back(0); }
	  }
	  ++hhNnm.at( odw ).at( odm + 1 );// in with the new
	  ++hhNnm.at( oow ).at( oom - 1 );
	}//end hhupdate
      }
      catch( exception& ex ){
	cout << ex.what() <<" in updating hhNnm..."<< endl;
	cout << "odm,odw = " << odm << "," << odw << endl;
	cout << "oom,oow = " << oom << "," << oow << endl;
	cout << "gender: " << v.at( who ).gender << endl;
	cout << "grown?: " << grown << endl;
	cout << "hhmaxm, hhmaxn = " << hhmaxm << "," << hhmaxn << endl; 
	cout << "sizes = " << hhNnm[0].size() <<","<< hhNnm.size() << endl; 
	exit(-1);
      }

      try {
	grown2 = false;
	if( (od+1) > (int) hhmax ){
	  hhNn.push_back( 0 );
	  grown2 = true;
	  ++hhmax;
	} 
	//	if ( hhNn.at(od) == 0 ){lgf << "zero at d " << od <<endl;}
	--hhNn.at( od );
	if ( hhNn.at(oo) == 0 ){lgf << "zero at o " << oo <<endl;}
	--hhNn.at( oo );
	++hhNn.at( (od+1) );
	++hhNn.at( (oo-1) );
      }
      catch( exception& ex){
	cout << ex.what() << " in updating hhNn..."<< endl;
	cout << grown2 << endl;
	cout << od+1 <<","<<hhmax<<","<<hhNn.size()<<endl;
      }

    }//end not dead !isX
  }//end while

  return j;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::counts( ofstream& lgf, parameters& parz ){//keeps tabs of population counts
  //this may not be quicker, but sure is simpler
  //hiv,1549 and size are stil counted at changes....
  int Sn(0), Ln(0), Dn(0), Tn(0), Sp(0), Lp(0), Dp(0), Tp(0);
  sizeyoung = 0;
  size1549 = 0;
  hiv1549 = 0;  hiv18 = 0; hiv15 = 0; noart15 = 0;
  nophhi = 0;nohhipt = 0;nohhart = 0; notecf = 0; nochhi = 0;
  size18 = 0; prev18 = 0; prevSP18 = 0; prevT18 = 0; prevTH18=0;        // adult population
  rL = 0;		       // resistant
  ARTnum=0; 			// ART number
  ARTdenom=0;			// ART eligible
  shivcount = 0;		// self-knowledge of HIV, including ART
  int nbreaks = 20;		// NB this is set in TBheader too
  cd4hist.assign(nbreaks, 0.0);	// CD4 histogram
  iptnos = 0; iptnosP = 0; iptnosA = 0; cipth = 0; // IPT counts
  tbcoprev = 0; hivcoprev = 0;		// HH co-prevs
  int tbcoprevdenom(0), hivcoprevdenom(0);
  // HIV-MAC counters
  GLRn = 0; GLRp = 0; prev15 = 0;
  GLRpo350 = 0; GL15 = 0;
  N15 = 0; h15u200 = 0; h15u200an = 0; N15n = 0;
  GLnoteli = 0; prevART = 0; prevu200 = 0; prevo200 = 0; Noteli = 0;
  prev15an1 = 0; prev15an2 = 0; prev15an3 = 0;
  GLpnoART = 0;
  GLRnCD41 = 0; GLRnCD42 = 0; GLRnCD43 = 0; // CD4 stratified counts
  GLRpCD41 = 0; GLRpCD42 = 0; GLRpCD43 = 0;
  // end HIVMAC

  // loop through population
  for( unsigned int uu = 0; uu < v.size(); ++uu){
    if( !v.at(uu).isX ){	// not dead
      if( v.at(uu).isrL ){ ++rL; }
      if( v.at(uu).phhi ){ ++nophhi; }
      if( v.at(uu).hhi ){ ++nochhi; }
      if( v.at(uu).hhipt ){ ++nohhipt; }
      if( v.at(uu).tecf ){ ++notecf; }
      if( v.at(uu).iptd && v.at(uu).hiv ){ ++cipth; }
      if( v.at(uu).hhart && v.at(uu).isART ){ ++nohhart; }//includes indirect
      if( v.at(uu).age < 12 && v.at(uu).age > 5 ){ ++sizeyoung; }
      if( v.at(uu).age < 50 && v.at(uu).age >= 15 ){ ++size1549; }
      if( v.at(uu).isIPT ){	// IPT stuff
	++iptnos;
	if( v.at(uu).hiv ){ ++iptnosP; }
	if( v.at(uu).isART ){ ++iptnosA; }
      }
      if( v.at(uu).age >=15 ){	// HIVMAC
	++N15;
	if( !v.at(uu).hiv ){++N15n;}

	if( v.at(uu).isD ){ 
	  ++prev15; 
	  if(!v.at(uu).hiv){++prev15n;}
	  if(v.at(uu).isART ){++prevART;}
	  if(v.at(uu).hiv && !v.at(uu).isART ){
	    if( v.at(uu).cd4 < 200 ){ // <200
	      ++prevu200; ++prev15an1;
	    } else{
	      if( v.at(uu).cd4 < 350 ){ // 200-350
		++prev15an2;
	      } else{		// >350
		++prev15an3;
	      }
	      ++prevo200;
	    }
	  } // hiv not art 200 decision
	}
	double guide = (parz.time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
	if( !v.at(uu).hiv || v.at(uu).cd4 > guide ){
	  ++Noteli;		// frac denominator
	}

	// prevalences guideline readiness
	if( v.at(uu).GL == 1 ){
	  ++GL15;
	  double guide = (parz.time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
	  if( !v.at(uu).hiv || v.at(uu).cd4 > guide ){
	    ++GLnoteli;
	  }
	  if( !v.at(uu).hiv ){ ++GLRn;} else {
	    ++GLRp;
	    if( v.at(uu).cd4 > 350 ){
	      ++GLRpo350;
	    }
	  }
	}		    // end GL==1
	if( v.at(uu).hiv ){	// all over 15
	  ++hiv15;
	  if( !v.at(uu).isART ){
	    ++GLpnoART;		// GLR HIV+ ART-
	    if( v.at(uu).cd4 < 200 ){
	      ++h15u200; ++hpop1;
	      if(v.at(uu).GL == 1){++GLRpCD41;} else {++GLRnCD41;}
	    } else{
	      if( v.at(uu).cd4 < 350 ){ // 200-350
		++hpop2;
		if(v.at(uu).GL == 1){++GLRpCD42;} else {++GLRnCD42;}
	      } else{		// >350
		++hpop3;
		if(v.at(uu).GL == 1){++GLRpCD43;} else {++GLRnCD43;}
	      }
	      ++h15u200an;
	    }
	  }
	  // NB this has been changed to match the new definition
	  // i.e. summing to HIV + not on ART
	} // end HIV+ and over 15
      }	  // end over 15
      if( v.at(uu).age >=18 ){	// adult stuff
	++size18;
	if( v.at(uu).hhtb ){ ++tbcoprevdenom; }
	if( v.at(uu).hhhiv ){ ++hivcoprevdenom; }
	if( v.at(uu).isD ){ 
	  ++prev18;
	  if( v.at(uu).isSmr ){ ++prevSP18; }
	  if( v.at(uu).hhtb ){ ++tbcoprev; }
	}
	if( v.at(uu).isT ){
	  ++prevT18;
	  if( v.at(uu).hiv ){ ++prevTH18;}
	}
      }
      if ( v.at(uu).hiv ){	// HIV stuff
	if( v.at(uu).shiv ){ ++shivcount; }
	if( v.at(uu).age >=18 ){ 
	  ++hiv18; 
	  if( v.at(uu).hhhiv ){ ++hivcoprev; }
	}
	if( v.at(uu).age < 50 && v.at(uu).age >= 15 ){ 
	  ++hiv1549; 
	  //	  if( v.at(uu).isART ){	// checking on ART only
	  // cd4 stuff
	  int cd4bin = 0;
	  for( int k = 0; k < nbreaks; ++k ){ // establish bin
	    if( v.at(uu).cd4  > k * parz.cd40 / nbreaks ){ cd4bin = k; }
	  }
	  ++cd4hist.at(cd4bin);
				//}//ART only
	}
	if( v.at(uu).isS ){ ++Sp;}
	if( v.at(uu).isD ){ ++Dp;}
	if( v.at(uu).isT ){ ++Tp;}
	if( v.at(uu).isL ){ ++Lp;}
	if( v.at(uu).isART ){ //on ART
	  ++ARTnum; 
	} else {
	  double guide = (time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
	  if( v.at(uu).cd4 < guide ){
	    ++ARTdenom;//eligible ART not on ART
	  }
	}
      } else {
	if( v.at(uu).isS ){ ++Sn;}
	if( v.at(uu).isD ){ ++Dn;}
	if( v.at(uu).isT ){ ++Tn;}
	if( v.at(uu).isL ){ ++Ln;}
      }//hiv
    }//X
  }//popn

  // HH clustering prevalences
  hivcoprev /= (hivcoprevdenom + .1);
  tbcoprev /= (tbcoprevdenom + .1);



  // cd4 stuff
  if( hiv1549 == 0 ){
    cd4hist.at(nbreaks-1) = 1; 
  } else {			// normalise
    for( int k = 0; k < nbreaks; ++k ){ cd4hist.at(k) /= hiv1549; }
  }

  popsP.at(0) = Sp;
  popsP.at(1) = Lp;
  popsP.at(2) = Dp;
  popsP.at(3) = Tp;
  popsN.at(0) = Sn;
  popsN.at(1) = Ln;
  popsN.at(2) = Dn;
  popsN.at(3) = Tn;
  //safety
  hiv = popsP.at(0) + popsP.at(1) + popsP.at(2) + popsP.at(3);
  size = popsN.at(0) + popsN.at(1) + popsN.at(2) + popsN.at(3) + hiv;


  // analysis for activation counts...zerod in scheduledevents,
  // counted in event
  if(noa>0){		// for recording activations
    propFI = 1.0 - (double)slow / noa;
    propIH = (double)noah / noa; 
    propViaF = (double)noaf / noa; 
  }
  if(noah>0){
    propFIP = 1.0 - (double)noahs / noah;
  }
  // and activation logging
  if(!parz.qlog && parz.time > parz.starttime ){
    lgf <<"\t "<<noa <<"("<<noah<<"+)"<<" activations, "<<noc<<" Rx completions"<<endl;
    lgf <<"\t "<< 100*propViaF <<" % activations via Fast route"<<endl;

    lgf << "\tav sigma(act)= "<<av2/(noa+1e-10)<<" / "<< noa <<endl;
    lgf << "\tav sigma(act+) = "<<av3/(noah+1e-10)<<" / "<< noah <<endl;
    lgf << "\tav tau prop (act+)= "<<av4/(noah+1e-10)<<" / "<< noah <<endl;
    lgf << "\tav meancd4 (act+)= "<<av5/(noah+1e-10)<<" / "<< noah <<endl;

    lgf << "\tSLOW = " << slow << endl;
  }

  //  cout <<"ARTeliinHIV=" <<(double)ARTdenom / (hiv + 1) << endl;
  // cout <<hiv<<", "<< ARTdenom <<", "<<ARTnum  << endl;

  return Sn;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int population::hhcounts( ofstream& lgf, parameters& parz ){//keeps tabs of population counts
  hhmaxn = 15; hhmaxm = 15; hhmax = 0;//last is actual max size;
  hhNn.assign( 15, 0 );//0 houses of size
  hhNnm.assign( 15, vector< unsigned int >(15,0) );//0 houses of size
  unsigned int coh(0);
  for( unsigned int ui = 0; ui < hh.size(); ++ui ){//hh loop
    //resizing stuff//////
    coh = (hh.at(ui).nwomen + hh.at(ui).nmen);
    try{
    while( coh > (hhNn.size()-1) ){//resize
      hhNn.push_back(0);
    }//end resize
    }
    catch( exception& ex ){
      cout << ex.what() << " in resizing hhNn!"<<endl;
      exit(-1);
    }
    try{
    while( hh.at(ui).nwomen >= hhmaxn ){//max women exceeded: extra row required
      ++hhmaxn;//grow max
      vector< unsigned int > tempv; 
      for ( unsigned int vi = 0; vi < hhNnm.at(0).size(); ++vi ){ 
	tempv.push_back(0); 
      }
      hhNnm.push_back( tempv );//extra row
    }//end n-resize
    }
    catch( exception& ex ){
      cout << ex.what() << " in resizing women!"<<endl;
      exit(-1);
    }
    try{
    while( hh.at(ui).nmen >= hhmaxm ){//max men exceeded: extra col required
      ++hhmaxm;//grow max
      for ( unsigned int vi = 0; vi < hhNnm.size(); ++vi ){ 
	hhNnm.at(vi).push_back(0);
      }
    }//end m-resize
    }
    catch( exception& ex ){
      cout << ex.what() << " in resizing men!"<<endl;
      exit(-1);
    }
    //resizing stuff//////

    //update counts
    try{
    ++hhNn.at( coh );
    ++hhNnm.at( hh.at(ui).nwomen ).at( hh.at(ui).nmen );
    if( coh > hhmax ){ 
      hhmax = coh;//update biggest hh size
    }
    }
    catch( exception& ex ){
      cout << ex.what() << " in updating counts!!"<<endl;
      cout << coh <<", " <<hhNn.size() << endl;
      cout << hh.at(ui).nwomen <<", "<<hh.at(ui).nmen << endl;
      cout << hhNnm.size() <<", "<<hhNnm.at(0).size() << endl;
      cout << hhmaxn <<", "<<hhmaxm << endl;
      exit(-1);
    }
  }//end hh loop

  return( (int)hhmax );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::HHpostinit( const gsl_rng *r, ofstream& lgf, parameters& parz ){//initializing properly
  nohh = hh.size();		// the number of households
  unsigned int ii(0);
  //recording max hh sizes
  hhmaxn = 0;
  hhmaxm = 0;
  hhmax = 0;
  for( ii = 0; ii < hh.size(); ++ii ){//max sizes
    if ( hh.at(ii).nwomen > hhmaxn ){ hhmaxn = hh.at(ii).nwomen; }
    if ( hh.at(ii).nmen > hhmaxm ){ hhmaxm = hh.at(ii).nmen; }
    if ( hh.at(ii).nmen + hh.at(ii).nwomen > hhmax ){  
      hhmax = hh.at(ii).nmen + hh.at(ii).nwomen;
    }
  }
  lgf <<"hhmaxn="<<hhmaxn<<",hhmaxm="<<hhmaxm<<",hhmax="<<hhmax<<endl;

  hhNnm.clear(); hhNn.clear();
  //setting up hhNnm -- matrix of sizes
  vector< unsigned int > tempv;
  for ( ii = 0; ii < (hhmaxm + 1); ++ii){//NB this is so the indices are literal nos
    tempv.push_back(0);
  }
  for( ii = 0; ii < (hhmaxn + 1); ++ii ){//same NB
    hhNnm.push_back( tempv );
  }
  //filling it up
  for( ii = 0; ii < hh.size(); ++ii ){
    ++hhNnm.at( hh.at(ii).nwomen ).at( hh.at(ii).nmen );
  }
  //totals
  for ( ii = 0; ii < (hhmax + 1); ++ii ){
    hhNn.push_back(0);
  }
  //filling it up
  for( ii = 0; ii < hh.size(); ++ii ){
    ++hhNn.at( hh.at(ii).nwomen + hh.at(ii).nmen );
  }
  vector< vector< double > > hhsizes( hhNnm.size(), vector< double >(hhNnm[0].size(),0.0) );
  //matrix of hh sizes - initialized to zero
  //copy the count version for below working
  for( unsigned int ii = 0; ii < hhsizes.size(); ++ii ){
    for (unsigned int jj = 0; jj < hhsizes.at(ii).size(); ++jj){
      //lgf << " - "<<ii<<","<<jj<<endl;
      hhsizes.at(ii).at(jj) = hhNnm.at(ii).at(jj) ;
    }
  }

  vector< double > hhtotsizes;//vector of same but only total hh sizes
  hhtotsizes.push_back(0);
  //renormalize
  int nohh0 = 0;
  for( unsigned int ii = 0; ii < hhsizes.size(); ++ii ){
    for (unsigned int jj = 0; jj < hhsizes.at(ii).size(); ++jj){
      //nohh0 += hhsizes.at(ii).at(jj);
      if( (ii+jj) != 0 ){//ignore 0	
	nohh0 += hhsizes.at(ii).at(jj);
	if ( ii+jj >= hhtotsizes.size() ){
	  hhtotsizes.push_back( hhsizes.at(ii).at(jj) );
	} else{
	  hhtotsizes.at( ii + jj ) += hhsizes.at(ii).at(jj);
	}
      }//ignore the zeros
    }
  }
  //lgf << nohh0 << " initial households..." << endl;
  for( unsigned int ii = 0; ii < hhsizes.size(); ++ii ){
    for (unsigned int jj = 0; jj < hhsizes.at(ii).size(); ++jj){
     hhsizes.at(ii).at(jj) /= nohh0 ;
    }
  }
  hhsizes.at(0).at(0) = 0;
  for( unsigned int ii = 0; ii < hhtotsizes.size(); ++ii ){hhtotsizes.at(ii) /= nohh0;}
  hhtotsizes.at(0) = 0;

  //compute weights for alternative method...
  //going to imagine buffer of 10% empty households
  //NB the above are proportions...
  f_weight.assign( hhtotsizes.size(), 0.0 );
  a_weight.assign( hhtotsizes.size(), 0.0 );
  f_weight.at(0) = hhtotsizes.at(1);// / 0.1 ;//buffer
  double temp = 0.0; //buffer
  a_weight.at(0) = 1 - temp;//buffer
  for( unsigned int  ii = 1; ii < hhtotsizes.size()  - 1; ++ii ){
    f_weight.at(ii) =  (ii+1) * hhtotsizes.at(ii+1) / ( hhtotsizes.at(ii) + 1e-10 );
    if ( hhtotsizes.at(ii) < 1e-10 ){ f_weight.at(ii) = 0; }
    temp += (1 - 0.0) * hhtotsizes.at(ii);//buffer
    a_weight.at(ii) = ( 1 - temp )  /  ( hhtotsizes.at(ii) + 1e-10 );
    if ( hhtotsizes.at(ii) < 1e-10 ){ a_weight.at(ii) = 0; }
  }
  f_weight.at( hhtotsizes.size()-1 ) = 0;
  a_weight.at( hhtotsizes.size()-1 ) = 0;

  //  lgf << "N = ";
  // for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
  //   lgf << hhtotsizes.at(ii) <<",";
  // }
  // lgf << endl;

  // lgf << "F = ";
  // for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
  //   lgf << f_weight.at(ii) <<",";
  // }
  // lgf << endl;

  // lgf << "A = ";
  // for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
  //   lgf << a_weight.at(ii) <<",";
  // }
  // lgf << endl;

  lgf <<"\tdone!" << endl;

  return(nohh);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::snapshot( ofstream& lgf, vector< vector< double > >& ibmstate ){
  // to go through and recrod the whole population
  // age, acat, gender, TBu, TBRx, smrp,HIV, ART, HHHIV,phhi,tecf,
  // dead,

  // stuff on mortality causes?
  int artdefdead(0), artdef(0);
  for( unsigned int ui = 0; ui < v.size(); ++ui ){
    if( v.at(ui).artdefr ){
      ++artdef;
      if( v.at(ui).isX ){++artdefdead;}
    }
  }

  lgf << "mortality among ART defaulters (all time) = " << 
    ( 100.00 * artdefdead ) / (artdef+1) << endl;


  // alternate way - going through households
  // age,acat,gender,TBu,TBRx,smrp,HIV,ART,HHHIV,phhi,tecf,LTBI,hhsize,hhid,cd4,isX,artdefr
  vector< double > tmp( 22, 0.0 );
  int i(0);
  lgf << "starting population snapshot: size=" << v.size()<<"..."<< endl;
  lgf << "\t\t: size hh=" << hh.size()<<":"<< endl;
  for( unsigned int ui = 0; ui < hh.size(); ++ui ){
    unsigned int hhsize = (hh.at(ui).nmen + hh.at(ui).nwomen);
    if( hhsize > 0 ){
      list< unsigned int >::iterator it; // go through list of members
      for( it = hh.at(ui).members.begin(); it != hh.at(ui).members.end(); ++it ){
	unsigned int who = *it;
	//if( true ){	//including dead
	if( !v.at(who).isX ){	// not dead
	  tmp.at(0) = (double)v.at(who).age ; // record age
	  tmp.at(1) = (double)v.at(who).acat ; // record agecat
	  tmp.at(2) = (double)v.at(who).gender ; // record sex
	  tmp.at(3) = (double)v.at(who).isD ; // record TBu
	  tmp.at(4) = (double)v.at(who).isT ; // record TBRx
	  tmp.at(5) = (double)v.at(who).isSmr ; // record smrp
	  tmp.at(6) = (double)v.at(who).hiv ; // record HIV
	  tmp.at(7) = (double)v.at(who).isART ; // record ART
	  tmp.at(8) = (double)v.at(who).hhhiv ; // record HHHIV
	  tmp.at(9) = (double)v.at(who).phhi ; // record phhi
	  tmp.at(10) = (double)v.at(who).tecf ; // record tecf
	  tmp.at(11) = (double)v.at(who).isL ; // record LTBI

	  tmp.at(12) = (double)hhsize ; // hh size
	  tmp.at(13) = (double)ui ; // hh id

	  tmp.at(14) = v.at(who).cd4; // cd4 count

	  tmp.at(15) = v.at(who).isX; // dead?
	  tmp.at(16) = v.at(who).artdefr; // ART defaulter

	  tmp.at(17) = v.at(who).isIPT; // on IPT?
	  tmp.at(18) = v.at(who).iptd; // ever on IPT?
	  tmp.at(19) = v.at(who).dTBact; // protection?
	  tmp.at(20) = v.at(who).isPI; // previous infection?
	  tmp.at(21) = v.at(who).tbrate; // rate of TB

	  ibmstate.push_back( tmp );	// record these states
	  ++i;

	} // not dead
      }	  // household members
    }	  // not empty hh
  }	  // all households
  lgf << "...finished population snapshot:"<<i<<" alive!"<<endl;
  return i;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int population::setup_weights(ofstream& lgf, parameters& parz, int kappa, map<int,int>& ghk, vector< double >& weights, vector< unsigned int>& where, unsigned int& noe, unsigned int& num, const gsl_rng *r ){
  //this sets up the weights for an event kappa from amongst the eligibles
  // it has been change to return an array with noe, num and so
  // now includes event-specific calculations of howmany_events

  unsigned int k(0), j(0);//noe, num
  list< unsigned int >::iterator li;
  // keeping size pretty constant
  weights.assign( v.size(), 0 );//zero the weights
  where.assign( v.size(), 0 );//zero the wheres 

  switch( kappa ){
  case 3:
    {

      double omr = parz.mr( time );//for geting out migration rate
      omr = ( omr < 0 ? -omr : 0 ); // only outmigration
      j = howmany_events( r, omr * size, parz.dt  );//how many 

      if( j > 0 ){
	for ( unsigned int uu = 0; uu < v.size(); ++uu ){
	  if (  !v.at(uu).isX && (v.at(uu).age > 18) ){
	    //living, some chance of event
	    where.at(k) = uu;
	    weights.at(k) = 1;	// random
	    ++k;//advance abbreviated counter
	  }
	}
      } // j>0

    }
    break;
    // Global hazard stuff
  case 4:			// resistant strain infection
    if( !parz.Rflag ){break;}
  case 1:			// Mtb. infection
  case 2:			// activation
  case -1:			// HIV infection
    {				// setting global hazards
      for ( unsigned int uu = 0; uu < v.size(); ++uu ){
	double value = v.at( uu ).indhaz.at(ghk[kappa]) - 1e-20; 
	if (  !v.at(uu).isX && value > 0 ){//living, some chance of event
	  where.at(k) = uu;
	  weights.at(k) = value;
	  ++k;//advance abbreviated counter
	}
      }

      j = howmany_events( r, ghaz.at(ghk[kappa]), parz.dt  );//how many 

      // j = 0;			// debug
      
    }
    break;
  case -2:			// Birth weights
  case -5:			// imr weights
    {// NB no fertility involved
      
      // choice here
      if( kappa == -2 ){				 // BR
	j = howmany_events( r, parz.br( time ) * size, parz.dt  );
	nob = j;					// for building hhs
      } else {					// OMR
        // including population-undershoot correction here
        relpop = (double)size / (size0+.1);
	double over = relpop - parz.rpop();
        over = ( over<-0.0005 ? -over  : 0);
        j = auxFloor( over * size );
        // imr
	double imr = parz.mr(time);
	imr = ( imr > 0 ? imr : 0 );
	j += howmany_events( r, imr, parz.dt  );//how many 
      }

      // j = 0;			// debug!!

      if( j > 0){	
	// from inside destination_hhs - dealing with weights over hhs
	if( parz.hhFLAG ){
	  for( unsigned int ui = 0; ui < hh.size(); ++ui ){
	    unsigned int nopeo = hh.at(ui).members.size();//how many in house
	    if( nopeo < a_weight.size() ){
	      weights.at(k) = a_weight.at( nopeo );
	      where.at(k) = ui;
	      ++k;			// advance abbreviated counter
	    }
	  }
	} else {			// no households
	  cout <<"Uh-oh! no households not allowed in setup_weights!"<<endl;
	  exit(-1);
	}
      }
      
    }
    break;

  case 11:			// correxion death
    {				// think about whether this happens after sch
      if( parz.rpopFN != "none" ){ // corx
	// relative population
	relpop = (double)size / (size0+.1);
	double over = relpop - parz.rpop();
	if( over > 0.0005 ){
	  j = auxFloor(  over * size ); // how many extra deaths needed
	} // else {j = 0;lgf <<"warning: population undershoot - correction not possible" << endl;
	  // lgf << size <<" "<<parz.popsize0<<" "<<relpop <<" "<<parz.rpop() << endl;}

	if( j > 0 ){
	  for ( unsigned int uu = 0; uu < v.size(); ++uu ){
	    if ( !v.at(uu).isX ){
	      where.at(k) = uu;
	      weights.at(k) = 1.0;//random
	      ++k;//advance abbreviated counter
	    }
	  }//population loop
	}  // j > 0

      }	// end flag check
      // NB need to be very careful comparing deaths across scenarios!
    }
    break;
    
  case -3:	//true ECF
    {
      if( parz.tECFflag && time > parz.tECFst && time < parz.tECFet ){
	
	for ( unsigned int uu = 0; uu < v.size(); ++uu ){
	  // if ( !v.at(uu).tecf && !v.at(uu).isX ){//a living ECF naive
          if ( !v.at(uu).tecf && !v.at(uu).phhi  &&!v.at(uu).isX ){//a living ECF + hh naive
	    where.at(k) = uu;
	    // weights.at(k) = parz.tECFhaz;//doesn't really matter: all same if!=0
            weights.at(k) = parz.tECFhaz * ((double)nochhi / size);//now infectious process
	    ++k;//advance abbreviated counter
	  }
	}//population loop
	// changed here too!
        j = howmany_events( r, k * parz.tECFhaz * ((double)nophhi / size), parz.dt  );//how many 
      } else {j=0;}
    }
    break;

  case -4:
    {
      // ART cessation from default
      for ( unsigned int uu = 0; uu < v.size(); ++uu ){
	if ( v.at(uu).isART && !v.at(uu).isX ){//a living on ART
	  where.at(k) = uu;
	  weights.at(k) = parz.artDR;//doesn't really matter: all same if!=0
	  ++k;//advance abbreviated counter
	}
      }//population loop
      j =  howmany_events( r, parz.artDR * ARTnum, parz.dt  );//how many 
    }
    break;


  default:
    cout << "unrecognized event in setup_weights!" << endl;
    exit(-1);
    break;
  }	      // end switch
  
  noe = k;			// noe=number eligible
  num = j;			// num=number of events
  return k;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int population::scheduledevents( ofstream& lgf, const gsl_rng *r, bool qlog, parameters& parz, map < string, int >& evk ){//for things going on in events
  //need r, qlog, and maybe parz
  //set a load of relevant counters
  // NB if you wanted to use the event key, you'd need the loops tother way
  // to avoid multiple look-ups

  noa=0; 			// no activations
  nob=0; 			// no births
  noaf=0; 			// no activations via fast.
  noafp=0; 			// no activations via fast. hiv+
  noafn=0; 			// no activations via fast. hiv-
  noc=0; 			// no Rx completions (succ or unsucc)
  nnd=0;			// no latently infected -- needed?
  noart=0;			// no commencing ART
  noah=0;			// no incident TB hiv+
  noahs=0;			// no incident TB hiv+ slow
  noaa=0;			// no incident TB art+
  noap=0;			// no incident TB smr+
  ntbd=0; 			// no TB deaths from Rx
  ntbdet=0; 			// no TB detections
  ntbsc=0; 			// no TB self-cures
  ntbdh=0; 			// no TB deaths hiv+
  ntbdeth=0;			// no TB detection hiv+
  ntbdetsp=0; 			// no TB detections smr+
  ntbdetf=0; 			// no tb detection "fast"
  kidi=0;adi=0;			// no u12 infected (not nec for first time), same for adults
  av2=0;		// computing mean tb tsi both L,(+) and incident TB
  av3=0; av4=0; av5=0;	       	// computing mean tb tsi for HIV, and HIV tsi
  prevhh=0.0; 			// mean prevalence TB in HH detections
  tbavd=0; 			// for mean time to TB outcome in detections
  tbav=0;			// for mean time to outcome in detect-avoids
  slow = 0;//debug
  propFI=-1; 			// proportion fast in incident
  propFIP=-1; 			// proportion fast in incident HIV+
  propViaF=-1; 			// proportion via fast route
  propFD=-1; 			// proportion fast in detected
  propIH=-1; 			// proportion HIV+ in incident
  propDH=-1;  			// proportion HIV+ in detected
  propDS=-1;			// proportion smr+ in detected
  foibar = 0; 			// mean FOI experienced in pop
  foik = 0;			// mean FOI in kids
  rfoi = 0;			// mean rFOI
  hhfoibar = 0; 		// mean FOI in the HH
  prevHHH=0;			// mean prevalence HIV in HH
  norx = 0;			// for mean time on Rx
  taurx = 0;			// for mean time on Rx
  dethhno = 0;			// the number in HH of dets (regardless HH)
  dethhno1549 = 0;		// the number in HH of dets (regardless HH) 1549
  nodtecf = 0;			// the number in detections who are tecf
  propECF = -1;			// proportion via ECF
  ntbdetprev = 0;		// no TB detections with previous TB
  iptends = 0;		     	// IPT completions

  int bbsstarts = 0;		// starts via bbs for log
  nod = 0; bdcount = 0; avAdeath = 0;
  double tau = 0;
  int j(0);
  //nod and tauav used for HIV deaths
  double test(0); 

  //  int seauton(0);
 //loop through events -/TB activations, Rx completions, and ART starts, HIV self-knows
  for ( unsigned int ui = 0; ui < v.size(); ++ui ){
    if ( !v.at(ui).isX ){

      //hiv death
      if ( v.at(ui).hiv  ){//infected, alive
      	if ( v.at(ui).tau > v.at(ui).taumax ){//must die!
      	  //cout << time << "," << v.at(ui).age << "," <<v.at(ui).tau <<"," <<v.at(ui).taumax <<","<<v.at(ui).isART << endl;
      	  tau += v.at(ui).taumax;
      	  j = event( ui, "death", evk, parz, r );//die!
      	  //j = event( ui, 6, parz, r );//die!
      	  nod++;
      	  if ( j == 0 ){ lgf << "**error enacting HIV death!**" <<  endl;}
      	}
      }

      // background death
      if(  (v.at(ui).age > v.at(ui).LE) ){ // time is spent
      	avAdeath += v.at(ui).age;
      	j = event( ui, "death", evk, parz, r );//die!
      	//j=event( ui, 6, parz, r );//die!
      	++bdcount;
      	if ( j == 0 ){ lgf << "**error enacting background death!**" <<  endl;}
      }

      // IPT completions
      if( v.at(ui).isIPT ){
	if( ( v.at(ui).hiv && v.at(ui).tonipt > parz.IPTdurnP ) ||
	    ( !v.at(ui).hiv && v.at(ui).tonipt > parz.IPTdurnN )){ // finish IPT
	  event( ui, -5, parz, r );//ipt end
	  ++iptends;
	  // no counting as yet
	}
      }

      // completing treatment
      if ( v.at(ui).isT ){++norx; taurx += v.at(ui).Rxtime;}
      if ( v.at(ui).isT && (v.at(ui).Rxtime > parz.RxL) ){
	double test = parz.Rxpfn(); double test2 = parz.Rxpd;
	if( v.at(ui).hiv ){ test = parz.Rxppfn(); test2 = parz.Rxppd; }
	if( gsl_rng_uniform(r) < test ){
	  event( ui, 10, parz, r );//completion-successful
	} else {		  // not successful
	  if ( gsl_rng_uniform(r) <  test2 ){ // unlucky
	    //die
	    if( v.at(ui).hiv ){++ntbdh;}
	    ++ntbd;
	    event( ui, 6, parz, r );//completion-unsucc-die
	  } else {		// lucky
	    event( ui, 7, parz, r );//treat like new infection - higher risk for recurrence
	  }
	}
	++noc;
      }

      // ART either routine or via household
      if( parz.ARTflag ){					 // is ART happening?

	if(!parz.hivmac){	// non-HIVMAC version of ART starts

	  // learning HIV status & linking to care
	  // sets threshold to guideline
	  if( !v.at(ui).shiv && v.at(ui).dshiv 
	      && (v.at(ui).artdefr == 0) ){ // ignorant but destined to learn
	    if( v.at(ui).cd4 <= v.at(ui).cd4st + parz.shivt ){ 
	      v.at(ui).shiv = 1; 	// linked
	      // those who are linked to care start at the guideline
	      double guide = (time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
	      v.at(ui).cd4st = guide; 	// guideline
	      if( parz.IPTshiv && gsl_rng_uniform(r) < parz.IPTcovPf() &&
		  time > parz.IPTst && !(!parz.IPTmultiple * v.at(ui).iptd ) && 
		  !v.at(ui).isT && !v.at(ui).isIPT && !v.at(ui).isD){
		event( ui, -4, parz, r );//put on IPT
	      }
	    } // learn
	  }
	  
	  // label those who were below CD4 thresholds before ART roll-out
	  if ( time >= parz.ARTst && time < (parz.ARTst + parz.dt) ){ 
	    // just passed start
	    if( v.at(ui).hiv && ( v.at(ui).cd4 < v.at(ui).cd4st )
		&& !v.at(ui).isART  && (v.at(ui).artdefr == 0)){ 
	      v.at(ui).bbs = 1;	// label as below-before-start
	    }
	  }
	  // apply rate to bbs's
	  if( v.at(ui).hiv && v.at(ui).bbs && (parz.artCU > 0)
	      && !v.at(ui).isART  && (v.at(ui).artdefr == 0)){ 
	    if( gsl_rng_uniform(r) < (1-exp(-parz.dt * parz.artCU)) ){
	      event( ui, -2, parz, r );//ART commencement
	      ++noart; ++bbsstarts;
	    }
	  }

	  // have now included second bit to get ART 2 linkers who've passed?
	  if( !v.at(ui).dart && v.at(ui).hiv && !v.at(ui).isART 
	      && (v.at(ui).artdefr == 0) && (time > parz.ARTst) ){ 
	    // dart now recording
	    // that no thresholds passed
	    if( v.at(ui).cd4 <= v.at(ui).cd4st ){ 
	      // HH now dealt with as shift at detection!
	      v.at(ui).dart = 1;	// threshold passed -- CHANGE
	      // criterion for starting ART
	      if( ( ((gsl_rng_uniform(r) < ( 1- parz.shivp )*parz.ARTcov(time) )
		     && !v.at(ui).shiv ) || v.at(ui).shiv ) 
		  //		&& !v.at(ui).ARTrefusenik
		  ){
		if( ( parz.IPTshiv || parz.IPTart ) && 
		    (time > parz.IPTst) && !(!parz.IPTmultiple*v.at(ui).iptd ) 
		    && (gsl_rng_uniform(r) < parz.IPTcovPf()) &&
		    !v.at(ui).isD && !v.at(ui).isT && !v.at(ui).isIPT ){
		  event( ui, -4, parz, r );//IPT commencement
		}
		event( ui, -2, parz, r );//ART commencement - including shiv
		++noart;
	      }
	    }
	  }

	  // extra bit for passed linkers
	  if( v.at(ui).shiv && v.at(ui).hiv && !v.at(ui).isART
	      && (v.at(ui).artdefr == 0) && (time > parz.ARTst)
	      && v.at(ui).cd4 <= v.at(ui).cd4st ){ 
	    event( ui, -2, parz, r );//ART commencement
	    ++noart;
	  } // end linkers catchup

	} else {		// HIVMAC version of ART starts
	  // Notes: This means the hivmac flag overrules other information. I don't think the sigma parameter is going to be used.
	  if( gsl_rng_uniform(r) < 1-exp(-v.at(ui).artstrate * parz.dt ) && time > parz.glrst ){
	    event( ui, -2, parz, r );//ART commencement - including shiv
	    ++noart;
	    if( v.at(ui).age >= 15 ){++noart15;}
	    // now includes IPT
	    if( ( parz.IPTshiv || parz.IPTart ) && 
		(time > parz.IPTst) && !(!parz.IPTmultiple*v.at(ui).iptd ) 
		&& (gsl_rng_uniform(r) < parz.IPTcovPf()) &&
		!v.at(ui).isD && !v.at(ui).isT && !v.at(ui).isIPT ){
	      event( ui, -4, parz, r );//IPT commencement
	    }
	  } // end test for starting

	} // end whether HIVMAC version or not

      }	// end ART flag


      // TB disease outcomes
      if( v.at(ui).isD && v.at(ui).tsi > v.at(ui).tbout ){
	//	    lgf << "debug:"<< v.at(ui).dTBsc <<" "<< v.at(ui).dTBd<<" "<<v.at(ui).dTBdet << endl;
	//outcome due:
	if( v.at(ui).dTBsc ){
	  tbav += v.at(ui).tbout;
	  event( ui, 3, parz, r );//cure
	  ++ntbsc;
	}//self-cure
	if( v.at(ui).dTBd ){
	  tbav += v.at(ui).tbout;
	  if( v.at(ui).hiv ){++ntbdh;}
	  event( ui, 6, parz, r );//die
	  ++ntbd;
	}//death
	if( v.at(ui).dTBdet ){
	  if( v.at(ui).prevTB > 1 ){ ++ntbdetprev; } // detected TB that's had it before
	  if( v.at(ui).hiv ){++ntbdeth;}
	  if ( v.at(ui).sigma < 5 ){++ntbdetf;}
	  if( v.at(ui).tecf ){++nodtecf;}
	  if ( v.at(ui).isSmr ){++ntbdetsp;}
	  tbavd += v.at(ui).tbout;
	  test += event( ui, 4, parz, r );//detect, append hh prev
	  ++ntbdet;
	}//detection
      }//tb outcome
    }//not dead
  }//pop loop

  //  lgf << "DEBUG:" << test << endl;

  if(ntbdet>0){		// recording for detection
    propFD = (double)ntbdetf / ntbdet;
    propDH = (double)ntbdeth / ntbdet;
    propDS = (double)ntbdetsp / ntbdet;
    propECF = (double)nodtecf / ntbdet;//no via ECF
    prevhh /= (dethhno + 1e-10);	    // dividing by the number of hh folk
    //    prevHHH /= dethhno;
  } else {prevhh = -1;}// prevHHH = -1;}
  prevHHH = ( dethhno1549 > 0 ? prevHHH/dethhno1549 : -1 );


  avAdeath /= (bdcount+1e-10);
  tauav = tau/(nod+1e-10);

  if(!qlog){
    lgf << "\tav Rxtau= "<< taurx / (norx+1e-10) << "/" << norx << endl; 
    lgf << "\tav t-to other= "<<tbav/(ntbd+ntbsc+1e-10)<<" / "<< (ntbd+ntbsc) <<endl;
    lgf << "\tav t-to detec= "<<tbavd/(ntbdet+1e-10)<<" / "<< ntbdet <<endl;
    lgf << "\tno self-cures = " << ntbsc << endl;
    lgf << "\tno TB deaths = " << ntbd << " ("<< ntbdh <<":w/hiv) "<< endl;
    lgf << "\tno TB detections = " << ntbdet <<" ("<<ntbdeth<<":w/hiv)" ;
    lgf <<" ("<<ntbdetprev<<":recurrent)" << endl;
    lgf << "\tno routine ART starts = " << noart << endl;
    lgf << "\t hh TB prev = " << prevhh <<"/"<<dethhno << endl;
    lgf << "\t hh HH prev = " << prevHHH <<"/"<<dethhno1549 << endl;

    lgf << "\t bbs ART starts = " << bbsstarts << endl;
    lgf <<"\t(" << nod << ") HIV deaths done...hiv="<< hiv <<"("<<hiv1549;
    lgf <<"/" << size1549 << ")"<< endl;
    lgf <<"\t" << bdcount <<" deaths from senescence...rate/10^3=";
    lgf<< 1000*bdcount / (parz.dt*size) <<", mean age at death = "<<avAdeath<< endl;
  }
  if ( nod > 0 ){if(!parz.qlog){lgf << "\taverage HIV-TSI at death=" << tauav <<endl;}
    if( iptends > 0 && !parz.qlog){
      lgf << "\t IPT completions = " << iptends << endl;//IPT finished
    }
    //    lgf << "\t learnt HIV status : " << seauton << endl;
  }
      

  return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int population::event( unsigned int who, int etype, parameters& parz, const gsl_rng *r ){//event doer
  //update person
  int j = -1;			// set to who, except in household events
  // where set to whoami
  if ( etype >= 0){
    switch ( etype ){//SEE ALSO: DUPLICATE
    case 0:// birth types handled separately....
      {
	int loc = (int)who;
	unsigned int whoami = v.size();
	int sex(0);
	if ( gsl_rng_uniform(r) > 0.5 ){sex = 1;}
	switch( sex ){
	case 0:// birth - girl 
	  v.push_back( person( 0,  0.0 , noevs ) );//it's a girl!
	  hh.at( loc ).add( whoami, 0 );//add a girl
	  break;
	case 1:// birth - boy
	  v.push_back( person( 1,  0.0 , noevs ) );//it's a boy! 
	  hh.at( loc ).add( whoami, 1 );//add a boy
	  break;
	}//end switch for gender
	//both
	v.at( whoami ).hiv = 0; v.at( whoami ).tau = 0; 
	v.at( whoami ).taumax = 99; v.at( whoami ).is1549 = 0;
	if( parz.DHflag ){ //assign sexual activity 
	  double test = gsl_rng_uniform(r);
	  if( test < parz.sbFracH  ){
	    v.at( whoami ).scat = 1;
	  } else {
	    if( test - parz.sbFracH < parz.sbFracM ){
	      v.at( whoami ).scat = 2;
	    } else { v.at( whoami ).scat = 3; }
	  }
	} // end DHflag
	size++;//one more! - maintained elsewhere
	v.at( whoami ).hhid = loc;//where are they
	v.at( whoami ).LE = parz.lifeexpectancy( r, 0.0 );
	j = (int)whoami;
      }
      break;
    case 8:			// imr, similar to above
      {
	int loc = (int)who;
	who = gsl_rng_uniform_int( r, v.size() );//who to duplicate
	v.push_back( v.at(who) );
	unsigned int whoami = v.size();//one over
	--whoami;//resets who to the doppelganger
	//add to household
	hh.at( loc ).add( whoami, v.at(whoami).gender );//add a ?
	v.at( whoami ).hhid = loc;//where are they
	j = (int)whoami; 
      }
      break;
    case 9:			// rMtb infection
      if( !parz.Rflag ){break;}	// else do below
    case 1://MTB infection (stoch)
      {
	// protection done like this rather than HR because of
	// // differential protection if now S
	// v.at(who).dTBact = 1;
	// if( v.at(who).isPI ){
	//   if( !v.at(who).isL && !v.at(who).isrL ){
	//     if(gsl_rng_uniform(r) < parz.PprotnS){v.at(who).dTBact = 0;}
	//   } else {
	//     if(gsl_rng_uniform(r) < parz.Pprotn){v.at(who).dTBact = 0;}
	//   }
	// }		   // protection from infection
        int skip(0);
	if( v.at(who).isPI ){
	  if( !v.at(who).isL && !v.at(who).isrL ){
	    if(gsl_rng_uniform(r) < parz.PprotnS){skip = 1;} // protected
	  } else {
	    if(gsl_rng_uniform(r) < parz.Pprotn){skip = 1;} // protected
	  }
	}		   // protection from infection
        if( skip == 0 ){
          v.at(who).dTBact = 1;
          v.at(who).isS = 0;//NB needed before isL=1 to distinguish reinfection
          v.at(who).sigma = 0.0;
          v.at(who).isL = 1;
          if( etype == 1 ){ v.at(who).isrL = 0;} else {v.at(who).isrL = 1;}
          v.at(who).isPI = 1;
          v.at(who).isSmr = 0;
          v.at(who).isT = 0;
	
          //record ARI
          if( v.at(who).age < 12.0 && v.at(who).age > 5.0 ){ ++kidi; } else { ++adi; }
        }
	//note this is everytime
	j = (int)who;
      }
      break;
    case 2://activation to smr+/- (sched)   
      // DEBUG!
      // cout <<parz.time<<","<<v.at(who).hiv<<","<<v.at(who).fast<<",";
      // cout<<v.at(who).cd4<<"," <<v.at(who).sigma<<","<<v.at(who).sigmamax  <<endl;
      {
	v.at(who).isD = 1;
	if( v.at(who).isrL ){v.at(who).isrD = 1;}
	v.at(who).isS = 0;
	v.at(who).isL = 0; v.at(who).isrL = 0; 
	v.at(who).isPI = 1;
	v.at(who).isSmr = 0;//change!
	v.at(who).isT = 0;
	v.at(who).tsi = 0;	// set clock
	v.at(who).prevTB += 1;	// has had TB n times
	if( gsl_rng_uniform(r) < v.at(who).probsmrpos( parz ) ){
	  v.at(who).isSmr = 1;//smear positive
	}
	//destined for death or sc?, or detection
	v.at(who).setTBoutcome( r, parz, time );
	
	// counting...zerod in scheduled events
	double recentthreshold = 5.0;// CHANGE RECENT DEFN HERE!
	av2 += v.at(who).sigma;  	// for mean sigma at activation
	if ( v.at(who).sigma > recentthreshold ){++slow;} 
	if ( v.at(who).hiv ){
	  ++noah;
	  if ( v.at(who).sigma > recentthreshold ){ ++noahs; } // slow HIV
	  if ( v.at(who).isART ){ ++noaa; }
      	  av3 += v.at(who).sigma; // mean tb-tsi in HIV
      	  av4 += v.at(who).tau / v.at(who).taumax;//mean hiv-tsi propn in HIV
      	  av5 += v.at(who).cd4;
	}
	if( v.at(who).isSmr ){++noap;}//smr+
	if( v.at(who).fast ){
	  ++noaf;
	  if( v.at(who).hiv ){ ++noafp; } else { ++noafn; }
	}//via fast
	++noa;
	j = (int)who;
      }
      break;
    case 3://self-cure (sched)  
      {
	v.at(who).isL = 1;
	if( v.at(who).isrD ){v.at(who).isrL = 1; } else {v.at(who).isrL = 0;}
	v.at(who).isS = 0;
	v.at(who).isPI = 1;
	v.at(who).isD = 0; v.at(who).isrD = 0;
	v.at(who).isSmr = 0;
	v.at(who).isT = 0;
	v.at(who).sigma = 0;//NB self cure is considered to be like infection
	j = (int)who;
      }
      break;
    case 4://detection (detect) -- really Rx
      {
	v.at(who).isS = 0;
	v.at(who).isPI = 1;
	v.at(who).isD = 0; v.at(who).isrD = 0;
	v.at(who).isSmr = 0;
	v.at(who).isT = 1;
	v.at(who).isIPT = 0;	// IPT set to zero when starting Rx
	v.at(who).Rxtime = 0;

	// ART 2 TB policy
	if( time > parz.ART2TBst && v.at(who).hiv && !v.at(who).isART 
	    && gsl_rng_uniform(r) < parz.ART2TB 
	    //&& !v.at(who).ARTrefusenik 
	    ){ 
	  // starting at ART st!
	  event( who, -2, parz, r );		//put on ART
	}

	//HH statistics & interventions
	if( parz.hhFLAG ){//Includes coverage!
	  int loc = v.at(who).hhid;//where?
	  list< unsigned int >::iterator hmi;
	  for( hmi = hh.at(loc).members.begin(); hmi != hh.at(loc).members.end(); ++hmi ){
	    //go through household
	    unsigned int hm = *hmi;//housemate
	    if( hm != who ){//let's not be doing owt to our man
	      ++dethhno;// counter for HH folk
	      if( v.at(hm).is1549 ){ ++dethhno1549; }
	      if( v.at(hm).isD ){//another undetected case 
		prevhh += 1.0;
	      }
	      if( v.at(hm).hiv && v.at(hm).is1549 ){//HIV+ hh member
		prevHHH += 1.0;
	      }

	      if( parz.HHIflag ){	// 
		//interventions!
		if( time > parz.HHIst && time < parz.HHIet && gsl_rng_uniform(r) < parz.HHIcov ){
		  // BC record
		  v.at(who).phhi = 1; v.at(who).hhi = 1; v.at(who).tshh = 0;  //does this once each loop which is wasteful
		  v.at(hm).phhi = 1; v.at(hm).hhi = 1; v.at(hm).tshh = 0; // has had a household intervention

		  // IPT
		  if( parz.HIPflag && !v.at(hm).isD && 
		      !( !parz.IPTmultiple * v.at(hm).iptd ) &&
		      !v.at(hm).isT && !v.at(hm).isIPT){//hh IPT
		    if( v.at(hm).age < 6 ){// && gsl_rng_uniform(r) < parz.TSTu16sens ){
		      v.at(hm).hhipt = 1;
		      event( hm, -4, parz, r );//put on IPT
		    }
		    if( v.at(hm).age >= 18 && v.at(hm).hiv &&
			!( !parz.IPTmultiple * v.at(hm).iptd ) ){
		      v.at(hm).hhipt = 1;
		      event( hm, -4, parz, r );//put on IPT
		    }
		  }

		  // ART: this is obvserved -
		  // doesn't mean they wouldn't have been on ART anyway
		  if( parz.HARflag && v.at(hm).hiv & !v.at(hm).isART ){//hh ART
		    //if( !v.at(hm).dart ){ v.at(hm).hhart = 1;}//will go on ART from HH
		    v.at(hm).hhart = 1;
		    v.at(hm).cd4st += parz.cd4h; // shift in start criterion
		    // v.at(hm).dart = 0;//predestine for ART
		    double guide = (time < parz.dCD4st ? parz.cd4g : parz.cd4g2);
		    if( v.at(hm).cd4 < guide ){ // immediate start
		      event( hm, -2, parz, r );		//put on ART
		    }
		  }

		  // CD
		  if( v.at(hm).isD ){//another undetected case - don't recurse
		    v.at(hm).isS = 0; v.at(hm).isL = 0; v.at(hm).isPI = 1;v.at(hm).isD = 0; v.at(hm).isrD = 0;
		    v.at(hm).isSmr = 0;v.at(hm).isT = 1;v.at(hm).Rxtime = 0;
		    ++ntbdet;	// NB there are othercounters not updated
		  }//case put on Rx

		}//end intervention
	      }//end intervention flag

	    }//not our man
	  }//hhloop
	}//hhflag
	j = (int)who;
      }
      break;
    case 10:      //Rx completion success (sched)
    case 7:      //Rx completion failure (sched)
    case 5:			// clearance!
      {
	v.at(who).isS = 1;
	v.at(who).isPI = 1;
	v.at(who).isD = 0; v.at(who).isrD = 0;
	if( etype == 5 ){ v.at(who).isL = 0; v.at(who).isrL = 0;}
	if( etype == 7 ){	 // failure like new infection - CHANGE!
	  v.at(who).sigma = 0.0;//TSTI
	}
	v.at(who).isSmr = 0;
	v.at(who).isT = 0;
	v.at(who).isIPT = 0;	// a check
	j = (int)who;
      }
      break;
    case 6://.->X -- death
      {
	if ( !v.at( who ).isX ){//test not dead - otherwise rejection screws populations

	  // recording, for HIVMAC, using counters in parz
	  if( v.at(who).age >= 15 ){
	    if( v.at(who).hiv ){ // HIV+
	      if( v.at(who).isART ){ // ART+
		++parz.artp15x;
		if( v.at(who).isD || v.at(who).isT ){++parz.X15tbap;}
	      } else { 		// ART-
		if( v.at(who).GL == 1 ){++parz.hivpgp15x;} 
		else{++parz.hivpgn15x;}
		++parz.artn15x; 
		if( v.at(who).isD || v.at(who).isT ){++parz.X15tban;}
	      }
	    } else {		// HIV-
	      ++parz.hivn15x;
	      if( v.at(who).isD || v.at(who).isT ){++parz.X15tbn;}
	    }
	  } else {
	    ++parz.u15x;	
	  }
	  // end HIVMAC stuff

	  v.at(who).isS = 0;
	  v.at(who).isL = 0; v.at(who).isrL = 0;
	  v.at(who).isPI = 0;
	  v.at(who).isD = 0; v.at(who).isrD = 0;
	  v.at(who).isSmr = 0;
	  v.at(who).isT = 0;
	  v.at(who).isX = 1;
	  v.at(who).isIPT = 0;

	  try{
	    size--;//one less
	    if (v.at( who ).is1549 == 1 ){//no longer 15-49
	      //--size1549;//thing now unnecessary as in counts
	      v.at( who ).is1549 = 0;
	    }
	  }
	  catch(exception& e){
	    cout << e.what() <<" in 10: household main event"<<endl;
	  }
	  
	  if( parz.hhFLAG ){
	    int loc = v.at(who).hhid; 
	    try{
	      //no longer in household
	      hh.at( loc ).remove( who, v.at(who).gender );
	      v.at(who).hhid = -1;
	    }
	    catch(exception& e){
	      cout << e.what() <<" in 10: household gubbins"<<endl;
	      cout << loc<<","<<who<<","<<endl;
	    }
	  }//end hhflag
	  
	}//end not dead block
	j = (int)who;
      }
      break;

    default:
      cout << "**unrecognized positive event type in population's event!**" << endl;
      j = -1;
      exit(-1);
      break;
    }//end switch
    // most negative event types involve changes to TB hazard
    // this means calculating the pre-change, back-dated hazard for use in the reset
  } else {
    switch ( etype ){
      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -1:			// HIV infection 
      {

	// HIVMAC record
	if( v.at(who).age >=15 ){
	  ++parz.ihiv15;
	  if( v.at(who).GL == 1 ){++parz.iHIVgp;}
	}
	if( v.at(who).age >=15 && v.at(who).age <50  ){++parz.ihiv1549;}

	//infect
	v.at( who ).hiv = 1;
	v.at( who ).Lhiv = 1;
	v.at( who ).tau = 0;
	v.at( who ).taumax = gsl_ran_weibull( r, parz.thiv_beta, parz.thiv_alpha );//survival time
	v.at( who ).taumax0 = v.at( who ).taumax; // unchanging copy
	v.at( who ).cd4 *= 0.75;	// initial drop
	// ART starts
	v.at( who ).cd4st = gsl_ran_weibull( r, parz.artTl, parz.artTk ); //individual's thsh
	v.at( who ).cd4st = 500 - v.at( who ).cd4st;
	// upward shift in start threshold
	double gr = parz.dCD4ep / (parz.dCD4et - parz.dCD4st);
	double theta2 = v.at(who).cd4st;
	theta2 += (time + v.at(who).taumax - parz.dCD4st) * gr;
	theta2 /= (1 + v.at(who).taumax * gr / parz.cd41); // putative new threshold
	// now check time condition to see if using
	if( time > parz.dCD4st - (1-theta2/parz.cd41)*v.at(who).taumax){ // ie some change
	  if( time < parz.dCD4et - (1-theta2/parz.cd41)*v.at(who).taumax ){ // linear
	    v.at(who).cd4st = theta2; 
	  } else {		// after linear phase
	    v.at(who).cd4st += parz.dCD4ep;
	  }
	}
	// end upward shift in start threshold
	if( v.at( who ).cd4st <  0){ v.at( who ).cd4st = 0;} // safety
	// destined to learn HIV status earlier than going onto ART?
	double deltat = v.at(who).taumax * (1 - v.at(who).cd4st / parz.cd41); //t might stART?
	if( gsl_rng_uniform(r) < parz.shivp * parz.ARTcov(time + deltat) ){
	  v.at(who).dshiv = 1;	// destined to learn
	} // look fwd

	// mass-ART refusenik ism (but not for normal ART)
	if( gsl_rng_uniform(r) < parz.partrefusenik ){
	  v.at(who).ARTrefusenik = 1; 
	} else {v.at(who).ARTrefusenik = 0; }

	if( v.at( who ).isX ){cout << "have HIV-infected a dead person!" << endl;}
	j = (int)who;
      }
      break;
      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -2:			// starting ART 
      {

	// cout<<time <<","<<v.at(who).cd4<<","<<v.at(who).cd4st<< endl;
	// FOR GWEN! output to stdout: time, cd4@ART, cd4 if not linked first.
	v.at(who).shiv = 1;	    // can't be on ART without knowing status
	v.at(who).isART = 1;
	v.at(who).taustart = v.at(who).tau; // record time started
	v.at(who).cd4st = v.at(who).cd4; //  record _actual_ cd4 started
	//resets time of death to something sensible
	//NB can't be in Rinv, o/w happens on post ART infection

	v.at(who).taumax = v.at(who).tau + 
	  parz.extratime( (v.at(who).LE - v.at(who).age), v.at(who).cd4 );

	// need to include setting of cd42 here
	v.at(who).cd42 = parz.ARTcd4asymptote( v.at(who).cd4 );
	if( parz.art750 ){ v.at(who).cd42 = 750;}


	//is their disease going to be like HIV- or HIV+?
	if( gsl_rng_uniform(r) > parz.tbart ){
	  v.at(who).Lhiv = 0;
	}
	j = (int)who;
      }
      break;

      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -3:			// //true ECF infection 
      {
	//need to record detections via ECF
	v.at(who).tecf = 1;//touched by ECF

	//NB need something to deal with those already prevalent...
	// this makes more likely to detect, and faster
	if( v.at(who).isD ){//currently undetected prevalent TURNED OFF if -1!!
	  //       if( v.at(who).dTBdet ){//already destined for detection
	  // 	v.at(who).tbout *= parz.tecfdF;//reduce the time they take to ncome fwd
	  //       } else {//reflip to see if detect, and 
	  double p = v.at(who).probtbdetect( parz, time );
	  //p = parz.tecfdOR * p / ( 1 + (parz.tecfdOR-1) * p);
	  if( gsl_rng_uniform(r) < p ){
	    v.at(who).dTBdet = 1;//destined
	    v.at(who).dTBsc = 0; 
	    v.at(who).dTBd = 0;
	    v.at(who).tbout = v.at(who).detectime( r, parz, time );
	  }
	  //      }//end redestin
	}//end prevalent

	j = (int)who;
      }
      break;

      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -4:			// starting IPT 
      {
	// initial hazard

	
	// changes
	v.at( who ).tonipt = 0;	// initialize the ipt time
	v.at( who ).iptd = 1;	// record having been IPTd
	v.at( who ).isIPT = 1;	// record on IPT
	// retiming

	j = (int)who;
	//if( posser == 0 ){cout << "etype="<<etype<<endl;}//double posser =
      }
      break;

      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -5:			// ending IPT 
      {
	
	// changes
	v.at( who ).tonipt = -1;	// initialize the ipt time
	v.at( who ).isIPT = 0;	// record on IPT
	// CLEARANCE
	if( v.at( who ).isL && v.at( who ).hiv &&
	    gsl_rng_uniform(r) < parz.IPTcprobP ){
	  event( who, 5, parz, r ); // clear
	} 
	if( v.at( who ).isL && !v.at( who ).hiv &&
	    gsl_rng_uniform(r) < parz.IPTcprobN ){
	  event( who, 5, parz, r ); // clear
	} 


	// retiming

	j = (int)who;
	//if( posser == 0 ){cout << "etype="<<etype<<endl;}//double posser =
      }
      break;

      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    case -6:			// defaulting ART
      {
	v.at( who ).isART = 0;	 // no longer on ART
	++v.at( who ).artdefr;   // did default
	v.at( who ).Lhiv = 1;    // like HIV
	v.at( who ).shiv = 0;    // no longer linked to care
	v.at( who ).isIPT = 0;	 // also default IPT


	// don't need to set cd4 as this is done in person's update_age
	v.at( who ).taumax = v.at( who ).tau 
	  + (v.at( who ).taumax0 - v.at( who ).taustart);
	  // + v.at( who ).cd4st * v.at( who ).taumax0 / (0.75*parz.cd40) ;
	// extra time is previous cd4 / gradient
	j = (int)who;
      }
      break;

      // ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    default:
      cout <<"** unrecognised negative event-type! **" << endl;
      j = -1;
      exit(-1);
    }
  } // end negative events
  return j;//j records success or failure
}

int population::event( unsigned int who, string evtype, map < string, int >& evk, parameters& parz, const gsl_rng *r ){//event overloaded version
  // polymorphic wrapper to above
  return event( who, evk[evtype], parz, r );
}

///////////////////////////////////////////////////////////////////
// NEW imports from metapopulation

int population::initialinfect( const gsl_rng *r, ofstream& lgf, parameters& parz ){//initial infections

  //initially at equilibrium -- worry about smear positivity; and the prob of progression and age
  int S0(0), L0(0), D0(0), T0(0), H0(0);//initial states
  //for the logic of how these values are derived, see IBMnotes
  //todo -- these need to depend on the inference parameters + other things to vary by community...
  //todo -- NB difference between RSA and ZA for the birth rate and HIV prevalence etc. -> somit in metapopulation identifying country

  rL = 0;			// resistant, in class

  double b(0.05), d(0.01);
  //d = 50.0 / 1e+5;//arbitrary
  d = (double)parz.D0 / (parz.D0 + parz.S0 + parz.L0 + parz.T0);
  cout << "d= "<< 100*d <<"%" << endl;
  double F = parz.beta * d *.67;//factor 2 for smr +
  //    F = 0.01;			 // debug!
  double x = 1 / ( 1 + F/b );
  //x = 0.7; 			// debug! -- could try this on F also
  double y = 1.0 - x;
  cout << "F="<<F<<", y="<< y << endl;
  double avtime = ( parz.tbd_betan ) * parz.tbsurvn_L * .88; // weibull with 2...
  // ( parz.CDRn * parz.tbd_betan + (1-parz.CDRn) ) * parz.tbsurvn_L ;
  int N = (int)v.size();
  D0 = auxFloor( d * N );
  //    d = 0; D0=0;
  // T0 = auxFloor( d * parz.CDRn * N * parz.RxL / parz.tbsurvn_L );
  T0 = 1+auxFloor( d * parz.CDRn * N * parz.RxL / avtime);
  //T0 = D0;			// NEW SAFER
  cout << "initial D:T=" << (double)D0 / T0 << endl;
  N -= ( D0 + T0 );
  S0 = auxFloor( x * N ); 
  L0 = auxFloor( y * N ); 

  H0 = auxFloor( parz.hiv0 * N / 100 );

  //cout << "initial infect="<<D0<<endl; // 
  //new version age-related prevalence...
  //actually setting up
  lgf << "\tsetting infection in community, (beta,d;S0,L0,D0,T0)=("<<parz.beta<<","<< d <<";";
  lgf <<S0<<","<<L0<<","<<D0<<","<<T0<<")" << endl;
  lgf << "HIV0=" << H0 << endl;
  int who(0), dover(0);
  double sig(0),sigm(0);

  vector< double > dweights( v.size() , 0.0 );//based on age categories
  for( unsigned int ii = 0; ii < dweights.size(); ++ii ){//weights for L
    dweights.at(ii) = 1 - exp( - F * v.at(ii).age );
  }
  while( L0 > 0){//latent
    who = find_one_event( r, dweights );
    if( !v.at(who).isL ){//not already infected
      --L0;//reduce
      // change for both rest!
      if( parz.Rflag && gsl_rng_uniform(r) < parz.Rfrac0 ){
	event( who, 9, parz, r ); //rmtb infec
      } else {
	event( who, 1, parz, r ); //mtb infec
      }
      if( gsl_rng_uniform(r) < y * parz.Pprotn ){ // proxy for previous
	v.at(who).dTBact = 0;
      }
      v.at(who).fast = 0; // crucial not to have a fast dump
      v.at(who).sigma = v.at(who).age - 1/F;//nb this can be negative for kids

      //v.at(who).sigma = 0.5*v.at(who).age;
      //age as proxy for how long infected
      
      v.at(who).sigma = v.at(who).age + log(1-gsl_rng_uniform(r)*(1-exp(-F*v.at(who).age))) / F; // should be prob infected at this age given infected
      // not this is really when first infected

      //v.at(who).sigma = gsl_rng_uniform(r) * v.at(who).age;
      // new attempt

      //so no immediate burst
      sig += v.at(who).sigma;
    }
  }

  double N18(0);
  for( unsigned int ii = 0; ii < dweights.size(); ++ii ){//weights for prev
    if( v.at(ii).age < 18 ){ 
      dweights.at(ii) = 0.01;
    } else {
      dweights.at(ii) = 1.0;
      ++N18;		// another over 18
    }
  }

  N18 /= N;			// now the fraction of pop over 18
  // D0 = auxFloor( D0 * N18 );
  // T0 = auxFloor( T0 * N18 ); // these are the numbers rescaled

  while( D0 > 0 ){//the active diseased
    who = find_one_event( r, dweights );
    if( !v.at(who).isD ){//not already diseased
      --D0;//reduce
      event( who, 2, parz, r ); //mtb act
      v.at(who).tsi = gsl_rng_uniform(r) * v.at(who).tbout; 
      // avoid inital wave
    }
  }
  // NB an initial wave 
  
  while( T0 > 0){//those on Rx
    who = find_one_event( r, dweights );
    if( !v.at(who).isD && !v.at(who).isT ){//not already onRx
      --T0;
      v.at(who).isPI = 1;
      v.at(who).isT = 1;
      v.at(who).isS = 0;
      v.at(who).isL = 0;
      v.at(who).isSmr = 0;
      v.at(who).Rxtime = gsl_rng_uniform(r) * parz.RxL;
    }
  }

  // HIV
  if( parz.DHflag && H0 > 0 ){			// dynamic HIV
    for( unsigned int ii = 0; ii < dweights.size(); ++ii ){//weights for prev
      dweights.at(ii) = v.at(ii).HIVprob( parz );
    }

    while( H0 > 0){//HIV infections
      who = find_one_event( r, dweights );
      if( !v.at(who).hiv  ){//not already infected
	--H0;
	event( who, -1, parz, r );
      }
    }
  }

//     } else { //no age mixing


  //safefty
  for( unsigned int ii = 0; ii < v.size(); ++ii ){//
    if( v.at(ii).isS ){ ++popsN.at(0);}
    if( v.at(ii).isL ){ ++popsN.at(1);}
    if( v.at(ii).isD ){ ++popsN.at(2);}
    if( v.at(ii).isT ){ ++popsN.at(3);}
    if( v.at(ii).isrL ){ ++rL;}
  }

  popsP.at(0) = 0;
  popsP.at(1) = 0;
  popsP.at(2) = 0;
  popsP.at(3) = 0;
  hiv = 0;

  lgf << "GOT(SLDT):"<<popsN.at(0)<<"," <<popsN.at(1);
  lgf <<","<< popsN.at(2)<<","<<popsN.at(3)<<endl;

  L0 = popsN.at(1);
  
  lgf << "\t...commmunity done!" << endl;
  lgf << sig/(L0+1e-6) <<" av, max " << sigm/(L0+1e-6) << endl;
  lgf << dover << " overage" << endl;
  lgf << rL << " resistant latent" << endl;

  return hiv;
}
////////////////////////////////////////////////////////////////////////////////////////////

int population::reinitialize( const gsl_rng *r, ofstream& lgf, parameters& parz  ){//init only

  //need something to clear the population here!....
  // zero population!
  nohh = 0; size = 0;
  hh.clear(); 
  v.clear();
  parz.time = parz.starttime;	// absolultely essential

  // ---------
  int j = 0;
  lgf << "initializeing  the baseline population state..." << endl;
  int hhbar(0), gend(0), agi(0); //household bar code and gender
  double ag(0);//age
  int chh(0);//current community, current household
  unsigned int  totpopcount(0), tothhcount(0);//total counts
  unsigned int ui(0), id(0), who(0), loc(0), chhcount(0);

  size0 = parz.popdata.size();
  parz.popsize0 = size0;

  for(  ui = 0; ui < parz.popdata.size(); ++ui ){//loop down data
    //what is the data 
    // not using 0th element
    hhbar = parz.popdata.at(ui).at(1);
    gend = parz.popdata.at(ui).at(2);
    agi = parz.popdata.at(ui).at(3);
    //smooth age out
    ag = agi +  gsl_rng_uniform(r);

    if( chh != hhbar ){//a new household!
      chh = hhbar;//adjust current household
      ++chhcount;
      ++tothhcount;
      loc = chhcount - 1;
      if(!parz.hhFLAG){ loc = 0; } else {
	//add a new household
	hh.push_back( household() );//add an empty household
	hh.at( loc ).hiv = 0;// make sure you zzero foie and hiv
	hh.at( loc ).hfoi = 0;// 
	nohh++;//update no hhs
      }
    }

    try{
      v.push_back( person( gend, ag, noevs ) );//add a person to this population
      who = v.size();//just off end of current pop
      --who;
      ++totpopcount;

  // sexual activity classes
  for( unsigned int ii = 0; ii < v.size(); ++ii ){//
  }


      //HIV stuff
      v.at( who ).hiv = 0; v.at( who ).tau = 0; v.at( who ).taumax = 99; v.at( who ).cd4 = parz.cd40;
      if( (ag > 15) && (ag < 50) ){
	v.at( who ).is1549 = 1;
      } else { v.at( who ).is1549 = 1;}
      if( parz.DHflag ){ //assign sexual activity 
	double test = gsl_rng_uniform(r);
	if( test < parz.sbFracH  ){
	  v.at(who).scat = 1;
	} else {
	  if( test - parz.sbFracH < parz.sbFracM ){
	    v.at(who).scat = 2;
	  } else { v.at(who).scat = 3; }
	}
      }

      size++;//one more!

      // life-expectancy
      v.at( who ).LE = parz.lifeexpectancy( r, ag );
      //household -- doesn't matter if hhFLAG off...though it should never be...?
      hh.at( loc ).add( who, gend );//add a girl
      v.at( who ).hhid = loc;	//seems forgot to add hhid recording...

      //age categories -- don't use ag after this!
      if( parz.AGEflag ){//assign to correct age class
	//double tempage(0);
	if( ag > parz.agebins.back() ){ 
	  ag = parz.agebins.back()-0.1;
	} //else { tempage = ag; }//safety for if off end of agebins
	v.at(who).acat = bisect_find( ag, parz.agebins );//which age category?
      }//age classes done 
    }//end try
    catch(exception& e){
      lgf << e.what() <<" in metapop initialization..!"<<endl;
      lgf << "ui= "<<ui <<endl;
      lgf << "totpopcount,tothhcount,id,loc,who= ";
      lgf<<totpopcount<<","<<tothhcount<<","<<id<<","<<loc<<","<<who << endl;
      exit(-1);
    }
    
    j = (int) ui;
    
  }//end ui data loop

  //log progress
  lgf << " total population count = "<< totpopcount << endl;
  lgf << " total household count = "<< tothhcount << endl;

  if( parz.hhFLAG ){//calculating household data, may only need doing first time
    HHpostinit( r, lgf, parz );
  }//end hhflag

  lgf << "...population re-initialisation completed!" << endl;
  return j;

}


////////////////////////////////////////////////////////////////////////////////////////////

int population::logstate( ofstream& lgf, int runno, int i ){
  int j = 0;//success/failure
  lgf << "[" << i <<"] passed: t="<< time <<":(SLDTN)=(";
  lgf << (popsN.at(0)+popsP.at(0)) << ",";
  lgf <<  (popsN.at(1)+popsP.at(1)) << "," ;
  lgf <<  (popsN.at(2)+popsP.at(2)) << "," ;
  lgf <<  (popsN.at(3)+popsP.at(3)) << ",";
  lgf << size << ")"<< endl;
  return ++j;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int population::record( ofstream& lgf, results& rez, int runno, int i, parameters& parz ){
  int j = 0;//success/failure
  if( !rez.light ){		  	// full records
    rez.time.at(i) = time;//overwrites for several runs, but hey
    rez.ibmN[runno].at(i) = size; // total population
    rez.ibmS[runno].at(i) = popsN.at(0) + popsP.at(0);
    rez.ibmL[runno].at(i) = popsN.at(1) + popsP.at(1);

    rez.ibmrL[runno].at(i) = rL; // resistant strain

    rez.ibmLP[runno].at(i) = ( hiv>0 ? (double)popsP.at(1) / hiv : -1 ); // positive latent
    rez.ibmLN[runno].at(i) = (double)popsN.at(1) / (size-hiv + 1e-10); // -ve latent
    rez.ibmD[runno].at(i) = popsN.at(2) + popsP.at(2);
    rez.ibmD2[runno].at(i) = (double)prev18/ size18; // now in over 18 population
    rez.ibmT2[runno].at(i) = (double)prevT18/ size18;
    rez.ibmT[runno].at(i) = popsN.at(3) + popsP.at(3);
    
    //cout << "time="<< time.at(i)<<", T="<<ibmT2[runno].at(i)<<endl;

    rez.ibmHIV[runno].at(i) = hiv1549; // now made 15-49!
    rez.ibmHIV2[runno].at(i) = (double)hiv18 / size18;// HIV18+
    rez.ibmFOI[runno].at(i) = foibar; 
    rez.ibmFHH[runno].at(i) = ( foibar > 0 ? hhfoibar / foibar : 0 ); 
    rez.ibmPHH[runno].at(i) = prevhh;
    rez.ibmPHHH[runno].at(i) = prevHHH; // HIV prevalence in HH HHs(ie inc T)
    rez.ibmPhh[runno].at(i) = tbcoprev;
    rez.ibmPHhh[runno].at(i) = hivcoprev; // HIV prevalence in hiv+ HHs

    //if( kidi == 0 && i > 0){ibmARI[runno].at(i) = ibmARI[runno].at(i-1);} else {
    rez.ibmARI[runno].at(i) = (double)kidi / ( sizeyoung * parz.dt ); 
    rez.ibmARIA[runno].at(i) = (double)adi / ( ( size - sizeyoung ) * parz.dt ); 
    // convert from number infected to annual risk
    rez.ibmARI[runno].at(i) = 1.0 - exp( -rez.ibmARI[runno].at(i) );
    rez.ibmARIA[runno].at(i) = 1.0 - exp( -rez.ibmARIA[runno].at(i) );

    rez.ibm1549[runno].at(i) = size1549;
    rez.ibmu12[runno].at(i) = sizeyoung;
    rez.propThiv[runno].at(i) = (double)prevTH18 / (prevT18 + 1e-6); 
    rez.ibmHHPI[runno].at(i) = nophhi;
    rez.ibmHHart[runno].at(i) = nohhart;
    rez.ibmHHipt[runno].at(i) = nohhipt;

    //tECF
    rez.ibmpECF[runno].at(i) = notecf;   // no. touched by ECF
    rez.ibmvECF[runno].at(i) = propECF; // via ECF
    
    //incidence - rate!
    rez.Inc[runno].at(i) = noa / parz.dt;//number of activations
    rez.ibmSP[runno].at(i) = (double)prevSP18 / prev18; // now in prev

    //record proportions resetting if needbe: Fast Inc, Fast Dets,Smr+ Dets,Inc HIV+,Dets HIV+
    rez.ibmFI[runno].at(i) = propFI;
    rez.ibmFIP[runno].at(i) = propFIP;
    rez.ibmFD[runno].at(i) = propFD;
    rez.ibmDS[runno].at(i) = propDS; // smr+ in detected TB
    rez.ibmIH[runno].at(i) = propIH; // HIV in incident TB
    rez.ibmDH[runno].at(i) = propDH;

    // on ART
    rez.ibmARTp[runno].at(i) = (double)ARTnum / (hiv + 1e-5);	// proportion on ART
    rez.ibmARTe[runno].at(i) = (double)ARTnum / (ARTnum + ARTdenom + 1e-5);
    rez.ibmSHIVcount[runno].at(i) = shivcount / (hiv + 1e-5);

    // proportion on ART - of eligible

    //these need changing...
    rez.TBdeaths[runno].at(i) = (double)ntbd / (size * parz.dt ); //TB deaths-includes H+
    rez.HIVdeaths[runno].at(i) = (double)nod / (size * parz.dt ); //HIV deaths
    rez.TBHIVdeaths[runno].at(i) = (double)ntbdh / (size * parz.dt ); //TB-HIV deaths
    rez.tbnotes[runno].at(i) = (double)ntbdet/(parz.dt * size);//notifications per year  
    rez.ibmPRec[runno].at(i) = (double)ntbdetprev / (ntbdet + 1e-10); // proportion recur

    // IRRs
    if( noa > 0 ){
      // IRR for TB|HIV
      rez.IRR[runno].at(i) = ((double)noah / noa) * ( ((double)size - (double)hiv ) / hiv ) ;	
      // IRR for TB | ART vs no ART
      rez.IRRa[runno].at(i) = ((double)noaa / ((double)noah + 1e-6) ) * (((double)hiv) / ARTnum) ;	
    } else {
      rez.IRR[runno].at(i) = -1;
      rez.IRRa[runno].at(i) = -1;
    }
  }		
				// end heavy
  // stratified incidences
  rez.strIn[runno].at(i) = (double)(noa - noah) / ((size-hiv)*parz.dt);
  rez.strIp[runno].at(i) = (double)(noah) / ((hiv + 1e-6)*parz.dt);
  rez.strIa[runno].at(i) = (double)(noaa) / ((ARTnum + 1e-6)*parz.dt);
  rez.strIpna[runno].at(i) = ((double)(noah)-(double)(noaa)) / ((hiv-ARTnum + 1e-6)*parz.dt);
  // strat F or not
  rez.strIf[runno].at(i) = (double)(noaf) / (noa + 1e-3);
  rez.strIpf[runno].at(i) = (double)(noafp) / (noah + 1e-3);
  rez.strInf[runno].at(i) = (double)(noafn)/(noa-noah + 1e-3);

  // IPT counts
  rez.IPTprev[runno].at(i) = iptnos / size;
  rez.IPTprevP[runno].at(i) = iptnosP / (hiv + 1e-6);
  rez.IPTprevA[runno].at(i) = iptnosA / (ARTnum + 1e-6);
  rez.IPTprevN[runno].at(i) = (iptnos - iptnosP) / (size - hiv );
  rez.IPTch[runno].at(i) = (double)cipth / (hiv+1);

  // cd4 counts
  for( unsigned int k = 0; k < rez.cd4hist.at(0).size(); ++k ){
    rez.cd4hist.at(i).at(k) = cd4hist.at(k);
  }


  // 6 month incidence
  int halfyr = auxCeil(0.5/parz.dt);
  if( i >= (parz.notime-halfyr) ){
    if( i == (parz.notime-halfyr) ){
      rez.Itmp = rez.Im;
      rez.aritmp = rez.ARIm;
    }
    rez.Im += noa / (parz.dt * halfyr * parz.noruns); // mean inc over last 6 mos
    rez.EP3.at(runno) += noa / (size * parz.dt * halfyr ); // mean inc over last 6 mos, each run
    rez.ARIm += (double)kidi / ( ( size - sizeyoung ) * parz.dt * parz.noruns ); 
    rez.EP2.at(runno) += (double)kidi / ( sizeyoung * parz.dt  * halfyr ); 
    rez.mortz.at(runno).at(0) += (double)bdcount / (size * parz.dt  * halfyr); // bckgrnd deaths
    rez.mortz.at(runno).at(1) += (double)ntbd / (size * parz.dt  * halfyr); // TB deaths
    rez.mortz.at(runno).at(2) += (double)nod / (size * parz.dt  * halfyr); // HIV deaths
    rez.mortz.at(runno).at(3) += (double)ntbdh / (size * parz.dt  * halfyr); // TB-HIV deaths
  }
  // NB the 'ARI' is 6mo mean FOI in kids
  // EP1 is prevalence
  // EP2 is kid ari 6 mo
  // EP3 is inc disease 6 mo
  // mortz are deaths: bckgrnd, TB, HIV, TB&HIV

  // recording some new info on end-point snapshots
  if( i == (parz.notime-1) ){	// only do at the end
    // easiest approach to take a snapshot of the population
    // at run 1 and output it for analysis in R
    snapshot( lgf, rez.ibmstate );   //record the population at the moment
    rez.Isqm += pow((rez.Im-rez.Itmp),2) * parz.noruns; // nrun*(Im-Itmp)=Ibar for one run, then divide
    double dl = (double)prev18/ size18; 
    rez.Dsqm += pow( dl, 2 ) / parz.noruns; // the prevalence squared
    rez.Dm += dl / parz.noruns; // the prevalence
    rez.EP1.at(runno) = dl;
    rez.ARIsqm += pow( (rez.ARIm - rez.aritmp), 2 ) * parz.noruns; // the ARI squared


  }

  j++;
  return j;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int population::demosnapshot( ofstream& lgf, results& rez, int runno, int i, parameters& parz ){
  //DEMOGRAPHIC SNAPSHOT by household
  int j(0);
  if ( runno == 0) {//first run only
    rez.hh0tot = 0; 
    vector< unsigned int > M( hh.size(), 0 ), W( hh.size(), 0 ), T( hh.size(), 0 );
    for ( unsigned int ii = 0; ii < hh.size(); ++ii ){
      M.at(ii) = hh.at(ii).nmen;
      W.at(ii) = hh.at(ii).nwomen;
      T.at(ii) = hh.at(ii).nwomen + hh.at(ii).nmen;
      if ( T.at(ii) == 0 ){ ++rez.hh0tot; }
    }
    rez.hhtot = (int)hh.size();

    // maxes and mins
    rez.Tmax = *max_element( T.begin(), T.end() );
    rez.Tmin = *min_element( T.begin(), T.end() );
    rez.Mmax = *max_element( M.begin(), M.end() );
    rez.Mmin = *min_element( M.begin(), M.end() );
    rez.Wmax = *max_element( W.begin(), W.end() );
    rez.Wmin = *min_element( W.begin(), W.end() );

    // zero P
    for(unsigned int u = 0; u < rez.P.size(); ++u){
      for(unsigned int v = 0; v < rez.P.at(u).size(); ++v){
	rez.P.at(u).at(v) = 0.0;
      }
    }

    //matrix of hh sizes n - women, m - men
    for ( unsigned int ii = 0; ii < hh.size(); ++ii ){
      if( hh.at(ii).nwomen < 10 && hh.at(ii).nmen < 10 ){
	++rez.P.at( hh.at(ii).nwomen ).at( hh.at(ii).nmen );
      }
    }
    rez.P.at(0).at(0) = 0;

    ////record age-mixing in households? todo?
    // ACTUALLY WRITE RESULTS
    j = rez.writedss( lgf, i, parz );
  }
  return j;
}
///////////////////////////////////////////////////////////////////////////////////
int population::event_sweep(   map < string, int >& swk, map < string, int >& evk, map < int, int >& ghk, const gsl_rng *r, ofstream& lgf, parameters& parz ){
  // this method is for setting up weights and then enacting events of a given type

  vector< double > weights( v.size() );// for choosing from among eligible
  vector< unsigned int > where( v.size(), 0); // worry about redeclarations?
  // these will need removing from ibmrun when complete : todo?
  unsigned int  noe(0), num(0);
  unsigned int who(0);
  int j(-1);
  list< unsigned int > event_list;
  int rans(0);			// total number of events
  map < string, int >::iterator kp;

  for( kp = swk.begin(); kp != swk.end(); ++kp ){ // loop event types

    // resets
    weights.assign( v.size(),0.0 );
    where.assign( v.size(), 0); 
    j = -1; noe = 0; num = 0;

    // setup using event type for setup_weights
    setup_weights(lgf, parz, kp->second, ghk, weights, where, noe, num, r);// where

    if(!parz.qlog){
      lgf << "\tstochastic-sw event "<< kp->first <<": " << num << " required from "<< noe << endl;
    }

    //int etype = evk[kp->first];	// event type for ::event
    int etype(0);
    // actual event sweep
    try{
      if( num > noe ){ num = noe; }//safety
      if ( num > 0 && noe > 0 ){//not no events loop
	event_list.clear();
	find_events( r, event_list, noe, num, weights );//decide who
	while ( !event_list.empty() ){//go through 
	  who = where.at( event_list.back() );//who's to be done
	  event_list.pop_back();//take off list
	  //j = event( who, etype, parz, r );//where and what
	  j = event( who, kp->first, evk, parz, r );//where and what
	  // NB j is who in all cases, including birth/imr
	  //double foiarg[] = { 0.0, FOI, 0.0, rFOI};
	  double foiarg[] = { 0.0, 0, 0.0, 0};
	  v.at( j ).update_haz( foiarg, parz, time );
	  //update individuals hazards, no worry hh - double hit rare
	  if (j < 0) {if(!parz.qlog){lgf << "\t\t ** problem in enacting event "<< etype<<"**"<< endl;}}
	}//end event-while
      }//end not no events
    }  // end try
    catch( exception& ex){
      cout << ex.what() << " in event "<<etype<<"!" << endl;
      cout << noe<<","<<num<<endl;
      cout << who<<","<<v.size() <<endl;
      cout << event_list.back()<<","<<where.size()<<endl;
      exit(-1);
    }
    rans += num;
  } // end loop over events

  return rans;
}

////////////////////////////////GUTS OF IBM///////////////////////////////////////////////////
//IBM RUNNER
int population::ibmrun( const gsl_rng *r, ofstream& lgf, parameters& parz, results& rez, int runno){

  bool qlog = parz.qlog;
  if(!qlog){lgf << "starting IBM..." << endl;}

  // EVENT CODINGS
  map < string, int > evk, swk;	// map to event and sw
  map < int, int > ghk;		// map between sw and gh
  // GH: (0,1=foi,2=a,3)
  // keys for events, setup_weights, and the ghaz global hazards
  evk.insert(make_pair<string,int>("birth",0));   // -- birth
  swk.insert(make_pair<string,int>("birth",-2)); 
  evk.insert(make_pair<string,int>("imr",8));	  // -- inmigration
  swk.insert(make_pair<string,int>("imr",-5));	  
  evk.insert(make_pair<string,int>("TBi",1));	  // -- TB infection
  swk.insert(make_pair<string,int>("TBi",1));	 
  ghk.insert(make_pair<int,int>(1,1));	 
  evk.insert(make_pair<string,int>("TBa",2));	  // -- TB activation
  swk.insert(make_pair<string,int>("TBa",2));	  
  ghk.insert(make_pair<int,int>(2,2));	  
  evk.insert(make_pair<string,int>("sc",3));	  // -- TB selfcure
  evk.insert(make_pair<string,int>("Rxs",4));	  // -- Rx start
  evk.insert(make_pair<string,int>("RxC",5));	  // -- Rx clearance
  evk.insert(make_pair<string,int>("death",6));   // -- death
  evk.insert(make_pair<string,int>("corx",6));    // -- death corx
  swk.insert(make_pair<string,int>("corx",11)); 
  evk.insert(make_pair<string,int>("omr",6));     // -- outmigration
  swk.insert(make_pair<string,int>("omr",3)); 
  evk.insert(make_pair<string,int>("RxF",7));	  // -- Rx failure
  evk.insert(make_pair<string,int>("HIVi",-1));   // -- HIV infection
  swk.insert(make_pair<string,int>("HIVi",-1));
  ghk.insert(make_pair<int,int>(-1,0));
  evk.insert(make_pair<string,int>("ARTs",-2));   // -- ART start
  evk.insert(make_pair<string,int>("ECFi",-3));   // -- ECF infection
  swk.insert(make_pair<string,int>("ECFi",-3)); 
  evk.insert(make_pair<string,int>("IPTs",-4));   // -- IPT start
  evk.insert(make_pair<string,int>("IPTe",-5));   // -- IPT end
  evk.insert(make_pair<string,int>("ARTdef",-6)); // -- ART default
  swk.insert(make_pair<string,int>("ARTdef",-4)); 
  evk.insert(make_pair<string,int>("Rxe",10));	  // -- Rx success
  // a resistant strain
  evk.insert(make_pair<string,int>("rTBi",9));	  // -- rTB infection
  swk.insert(make_pair<string,int>("rTBi",4));	 
  ghk.insert(make_pair<int,int>(4,3));	 

  int  i;
  double nextartround(3000), nextPACFround(3000), nextHIVround(3000), nextIPTround(3000);	// if there's going to be TnT/PACF/HIV screening/mass IPT
  if( parz.massart ){ nextartround = parz.artstart; } // time of first round
  if( parz.PACFlag ){ nextPACFround = parz.PACFst; } // time of first round
  if( parz.HIVflag ){ nextHIVround = parz.HIVst; } // time of first round
  if( parz.IPTflag ){ nextIPTround = parz.IPTst; } // time of first round

  //loop over runs...
  if(!qlog){
    lgf << "starting IBM run "<< runno <<"..." << endl;
    cout << "starting IBM run "<< runno <<"..." << endl;
  }

  time = parz.starttime;
  counts( lgf, parz );//popsP,popsN, maybe others

  //-----------------------starting time loop..
  if(!qlog){    lgf << "notimes = " << parz.notime << endl;}
  for(i = 0; i < parz.notime; ++i){//timeloop

    //    if( i == 150 ){cout << "RNG state @(i=150) :"<< gsl_rng_uniform(r) << endl; }//if( i == 200 ){cout << "RNG state @(i=200) :"<< gsl_rng_uniform(r) << endl; }//cout << i <<endl;//debug

    //TIMES & AGES
    update_ages( r, parz );//update ages 

    //SCHEDULED EVENTS
    if(!qlog){lgf << "\tscheduled events..." << endl;}
    // NB before stoch events so zeros
    scheduledevents( lgf, r, qlog, parz, evk );//this also zeros counters

    //STOCHASTIC EVENTS
    if(!qlog){lgf << "\t stochastic events..." << endl;}
    // updating hazards for events using global hazards
    update_hazards( lgf, parz );	   //update hazard!

    // events requiring setup_weights
    event_sweep( swk, evk, ghk, r, lgf, parz );

    //HOUSEHOLD MOVEMENT
    if ( parz.hhFLAG ){		// building new hhs
      build_new_households( parz, nob, r );
      if(!qlog){lgf << "\t(" << newhhcount <<") households added..." << endl; }
    }

    try{
      if( parz.hhFLAG ){
	hhcounts( lgf, parz );//do the hh counts 
	if( !qlog ){ lgf << "\t  hh counts updated..." << endl; }
	//household movements
	int num = gsl_ran_poisson( r, parz.hhmoverate * parz.dt * size );
	hhmovement( r, lgf, num, parz, false );//nb no action on weights...
      }
    }
    catch(exception& A){
      cout << A.what() << "in household movement!" << endl;
      exit(-1);
    }

    //INTERVENTIONS    -- todo: stick in own method?
    try{

      if( parz.IPTflag && time > nextIPTround && time < parz.IPTet ){
	nextIPTround += parz.IPTT;
	if(!qlog){lgf << "mass IPT round! (at t="<< time <<")" <<endl;}
	int iptcount(0);
	for( unsigned int ug = 0; ug < v.size(); ++ug ){
	  if( !v.at(ug).isIPT && !v.at(ug).isX && 
	      !v.at(ug).isD && !v.at(ug).isT &&
	      !( !parz.IPTmultiple * v.at(ug).iptd ) 
	      ){
	    if( gsl_rng_uniform(r) < parz.IPTcov && 
		gsl_rng_uniform(r) < (parz.TSTsens ? v.at(ug).isL : 1-parz.TSTspec)){ // positive TST
	      event( ug, -4, parz, r );//PI remains 1
	      ++iptcount;
	    }
	  }//eligible
	}//popsize
	if(!qlog){lgf << "successfully, instantly IPTd: " << iptcount << endl;}
      }


      //initial sweep of ART
       if( parz.massart && time > nextartround ){
	 nextartround += parz.artT; // time for next round
	 if(!qlog){lgf << "mass round of test/ART! (at t="<< time <<")" <<endl;}
      	int artcount(0);
      	for( unsigned int ug = 0; ug < v.size(); ++ug ){
      	  if( v.at(ug).hiv && !v.at(ug).isART && 
	      !v.at(ug).ARTrefusenik && !v.at(ug).isX ){//hiv+&alive
       	    if( gsl_rng_uniform(r) < parz.artP ){
	      event( ug, -2, parz, r );		//put on ART
      	      ++artcount;
       	    } // else { 		// chance of becoming permanent refusenik
	    //   if( gsl_rng_uniform(r) < parz.partrefusenik ){
	    // 	v.at(ug).ARTrefusenik = 1; 
	    //   }
	    // } 
      	  } // NB changed refusenik to at HIV infection
      	}
      	if(!qlog){lgf << "test/ART sweep found: " << artcount << endl; }
      }

       // periodic ACF
       if( parz.PACFlag && time > nextPACFround && time < parz.PACFet ){
	 nextPACFround += parz.PACFT; // time for next round
	 if(!qlog){lgf << "PACF round! (at t="<< time <<")" <<endl;}
      	int acfcount(0);
      	for( unsigned int ug = 0; ug < v.size(); ++ug ){
	  if( v.at(ug).isD ){ // undetected TB
	    double test(0);
	    if( !v.at(ug).isSmr ){ 
	      test = parz.PACFeffsn; 
	    } else {
	      test = parz.PACFeffsp;
	    }
	    if( v.at(ug).hiv && !v.at(ug).isART ){ // include ART?
	      double fac = parz.PACFhivOR;
	      test = fac * test / (1 + (fac-1) * test);
	    }
	    if( v.at(ug).isART ){ test = parz.PACFeffart;}

	    if( gsl_rng_uniform(r) < test ){
	      event( ug, 4, parz, r );		//detect!
	      ++acfcount;
	    }
	  }
      	}
      	if(!qlog){lgf << "ACF sweep found: " << acfcount << endl; }
       }


       // periodic HIV screening
       if( parz.HIVflag && time > nextHIVround && time < parz.HIVet ){
	 nextHIVround += parz.HIVT; // time for next round
	 if(!qlog){lgf << "HIV screening round! (at t="<< time <<")" <<endl;}
      	int hivcount(0);
      	for( unsigned int ug = 0; ug < v.size(); ++ug ){
	  if( v.at(ug).hiv && !v.at(ug).shiv ){ // undetected TB
	    if( gsl_rng_uniform(r) < parz.HIVeff ){
	      v.at(ug).shiv = 1; // now linked to care - will increase ART coverage, and linked coverage
	      ++hivcount;
	    }
	  }
      	}
      	if(!qlog){lgf << "HIV sweep found: " << hivcount << endl; }
       }


    }
    catch(exception& A){
      cout << A.what() << " in interventions!" << endl;
      exit(-1);
    }


    //RECORDING RESULTS
    counts( lgf, parz );//popsP,popsN, maybe others
    record( lgf, rez, runno, i, parz );  
    if( !qlog ){ logstate( lgf, runno, i ); }//record state in logfile 

    }//end timeloop i
    if(!qlog){lgf << "\t\tIBM run: "<< runno << ",done" << endl;}
    if(!qlog){lgf << "\tfraction dead = " << ( (double) size ) / v.size() << endl;}
    // if(!qlog){lgf <<"\t"<< hh.size() << " households in total..." << endl;}
    cout << size << endl;
  return runno;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//end population member fun
