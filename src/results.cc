//#######################
#include "results.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//RESULTS MEMBER FUNS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int results::analyze( ofstream& lgf, parameters& parz ){//averaging over runs, complicated by NA
  
  int n = parz.noruns;
  if(  !light && n > 1 ){	// heavy recording too

    // double ones...
    // non-zero HIV prevalence
    collapse_result_averager( ibmLP, ibmHIV );
    collapse_result_averager( ibmARTp, ibmHIV );
    collapse_result_averager( ibmARTe, ibmHIV );
    collapse_result_averager( ibmSHIVcount, ibmHIV );
    // non-zero notifs
    collapse_result_averager( ibmPHH, tbnotes ); 
    collapse_result_averager( ibmPHHH, tbnotes );  
    collapse_result_averager( ibmFD, tbnotes );
    collapse_result_averager( ibmDS, tbnotes );
    collapse_result_averager( ibmDH, tbnotes );
    collapse_result_averager( ibmvECF, tbnotes );
    collapse_result_averager( ibmPRec, tbnotes );
    // non-zero hiv-TB incidence - NB need doing before ibmIH
    collapse_result_averager( ibmFIP, ibmIH );
    // non-zero incidence 
    collapse_result_averager( ibmSP, Inc );
    collapse_result_averager( ibmFI, Inc );
    collapse_result_averager( ibmIH, Inc ); // ibmIH
    collapse_result_averager( IRR, Inc );
    collapse_result_averager( IRRa, Inc );

    // new averaging - single dependence
    collapse_result_averager( ibmN );
    collapse_result_averager( ibmS );
    collapse_result_averager( ibmL );
    collapse_result_averager( ibmrL ); // resistant latent prevalence
    collapse_result_averager( ibmLN );
    collapse_result_averager( ibmD );
    collapse_result_averager( ibmD2 );
    collapse_result_averager( ibmT );
    collapse_result_averager( ibmT2 );
    collapse_result_averager( ibmHIV );
    collapse_result_averager( ibmHIV2 );
    collapse_result_averager( ibmFOI );
    collapse_result_averager( ibmFHH );
    collapse_result_averager( ibmARI );
    collapse_result_averager( ibmARIA );
    collapse_result_averager( ibmARIA );
    collapse_result_averager( ibm1549 );
    collapse_result_averager( ibmu12 );
    collapse_result_averager( propThiv );
    collapse_result_averager( TBdeaths );
    collapse_result_averager( HIVdeaths );
    collapse_result_averager( TBHIVdeaths );
    collapse_result_averager( tbnotes );
    collapse_result_averager( ibmHHPI );
    collapse_result_averager( ibmHHart );
    collapse_result_averager( ibmHHipt );
    collapse_result_averager( ibmpECF );
    collapse_result_averager( Inc );
    // stratified
    collapse_result_averager( strIn );
    collapse_result_averager( strIp );
    collapse_result_averager( strIa );
    collapse_result_averager( strIpna );
    collapse_result_averager( strIpf );
    collapse_result_averager( strInf );
    collapse_result_averager( strIf );
    // IPT prevalences
    collapse_result_averager( IPTprev );
    collapse_result_averager( IPTprevP );
    collapse_result_averager( IPTprevA );
    collapse_result_averager( IPTprevN );
    collapse_result_averager( IPTch );

    // smooth strIs and use for IRR!
    vector< double > tempp, tempn, tempa, temppna;
    MAsmooth( strIn[0], tempn, 10 );
    MAsmooth( strIp[0], tempp, 10 );
    MAsmooth( strIa[0], tempa, 10 );
    MAsmooth( strIpna[0], temppna, 10 );

    // cout << "test="<<tempp.at(tempp.size()-1) << endl;
    // it seems these are calculated here, 
    // despite gathering data along the way
    for ( int i = 0; i < (int)IRR[0].size(); ++i ){
      IRR[0].at(i) = tempp.at(i) / tempn.at(i);
      IRRa[0].at(i) = tempa.at(i) / (temppna.at(i) + 1e-6);
    }

  }//end n > 0
  return n;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int results::write( ofstream& lgf, parameters& parz, time_t& writetime ){

  lgf <<  "beginning to write results at writetime="<< writetime << endl;

  int j = 0;//recording success
  ofstream res_strm;//the result file text
  bool uniqueoutput = true;	// timestamping
  ostringstream oss;
  lgf << "opening endpoint file(6) ..." << endl;
  if( uniqueoutput ){
    oss << parz.wkdir << parz.jobid  << "EP.dat";
    res_strm.open(  oss.str().c_str() );
  } else {
    res_strm.open(  (parz.wkdir + "EP.dat").c_str() );
  }
  if(!res_strm){lgf << "EP(6) file failed to open..." << endl;} else {
    lgf << "EP(6) file opened..." << endl;
    res_strm << "noruns,avD,avDsq,avARI,avARIsq,avIbar,avIbar2" << endl;
    res_strm << parz.noruns << ",";
    res_strm << Dm <<"," <<Dsqm << ",";
    res_strm << ARIm <<"," << ARIsqm << ",";
    res_strm << Im <<"," << Isqm << endl;
    // if( Ibarsq.size() > 1 ){
    //   res_strm << Ibarsq.at(0) <<"," << Ibarsq.at(1) << endl;
    // } else {
    //   res_strm << Ibarsq.at(0) <<"," << "NA" << endl;
    // }
    res_strm.close();
    lgf << "closing EP(6)..." << endl;
  }

  lgf << "opening endpoint file EP1 ..." << endl;
  if( uniqueoutput ){
    oss.str("");		// clear
    oss << parz.wkdir << parz.jobid  << "EP1.dat";
    cout << oss.str().c_str() << endl;
    res_strm.open(  oss.str().c_str() );
  } else {
    res_strm.open(  (parz.wkdir + "EP1.dat").c_str() );
  }
  if(!res_strm){lgf << "EP1 file failed to open..." << endl;} else {
    lgf << "EP1 file opened..." << endl;
    for(int i = 0; i < parz.noruns; ++i){ // record
      res_strm << EP1.at(i) << endl;
    }
    res_strm.close();
    lgf << "closing EP1..." << endl;
  }

  lgf << "opening endpoint file EP2 ..." << endl;
  if( uniqueoutput ){
    oss.str("");		// clear
    oss << parz.wkdir <<  parz.jobid  << "EP2.dat";
    res_strm.open(  oss.str().c_str() );
  } else {
    res_strm.open(  (parz.wkdir + "EP2.dat").c_str() );
  }
  if(!res_strm){lgf << "EP2 file failed to open..." << endl;} else {
    lgf << "EP2 file opened..." << endl;
    for(int i = 0; i < parz.noruns; ++i){ // record
      res_strm << EP2.at(i) << endl;
    }
    res_strm.close();
    lgf << "closing EP2..." << endl;
  }


  lgf << "opening endpoint file EP3 ..." << endl;
  if( uniqueoutput ){
    oss.str("");		// clear
    oss << parz.wkdir << parz.jobid << "EP3.dat";
    res_strm.open(  oss.str().c_str() );
  } else {
    res_strm.open(  (parz.wkdir + "EP3.dat").c_str() );
  } 
  if(!res_strm){lgf << "EP3 file failed to open..." << endl;} else {
    lgf << "EP3 file opened..." << endl;
    for(int i = 0; i < parz.noruns; ++i){ // record
      res_strm << EP3.at(i) << endl;
    }
    res_strm.close();
    lgf << "closing EP3..." << endl;
  }


  lgf << "opening endpoint file EPm ..." << endl;
  if( uniqueoutput ){
    oss.str("");		// clear
    oss << parz.wkdir << parz.jobid << "EPm.dat";
    res_strm.open(  oss.str().c_str() );
  } else {
    res_strm.open(  (parz.wkdir + "EPm.dat").c_str() );
  } 
  if(!res_strm){lgf << "EPm file failed to open..." << endl;} else {
    lgf << "EPm file opened..." << endl;
    for(int i = 0; i < parz.noruns; ++i){ // record
      res_strm << mortz.at(i).at(0) <<"\t"<< mortz.at(i).at(1) <<"\t";
      res_strm << mortz.at(i).at(2) <<"\t"<< mortz.at(i).at(3) << endl;
    }
    res_strm.close();
    lgf << "closing EPm..." << endl;
  }


  if( !light ){

    // NB anything with a -1 in should not be smoothed!

    //smooth notifs!
    MAsmooth( tbnotes[0] , 3 );	// notifications
    MAsmooth( Inc[0], 3 );	// TB incidence
    MAsmooth( ibmARI[0] , 3 );	// ARI
    MAsmooth( ibmARIA[0] , 3 );	// ARIA
    //NB - these can make shocks look like they propagate backwards!!!
    // MAsmooth( TBdeaths[0], 10 );	// TB deaths
    // MAsmooth( HIVdeaths[0], 10 );	// HIV deaths
    // MAsmooth( TBHIVdeaths[0], 10 );	// TB/HIV deaths
    MAsmooth( strIn[0], 10 ); //see also analyze!
    MAsmooth( strIp[0], 10 );
    MAsmooth( strIa[0], 10 );
    MAsmooth( strIpf[0], 10 );MAsmooth( strInf[0], 10 );MAsmooth( strIf[0], 10 );
    // already commented below

    //  MAsmooth( propThiv[0] , 3 );
    //   MAsmooth( ibmFI[0] , 6 ); //MAsmooth not appropriate
    //   MAsmooth( ibmDS[0] , 10 );
    //   MAsmooth( ibmSP[0] , 6 );
    //  MAsmooth( ibmSP[0], 10 );	// smr+
    //  MAsmooth( ibmPHHH[0], 5 );	// HIV in HH 1549 -1!
    //  MAsmooth( ibmPHH[0], 10 );	// TB prev in HH -1!
    //  MAsmooth( ibmFI[0], 15 );	// recent in incident TB -1!

    double N;
    
    //TIMESERIES
    lgf << "opening results file 1 ..." << endl;
    res_strm.open( (parz.wkdir + "results1.dat").c_str() );
    if(!res_strm){lgf << "results1 file failed to open..." << endl;} else {
      lgf << "results file opened..." << endl;
      int i;
      //double N0 = ( ibmS[0].at(0) + ibmL[0].at(0) + ibmD[0].at(0) + ibmT[0].at(0) );
      double N0 = parz.popsize0;
      for(i = 0; i < notime; ++i){//t,deltaN,L,hiv,midage,hivinT
	//N = ( ibmS[0].at(i) + ibmL[0].at(i) + ibmD[0].at(i) + ibmT[0].at(i) );
	N = ibmN[0].at(i);
	res_strm <<  time.at(i) <<"\t" << (N-N0)/N0 << "\t" << ibmL[0].at(i)/N << "\t";
	res_strm <<  ibmHIV[0].at(i)/ ibm1549[0].at(i) <<"\t"<<ibm1549[0].at(i)/N;
	res_strm << "\t" << propThiv[0].at(i) << "\t" << ibmu12[0].at(i)/N;
	res_strm<<"\t"<< ibmHIV2[0].at(i)<<"\t";
	if( ibmLP[0].at(i) > 0 ){ res_strm<<ibmLP[0].at(i)<<"\t"; } else {res_strm<<"NaN\t";}
	res_strm<<ibmLN[0].at(i)<<"\t"<<ibmrL[0].at(i)/N;
	res_strm <<endl;
      }
      res_strm.close();
      lgf << "closing results 1..." << endl;
      if ( i == notime ){ j += 1;}
    }

    // debug
    // cout << "(run,pop)=("<< endl;
    // for( unsigned int ii = 0; ii < ibmN.size(); ++ii ){
    //   cout << ibmN.at(ii).at(notime-1)<<",";
    // }
    // cout<<")"<<endl;

    lgf << "opening results file 2 ..." << endl;
    res_strm.open( (parz.wkdir + "results2.dat").c_str() );
    if(!res_strm){lgf << "results2 file failed to open..." << endl;} else {
      lgf << "results2 file opened..." << endl;
      for(int i = 0; i < notime; ++i){//t,fast,smr+prev,smr+detections
  	res_strm <<  time.at(i) << "\t";
	if( ibmFI[0].at(i) > 0 ){ res_strm << ibmFI[0].at(i)  << "\t"; } else {res_strm << "NaN\t"; }
	if( ibmSP[0].at(i) > 0 ){ res_strm << ibmSP[0].at(i)  << "\t"; } else {res_strm << "NaN\t"; }
	if( ibmDS[0].at(i) > 0 ){ res_strm << ibmDS[0].at(i)  << "\t"; } else {res_strm << "NaN\t"; }
	if( parz.hhFLAG  == 1 ){ res_strm << ibmFHH[0].at(i) <<"\t" ; } else {res_strm << "NaN\t"; }
	if( parz.hhFLAG  == 1 && ibmPHHH[0].at(i) > 0 ){ res_strm << ibmPHHH[0].at(i) <<"\t" ; } else {res_strm << "NaN\t"; }
	if( ibmIH[0].at(i) > 0 ){ res_strm << ibmIH[0].at(i)<<"\t" ; } else {res_strm << "NaN\t" ; }
	if( ibmPRec[0].at(i) > 0 ){ res_strm << ibmPRec[0].at(i); } else {res_strm << "NaN"; }
	res_strm << "\t";
	if( IRR[0].at(i) > 0 ){ res_strm << IRR[0].at(i); } else {res_strm << "NaN"; }
	res_strm << "\t";
	if( IRRa[0].at(i) > 0 ){ res_strm << IRRa[0].at(i); } else {res_strm << "NaN"; }
	res_strm << endl;

      }
      res_strm.close();
      lgf << "closing results2..." << endl;
    }


    lgf << "opening results file 3 ..." << endl;
    res_strm.open( (parz.wkdir + "results3.dat").c_str() );
    if(!res_strm){lgf << "results3 file failed to open..." << endl;} else {
      lgf << "results3 file opened..." << endl;
      for(int i = 0; i < notime; ++i){//t,FOI,ARI
  	res_strm <<  time.at(i) << "\t" << ibmARIA[0].at(i) << "\t"; 
	res_strm << ibmARI[0].at(i) << "\t";
	if( ibmPHH[0].at(i) >= 0 && parz.hhFLAG ){ res_strm << ibmPHH[0].at(i)<<"\t" ; } else {res_strm << "NaN\t"; }
	if( ibmFIP[0].at(i) > 0 ){ res_strm << ibmFIP[0].at(i)  << "\t"; } else {res_strm << "NaN\t"; }
	res_strm << strIf[0].at(i) <<"\t"<< strInf[0].at(i) <<"\t"<< strIpf[0].at(i);
	res_strm << "\t";
	if( ibmPHHH[0].at(i) >= 0 && parz.hhFLAG ){ res_strm << ibmPHHH[0].at(i)<<"\t" ; } else {res_strm << "NaN\t"; }
	if( parz.hhFLAG ){ res_strm << ibmPhh[0].at(i)<<"\t"<< ibmPHhh[0].at(i) ; }
	res_strm << endl;
      }
      res_strm.close();
      lgf << "closing results3..." << endl;
    }


    lgf << "opening results file 4 ..." << endl;
    res_strm.open( (parz.wkdir + "results4.dat").c_str() );
    if(!res_strm){lgf << "results4 file failed to open..." << endl;} else {
      lgf << "results4 file opened..." << endl;
      for(int i = 0; i < notime; ++i){//t,TB prev, notifs, Rx, Inc
	N = ( ibmS[0].at(i) + ibmL[0].at(i) + ibmD[0].at(i) + ibmT[0].at(i) );
  	res_strm <<  time.at(i) << "\t" << 1e+5*ibmD2[0].at(i) << "\t"; 
	res_strm << 1e+5*tbnotes[0].at(i) << "\t";
	res_strm << 1e+5*ibmT2[0].at(i) <<"\t"<< 1e+5*Inc[0].at(i)/N << "\t";
	// bung cd4 stuff here
	for( unsigned int k = 0; k < cd4hist.at(0).size(); ++k ){
	  res_strm << cd4hist.at(i).at(k) << "\t";
	}
	res_strm << endl;
      }
      res_strm.close();
      lgf << "closing results4..." << endl;
    }


    lgf << "opening results file 5 ..." << endl;
    res_strm.open(  (parz.wkdir + "results5.dat").c_str() );
    if(!res_strm){lgf << "results5 file failed to open..." << endl;} else {
      lgf << "results5 file opened..." << endl;
      for(int i = 0; i < notime; ++i){//t,HHPI,HHart,HHipt,p&v-ECF?
	N = ( ibmS[0].at(i) + ibmL[0].at(i) + ibmD[0].at(i) + ibmT[0].at(i) );
  	res_strm <<  time.at(i) << "\t" << ibmHHPI[0].at(i)/N << "\t";
  	res_strm <<  ibmHHart[0].at(i)/N << "\t" << ibmHHipt[0].at(i)/N << "\t";
  	res_strm <<  ibmpECF[0].at(i)/N << "\t";
	if( ibmvECF[0].at(i) > 0 ){ res_strm << ibmvECF[0].at(i); } else { res_strm << "NaN"; }
  	res_strm << "\t" <<ibmARTp[0].at(i) << "\t" << ibmARTe[0].at(i) << "\t";
  	res_strm <<TBdeaths[0].at(i) <<"\t"<<HIVdeaths[0].at(i) <<"\t"<<TBHIVdeaths[0].at(i) <<"\t";
	res_strm << strIn[0].at(i) <<"\t"<< strIp[0].at(i) <<"\t"<< strIa[0].at(i)<<"\t"<<ibmSHIVcount[0].at(i);
	res_strm <<"\t"<< IPTprev[0].at(i) <<"\t"<< IPTprevP[0].at(i) <<"\t"<< IPTprevA[0].at(i) <<"\t" << IPTprevN[0].at(i)<<"\t" << IPTch[0].at(i);
	res_strm << endl;
      }
      res_strm.close();
      lgf << "closing results5..." << endl;
    }



    bool bigSS = parz.bigSS;		// switch for turning this off
    if( bigSS ){			// can be very big for many runs
      lgf << "opening population snapshot ..." << endl;
      res_strm.open(  (parz.wkdir + "PopSnapshot.csv").c_str() );
      if(!res_strm){lgf << "population snapshot file failed to open..." << endl;} else {
	lgf << " file opened population snapshot..." << endl;
	res_strm <<"age,acat,gender,TBu,TBRx,smrp,HIV,ART,HHHIV,phhi,tecf,LTBI,hhsize,hhid,cd4,isX,artdefr,IPT,pIPT,dTB,PI,tbrate" << endl;
	for(unsigned int i = 0; i < ibmstate.size(); ++i){//age, acat, gender, TBu, TBRx, smrp,HIV, ART, HHHIV,phhi,tecf
	  for(unsigned int j = 0; j < ibmstate.at(i).size(); ++j ){
	    res_strm << ibmstate.at(i).at(j);
	    if( j < (ibmstate.at(i).size()-1) ){ // not the last in row
	      res_strm << ",";		     // append c in csv
	    }
	  }
	  res_strm << endl;		// new line
	}
	res_strm.close();
	lgf << "closing population snapshot..." << endl;
      } // else-open
    } else {lgf << "big snap-shot turned off..." << endl;}

  } else {
    lgf << "Heavy results turned off!..." << endl;
  }

  return j;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int results::writedss( ofstream& lgf, int i, parameters& parz ){
  //DEMOGRAPHIC SNAPSHOT
  int j(0);
  ofstream res_strm;//the result file text
  lgf << "\topening demography results..." << endl;
  ostringstream oss;
  oss << parz.wkdir + "dss";
  oss << i;
  oss << ".dat";
  string name( oss.str() );
  res_strm.open( name.c_str() );
  if(!res_strm){lgf << "**demog results file failed to open...**" << endl;} else {
    lgf << "\tdemog results file opened..." << endl;
    res_strm << ";max hh size=" << Tmax << endl;
    res_strm << ";max n size = " << Wmax << endl;
    res_strm << ";max m size = " << Mmax << endl;
    res_strm << ";min hh size = " << Tmin << endl;
    res_strm << ";min n size = " << Wmin << endl;
    res_strm << ";min m size = " << Mmin << endl;
    res_strm << ";hhs of size 0 = " << hh0tot << endl;
    res_strm << ";hhs in total = " << hhtot << endl;

    double mean(0);
    for ( unsigned int a = 0; a < P.size(); ++a ){
      for( unsigned int b = 0; b < P.at(a).size(); ++b ){
	res_strm << a << ","<< b <<","<< P.at(a).at(b) / (hhtot - hh0tot) << endl;//mean excluding hhs of size 0
	mean += (a+b)*P.at(a).at(b) / (hhtot - hh0tot);
      }
    }
    res_strm.close();
    lgf << "mean no/hh (excl0hhs) = " << mean << endl;
    lgf << "\tclosing demog results..." << endl;
    j += 1;
  }
  return j;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//end results member funs 
