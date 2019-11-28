//#######################
#include "parameters.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PARAMETER MEMBER FUNS
//THIS CLASS WILL INCLUDE GENERALIZED PARAMETERS: HIV INTERPOLATION FUNS, BASELINE HOUSEHOLD STRUCTURE, HIV DEATH HAZARDS ETC
int parameters::computerest( void ){
  // set baseline cd4
  cd4bl = cd40;//1000;
  cd41 = 0.75 * cd40;

  IPTcovP1 = IPTcovP;		// the first coverage from the st

  if( DHflag ){			// sexual activity fractions
    sbFrac1 = sbFracH;
    sbFrac2 = sbFracM * (1-sbFracH);
    sbFrac3 = (1-sbFracM) * (1-sbFracH);
  }


  // divide by 100...
  artDR /= 100;
  artCU /= 100;

  // set time!
  time = starttime;

  //for computing the dependent parms
  int i = 0;//testing success
  popsize0 = S0 + L0 + D0 + T0;//initial population size
  notime = (int) ceil( (stoptime - starttime)/dt );//number of times

  //max of two effects for when both apply
  maxdOR = ( hhdOR > tecfdOR ? hhdOR : tecfdOR );
  maxdF = ( hhdF < tecfdF ? hhdF : tecfdF );//NB this is actually the min

  if( hhFLAG ){//depracated - now in hhpostinit...
     //mean household size
     nophh0 = 0;
//      double temptot(0);
//      for( unsigned int ui = 0; ui < hhtotsizes.size(); ++ui ){
//        nophh0 += hhtotsizes.at(ui) * ui;
//        temptot += hhtotsizes.at(ui);
//      } 
//      nophh0 /= temptot + 1e-10;
//      cout << nophh0 << " /HH!" << endl;
  }

   if( AGEflag ){
     nbins = agedata.size();//set the number of age categories...
   } else {
     nbins = 0;
   }

  nohh0 = auxFloor( 1.05 * popsize0 / nophh0 );//nohh0

  //WARN of inconsistency
  if( ( HIPflag || HARflag ) ){//HH sub-intervention on
    cout << "HIPflag,HARflag = " << HIPflag << "," <<HARflag << endl;
    if (HHIflag == 0 || hhFLAG == 0 ){//no hh or hh inter
      cout <<"BAILING! inconsistent HH interventions (1)!" << endl;
      exit(-1);
    }
  } 
  if ( !hhFLAG && HHIflag ){//no households, but HHI interventions
    cout << "BAILING! inconsistent HH interventions (2)!" << endl;
    exit(-1);
  } 
  if( !hhFLAG && ART2HHhx ){	// no households kind of HHI
    cout << "BAILING! inconsistent HH interventions (3)!" << endl;
    exit(-1);
  }
  if( !HARflag && ART2HHhx ){	// enhanced but not basic HH
    cout << "BAILING! advanced HH art w/o basic!" << endl;
    exit(-1);
  }
  
  //HIV incidence other parameters 
  h_alpha = h_peakiness;
  if( h_alpha <= 1 ){cout << "HIV incidence peakiness needs to be >1 !" << endl; exit(-1);}
  h_beta = ( h_peaktime - h_t0 ) / ( h_alpha - 1 );

  // normalize life data
  for( unsigned int jj = 1; jj < UNLTb[0].size(); ++jj ){
    double temp = UNLTb.at(1).at(jj); // the top element
    for( unsigned int ii = 1; ii < UNLTb.size(); ++ii ){
      UNLTb.at(ii).at(jj) /= temp; // normalize
    }
  }


  i = 1;//??change
  //  if ( (popsize0 > 0) && ( notime > 0) && ( deltan > 0 ) && ( deltap > 0 ) ){i += 1;}
  return i; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int parameters::read( ofstream& lgf, string& wdir ){//reading parameters in
  //READING IN PARAMETERS
  string line;
  istringstream ss;
  ifstream parfile;
  int j = 0;
  wkdir = wdir;			   // the parameter working directory
  parfile.open( (wkdir + "parz.dat").c_str() ); // convert to c-string and open
  if(!parfile){
    lgf << "failed to open parameter file...(bail)" << endl;
    exit(-1);
  } else {
     
    //run parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> noruns >> starttime >> stoptime >> dt;
    ss >> S0 >> L0 >> D0 >> T0;//ICs   
    lgf << "run parameters read:" << endl;
    lgf << "\tparz.noruns = " << noruns << endl;
    lgf << "\tparz.starttime = " << starttime << endl;
    lgf << "\tparz.stoptime = " << stoptime << endl;
    lgf << "\tparz.dt = " << dt << endl;
    lgf << "\tparz.S0 = " << S0 << endl;
    lgf << "\tparz.L0 = " << L0 << endl;
    lgf << "\tparz.D0 = " << D0 << endl;
    lgf << "\tparz.T0 = "<< T0 << endl;


    //population parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> fertage_alpha >> fertage_beta >> mf_deathratio >> AGEflag;
    lgf << "population parameters read:" << endl;
    lgf << "\tparz.fertage_alpha =" << fertage_alpha << endl;
    lgf << "\tparz.fertage_beta =" << fertage_beta << endl;
    lgf << "\tparz.mf_deathratio =" << mf_deathratio << endl;
    lgf << "\tparz.AGEflag =" << AGEflag << endl;


    //non-HIV TB parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> beta >>  tbsurvn_k >> tbsurvn_L >> tbmortn >>  r0 >> r1 >> smrfac >> fscale >> smrnmort >> cdrFN;
    lgf << "non-HIV TB parameters read in:" << endl;
    lgf << "\tparz.beta = " <<  beta << endl;
    lgf << "\tparz.tbsurvn_k = " <<  tbsurvn_k  << endl;
    lgf << "\tparz.tbsurvn_L = " <<  tbsurvn_L << endl;
    lgf << "\tparz.tbmortn = " <<  tbmortn << endl;
    lgf << "\tparz.r0 = " <<  r0 << endl;
    lgf << "\tparz.r1 = " <<  r1 << endl;
    lgf << "\tparz.smrfac = " <<  smrfac << endl;
    lgf << "\tparz.fscale = " <<  fscale << endl;
    lgf << "\tparz.smrnmort = " <<  smrnmort << endl;
    lgf << "\tparz.cdrFN = " <<  cdrFN << endl;

    //more non-HIV TB 
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> CDRn >> RxL >> Rxp >> Rxpd >> Iprotn >> IprotnS >> Pprotn >> PprotnS;
    ss >> fSmr >> deltanSmr >> tbd_betan;
    lgf << "\tparz.CDRn = " << CDRn << endl;
    lgf << "\tparz.RxL = " << RxL << endl;
    lgf << "\tparz.Rxp = " << Rxp << endl;
    lgf << "\tparz.Rxpd = " << Rxpd << endl;
    lgf << "unused\tparz.Iprotn = " << Iprotn << endl;
    lgf << "unused\tparz.IprotnS = " << IprotnS << endl;
    lgf << "\tparz.Pprotn = " << Pprotn << endl;
    lgf << "\tparz.PprotnS = " << PprotnS << endl;
    lgf << "\tparz.fSmr = " << fSmr << endl;
    lgf << "\tparz.deltanSmr = " << deltanSmr << endl; 
    lgf << "\tparz.tbd_betan = " << tbd_betan << endl;


    //HIV parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> thiv_alpha >> thiv_beta >> mhiv_alpha >> mhiv_beta >> fhiv_alpha >> fhiv_beta >> fm_hivratio >> cd40;
    lgf << "HIV parameters read:" << endl;
    lgf << "\tparz.thiv_alpha = "<< thiv_alpha << endl;
    lgf << "\tparz.thiv_beta = "<< thiv_beta << endl;
    lgf << "\tparz.mhiv_alpha = "<< mhiv_alpha << endl;
    lgf << "\tparz.mhiv_beta = "<< mhiv_beta << endl;
    lgf << "\tparz.fhiv_alpha = " << fhiv_alpha << endl;
    lgf << "\tparz.fhiv_beta = " << fhiv_beta << endl;
    lgf << "\tparz.fm_hivratio =" << fm_hivratio << endl;
    lgf << "\tparz.cd40 =" << cd40 << endl;
    
    //HIV incidence parameters:
    getline( parfile, line );//comments
    getline( parfile, line );
    ss.clear();
    ss.str(line);
    ss >> h_t0 >> h_peak >> h_peakiness >> h_peaktime >> h_theta >> decline2010;
    ss >> hivincF >> hivincFN; 
    lgf << "HIV incidence parameters read:" << endl;
    lgf << "\tparz.h_t0 = "<< h_t0 << endl;
    lgf << "\tparz.h_peak = "<< h_peak << endl;
    lgf << "\tparz.h_peakiness = " << h_peakiness << endl;
    lgf << "\tparz.h_peaktime = " << h_peaktime << endl;
    lgf << "\tparz.h_theta = "<< h_theta << endl;
    lgf << "\tparz.decline2010 = "<< decline2010  << endl;
    lgf << "\tparz.hivincF = "<< hivincF  << endl;
    lgf << "\tparz.hivincFN = "<< hivincFN  << endl;

    //HIV-TB parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss  >> hivdurf >> f >> g >> rho >> CDRp >> Rxpp >> Rxppd;
    ss >> deltapSmr >>  tbd_betap >> H2Lhr;
    lgf << "HIV/TB parameters read:" << endl;
    //    lgf << "\tparz.probSCp = "<< probSCp << endl;
    lgf << "\tparz.hivdurf = "<< hivdurf << endl;
    lgf << "\tparz.f = " << f << endl;
    lgf << "\tparz.g = " << g << endl;
    lgf << "\tparz.rho = "<< rho << endl;
    lgf << "\tparz.CDRp = "<< CDRp << endl;
    lgf << "\tparz.Rxpp = " << Rxpp << endl;
    lgf << "\tparz.Rxppd = " << Rxppd << endl;
    lgf << "\tparz.deltapSmr = " << deltapSmr << endl;
    lgf << "\tparz.tbd_betap = " << tbd_betap << endl;
    lgf << "\tparz.H2Lhr = " << H2Lhr << endl;



    //household parameters
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> hhFLAG >> nophh0 >> hhmoverate >> betaH >> hhhivRR;
    lgf << "household parameters read:" << endl;
    lgf << "\tparz.hhFLAG = " << hhFLAG << endl;
    lgf << "\tparz.nophh0 = " << nophh0 << endl;
    lgf << "\tparz.hhmoverate = " << hhmoverate << endl;
    lgf << "\tparz.betaH = " << betaH << endl;
    lgf << "\tparz.hhhivRR = " << hhhivRR << endl;

    //intervention parameters
    getline( parfile, line);//comments
    getline( parfile, line);//ECF
    ss.clear();
    ss.str(line);
    ss >> ECFflag >> ECFst >> ECFet >> CDRn2 >> CDRp2;
    ss >> tbd_betan2  >> tbd_betap2;
    getline( parfile, line);//comments
    getline( parfile, line);//Rx improvement
    ss.clear();
    ss.str(line);
    ss >> Rxflag >> Rxtime >> Rxp2 >> Rxpp2;
    getline( parfile, line);//comments
    getline( parfile, line);//IPT interventions
    ss.clear();
    ss.str(line);
    ss >> IPTflag >> IPTst >> IPTet >> IPTcov >> IPTT;
    ss >> IPTart >> IPTshiv >> IPTcovP;
    ss >> IPTcovP2 >> IPTcovPt1 >> IPTcovPt2 >> TSTspec >> TSTsens;
    getline( parfile, line);//comments
    getline( parfile, line);//IPT characteristics
    ss.clear();
    ss.str(line);
    ss >> IPTpprot >> IPTdurnN;
    ss >> IPTdurnP >> IPThrN >> IPThrP >> IPThrA;
    ss >> IPTcprobN >> IPTcprobP >> IPTmultiple;
    getline( parfile, line);//comments
    getline( parfile, line);//ART
    ss.clear();
    ss.str(line);
    ss >> ARTflag >> ARTst >> ARTrt >> ARTcovmax  >> tbart >> cd4h >> cd4g >> cd4g2 >> ARTf >> cd42a >> cd42b >> cd42c >> expET >> artTk >> artTl >> dCD4st >> dCD4et >> dCD4ep >> art750 >> artDR >> artCU;
    getline( parfile, line);//comments
    getline( parfile, line);//HHI
    ss.clear();
    ss.str(line);
    ss >> HHIflag >> HIPflag >> HARflag >> HHIst >> HHIet >> HHIcov >> TSTu16sens >> hhdOR >>hhdF >> deltaSmrI ;
    //and now the real ECF!
    getline( parfile, line);//comments
    getline( parfile, line);//tECF
    ss.clear();
    ss.str(line);
    ss >> tECFflag >> tECFst >> tECFet >> tecfdOR >> tecfdF >> tECFhaz ;//need coverage/rate


    lgf <<"\t\t---"<<endl;
    lgf << "intervention parameters read:" << endl;
    lgf << "\tparz.ECFflag = " << ECFflag << endl;
    lgf << "\tparz.ECFst = " << ECFst << endl;
    lgf << "\tparz.ECFet = " << ECFet << endl;
    lgf << "\tparz.CDRn2 = " << CDRn2 << endl;
    lgf << "\tparz.CDRp2 = " << CDRp2 << endl;
    lgf << "\tparz.tbd_betap2 = " << tbd_betap2 << endl;
    lgf << "\tparz.tbd_betan2 = " << tbd_betan2 << endl;

    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.Rxflag = " << Rxflag << endl;
    lgf << "\tparz.Rxtime = " << Rxtime << endl;
    lgf << "\tparz.Rxp2 = " << Rxp2 << endl;
    lgf << "\tparz.Rxpp2 = " << Rxpp2 << endl;

    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.IPTflag = " << IPTflag << endl;
    lgf << "\tparz.IPTst = " << IPTst << endl;
    lgf << "\tparz.IPTet = " << IPTet << endl;
    lgf << "\tparz.IPTcov = " << IPTcov << endl;
    lgf << "\tparz.IPTT = " << IPTT << endl;
    lgf << "\tparz.IPTart = " << IPTart << endl;
    lgf << "\tparz.IPTshiv = " << IPTshiv << endl;
    lgf << "\tparz.IPTcovP = " << IPTcovP << endl;
    lgf << "\tparz.IPTcovP2 = " << IPTcovP2 << endl;
    lgf << "\tparz.IPTcovPt1 = " << IPTcovPt1 << endl; 
    lgf << "\tparz.IPTcovPt2 = " << IPTcovPt2 << endl;
    lgf << "\tparz.TSTspec = " << TSTspec << endl;
    lgf << "\tparz.TSTsens = " << TSTsens << endl;

    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.IPTpprot = " << IPTpprot << endl;
    lgf << "\tparz.IPTdurnN = " << IPTdurnN << endl;
    lgf << "\tparz.IPTdurnP = " << IPTdurnP << endl;
    lgf << "\tparz.IPThrN = " << IPThrN << endl;
    lgf << "\tparz.IPThrP = " << IPThrP << endl;
    lgf << "\tparz.IPThrA = " << IPThrA << endl;
    lgf << "\tparz.IPTcprobN = " << IPTcprobN << endl;
    lgf << "\tparz.IPTcprobP = " << IPTcprobP << endl;
    lgf << "\tparz.IPTmultiple = " << IPTmultiple << endl;


    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.ARTflag = " << ARTflag << endl;
    lgf << "\tparz.ARTst = " << ARTst << endl;
    lgf << "\tparz.ARTrt = " << ARTrt << endl;
    lgf << "\tparz.ARTcovmax = " << ARTcovmax << endl;
    lgf << "\tparz.tbart = " << tbart << endl;
    lgf << "\tparz.cd4h = " << cd4h << endl;
    lgf << "\tparz.cd4g = " << cd4g << endl;
    lgf << "\tparz.cd4g2 = " << cd4g2 << endl;
    lgf << "\tparz.ARTf = " << ARTf << endl;
    lgf << "\tparz.cd42a = " << cd42a << endl;
    lgf << "\tparz.cd42b = " << cd42b << endl;
    lgf << "\tparz.cd42c = " << cd42c << endl;
    lgf << "\tparz.expET = " << expET << endl;
    lgf << "\tparz.artTk = " << artTk << endl;
    lgf << "\tparz.artTl = " << artTl << endl;
    lgf << "\tparz.dCD4st = " << dCD4st << endl;
    lgf << "\tparz.dCD4et = " << dCD4et << endl;
    lgf << "\tparz.dCD4ep = " << dCD4ep << endl;
    lgf << "\tparz.art750 = " << art750 << endl;
    lgf << "\tparz.artDR = " << artDR << endl;
    lgf << "\tparz.artCU = " << artCU << endl;

    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.HHIflag = " << HHIflag << endl;
    lgf << "\tparz.HIPflag = " << HIPflag << endl;
    lgf << "\tparz.HARflag = " << HARflag << endl;
    lgf << "\tparz.HHIst = " << HHIst << endl;
    lgf << "\tparz.HHIet = " << HHIet << endl;
    lgf << "\tparz.HHIcov = " << HHIcov << endl;
    lgf << "unused\tparz.TSTu16sens = " << TSTu16sens << endl;
    lgf << "\tparz.hhdOR = " << hhdOR << endl;
    lgf << "\tparz.hhdF = " << hhdF << endl;
    lgf << "\tparz.deltaSmrI = " << deltaSmrI << endl;

    lgf <<"\t\t---"<<endl;
    lgf << "\tparz.tECFflag = " << tECFflag << endl;
    lgf << "\tparz.tECFst = " << tECFst << endl;
    lgf << "\tparz.tECFet = " << tECFet << endl;
    lgf << "\tparz.tecfdOR = " << tecfdOR << endl;
    lgf << "\tparz.tecfdF = " << tecfdF << endl;
    lgf << "\tparz.tECFhaz = " << tECFhaz << endl;

    //file names for external data
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> birthFN >> deathFN >> migrateFN >> amFN;
    ss >> synpop0 >> LTbFN >> abFN >> rpopFN;
    lgf << "file names for external data read:" << endl;
    lgf << "\tparz.birthFN = " << birthFN << endl;
    lgf << "\tparz.deathFN = " << deathFN << endl;
    lgf << "\tparz.migrateFN = " << migrateFN << endl;
    lgf << "\tparz.amFN = " << amFN << endl;
    lgf << "\tparz.synpop0 = " << synpop0 << endl;
    lgf << "\tparz.LTbFN = " << LTbFN << endl;
    lgf << "\tparz.abFN = "<< abFN  << endl;
    lgf << "\tparz.rpopFN = "<< rpopFN  << endl;

    //other intervention parms
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> ART2TB >> ART2TBst >> ART2HHhx >> ARThxT  >> shivt >> shivp;
    lgf << "extra interventions:" << endl;
    lgf << "\tparz.ART2TB = " << ART2TB << endl;
    lgf << "\tparz.ART2TBst = " << ART2TBst << endl;
    lgf << "unused\tparz.ART2HHhx = " << ART2HHhx << endl;
    lgf << "unused\tparz.ARThxT = " << ARThxT << endl;
    lgf << "\tparz.shivt = " << shivt << endl;
    lgf << "\tparz.shivp = " << shivp << endl;
    // will needs some parms to do with HIV self-knowlege/detection?

    //mass ART
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> massart >> artstart >> artT >> artP >> partrefusenik;
    lgf << "extra interventions:" << endl;
    lgf << "\tparz.massart = " << massart << endl;
    lgf << "\tparz.artstart = " << artstart << endl;
    lgf << "\tparz.artT = " << artT << endl;
    lgf << "\tparz.artP = " << artP << endl;
    lgf << "\tparz.partrefusenik = " << partrefusenik << endl;

    //periodic HIV screening
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> HIVflag >> HIVst >> HIVet >> HIVT >> HIVeff;
    lgf << "periodic HIV screening:" << endl;
    lgf << "\tparz.HIVflag = " << HIVflag << endl;
    lgf << "\tparz.HIVst = " << HIVst << endl;
    lgf << "\tparz.HIVet = " << HIVet << endl;
    lgf << "\tparz.HIVT = " << HIVT << endl;
    lgf << "\tparz.HIVeff = " << HIVeff << endl;


    //periodic ACF
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> PACFlag>> PACFst >> PACFet >> PACFT >> PACFeffsp >> PACFeffsn >> PACFeffart >>PACFhivOR;
    lgf << "periodic ACF interventions:" << endl;
    lgf << "\tparz.PACFlag = " << PACFlag << endl;
    lgf << "\tparz.PACFst = " << PACFst << endl;
    lgf << "\tparz.PACFet = " << PACFet << endl;
    lgf << "\tparz.PACFT = " << PACFT << endl;
    lgf << "\tparz.PACFeffsp = " << PACFeffsp << endl;
    lgf << "\tparz.PACFeffsn = " << PACFeffsn << endl;
    lgf << "\tparz.PACFeffart = " << PACFeffart << endl;
    lgf << "\tparz.PACFhivOR = " << PACFhivOR << endl;


    // resistant strains 
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> Rflag >> Rfrac0;
    lgf << "resistance:" << endl;
    lgf << "\tparz.Rflag = " << Rflag << endl;
    lgf << "\tparz.Rfrac0 = " << Rfrac0 << endl;


    // dynamic HIV
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> DHflag >> sbFracH >> sbFracM;
    ss >> sbc1 >> sbc2 >> sbc3 >> sba1 >> sba2 >> sba3 >> sbA;
    ss >> lam >> riA >> riE >> durA >> durE >> ARTeps >> hiv0 ;
    lgf << "dynamic HIV:" << endl;
    lgf << "\tparz.DHflag = " << DHflag << endl;
    if( DHflag ){
      lgf << "\tparz.sbFracH = " << sbFracH << endl;
      lgf << "\tparz.scFracM = " << sbFracM << endl;
      lgf << "\tparz.sbc1 = " <<  sbc1 << endl;
      lgf << "\tparz.sbc2 = " <<  sbc2 << endl;
      lgf << "\tparz.sbc3 = " <<  sbc3 << endl;
      lgf << "\tparz.sba1 = " <<  sba1 << endl;
      lgf << "\tparz.sba2 = " <<  sba2 << endl;
      lgf << "\tparz.sba3 = " <<  sba3 << endl;
      lgf << "\tparz.sbA = " <<  sbA << endl;
      lgf << "\tparz.lam = " <<  lam << endl;
      lgf << "\tparz.riA = " <<  riA << endl;
      lgf << "\tparz.riE = " <<  riE << endl;
      lgf << "\tparz.durA = " <<  durA << endl;
      lgf << "\tparz.durE = " <<  durE << endl;
      lgf << "\tparz.ARTeps = " <<  ARTeps << endl;
      lgf << "\tparz.hiv0 = " <<  hiv0 << endl;
    }

    // HIVMAC methods
    // artsig, pscale, pn, kq, gg
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> hivmac >> artsig >> pscale >> pn >> pmaxval >> kq >> gg >> gi >> glrst;
    ss >> nglstr >> hseed >> SQ >> SQGLprop;
    lgf << "HIVMAC extras:" << endl;
    lgf << "\tparz.hivmac = " << hivmac << endl;
    lgf << "\tparz.artsig = " << artsig << endl;
    lgf << "\tparz.pscale = " << pscale << endl;
    lgf << "\tparz.pn = " << pn << endl;
    lgf << "\tparz.pmaxval = " << pmaxval << endl;
    lgf << "\tparz.kq = " << kq << endl;
    lgf << "\tparz.gg = " << gg << endl;
    lgf << "\tparz.gi = " << gi << endl;
    lgf << "\tparz.glrst = " << glrst << endl;
    lgf << "\tparz.nglstr = " << nglstr << endl;
    lgf << "\tparz.hseed = " << hseed << endl;
    lgf << "\tparz.SQ = " << SQ << endl;
    lgf << "\tparz.SQGLprop = " << SQGLprop << endl;


    // switches 
    getline( parfile, line);//comments
    getline( parfile, line);
    ss.clear();
    ss.str(line);
    ss >> qlog >> lightFLAG >> snaps >> bigSS;
    lgf << "switches:" << endl;
    lgf << "\tparz.qlog = " << qlog << endl;
    lgf << "\tparz.lightFLAG = " << lightFLAG << endl;
    lgf << "\tparz.snaps = " << snaps << endl;
    lgf << "\tparz.bigSS = " << bigSS << endl;

    j += 1;
  }
  parfile.close();



  //READING IN BIRTH/DEATH/MIGRATION DATA
  double tempa;
  string fname;
  // fname = wkdir + birthFN;
  //fname = "data/" + birthFN;
  fname = birthFN;
  //births
  ifstream popfile ( fname.c_str() );//
  //  popfile.open ("data/birthrate.dat");
  if ( popfile.is_open() ){
    getline( popfile, line );//yrz
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){br_yrz.push_back( tempa );}
    getline( popfile, line );//data
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){br_dat.push_back( tempa );}
    popfile.close();
    lgf << "birthrates read..." << endl;
    if ( br_dat.size() != br_yrz.size() ){ lgf << "**birth-rate misread!**" << endl;} else {j += 1;}
  } else {
    lgf << "**failed to open birthrate.dat!**" << endl;
    cerr << "warning: failed to find birthrate data!" << endl;
    exit(-1);
  }
  //precomputing splines...
  lgf << "precomputing br splines..." << endl;
  br_spl.assign( (br_dat.size()-1), 0 );//to give it length
  splinetable( br_yrz, br_dat, br_spl );
  lgf <<"\tdone!" << endl;


  //deaths
  // fname = wkdir + deathFN;
  //fname = "data/" + deathFN;
  fname = deathFN;
  popfile.open ( fname.c_str() );
  if ( popfile.is_open() ){
    getline( popfile, line );//yrz
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){dr_yrz.push_back( tempa );}
    getline( popfile, line );//data
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){dr_dat.push_back( tempa );}
    popfile.close();
    lgf << "deathrates read..." << endl;
    if ( dr_dat.size() != dr_yrz.size() ){ lgf << "**death-rate misread!**" << endl;} else {j += 1;}
  } else {
    lgf << "**failed to open deathrate.dat!**" << endl;
    cerr << "warning: failed to find deathrate data!" << endl;
    exit(-1);
  }
  //precomputing splines...
  lgf << "precomputing dr splines..." << endl;
  dr_spl.assign( (dr_dat.size()-1), 0 );//to give it length
  splinetable( dr_yrz, dr_dat, dr_spl );
  lgf <<"\tdone!" << endl;

  // HIV incidence if read from file
  if ( hivincF ){		// yes, using file rather than parameters
    fname = hivincFN;
    popfile.open ( fname.c_str() );
    if ( popfile.is_open() ){
      getline( popfile, line );//yrz
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){h1549_yrz.push_back( tempa );}
      getline( popfile, line );//data
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){h1549_dat.push_back( tempa );}
      popfile.close();
      lgf << "hiv incidence read..." << endl;
      if ( h1549_dat.size() != h1549_yrz.size() ){ lgf << "**hiv inc misread!**" << endl;} else {j += 1;}
    } else {
      lgf << "**failed to open hiv incidence data!**" << endl;
      cerr << "warning: failed to find hiv incidence data!" << endl;
      exit(-1);
    }
    //precomputing splines...
    lgf << "precomputing hiv incidence splines..." << endl;
    h1549_spl.assign( (h1549_dat.size()-1), 0 );//to give it length
    splinetable( h1549_yrz, h1549_dat, h1549_spl );
    lgf <<"\tdone!" << endl;
  }



  // CDR timeseries if one
  if ( cdrFN != "none" ){      	// yes, using file rather than parameters
    fname = cdrFN;
    popfile.open ( fname.c_str() );
    if ( popfile.is_open() ){
      getline( popfile, line );//yrz
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){cdr_yrz.push_back( tempa );}
      getline( popfile, line );//data
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){cdr_dat.push_back( tempa );}
      popfile.close();
      lgf << "CDR timeseries read..." << endl;
      if ( cdr_dat.size() != cdr_yrz.size() ){ lgf << "**CDR misread!**" << endl;} else {j += 1;}
    } else {
      lgf << "**failed to open CDR data!**" << endl;
      cerr << "warning: failed to find CDR  data!" << endl;
      exit(-1);
    }
    //precomputing splines...
    lgf << "precomputing CDR splines..." << endl;
    cdr_spl.assign( (cdr_dat.size()-1), 0 );//to give it length
    splinetable( cdr_yrz, cdr_dat, cdr_spl );
    lgf <<"\tdone!" << endl;
  }

  // relpop timeseries if one
  if ( rpopFN != "none" ){       // yes, using file rather than parameters
    fname = rpopFN;
    popfile.open ( fname.c_str() );
    if ( popfile.is_open() ){
      getline( popfile, line );//yrz
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){rpop_yrz.push_back( tempa );}
      getline( popfile, line );//data
      ss.clear();
      ss.str(line);
      while ( ss >> tempa ){rpop_dat.push_back( tempa );}
      popfile.close();
      lgf << "relative population timeseries read..." << endl;
      if ( cdr_dat.size() != rpop_yrz.size() ){ lgf << "**rel pop misread!**" << endl;} else {j += 1;}
    } else {
      lgf << "**failed to open rel pop data!**" << endl;
      cerr << "warning: failed to find rel pop data!" << endl;
      exit(-1);
    }
    //precomputing splines...
    lgf << "precomputing relative population splines..." << endl;
    rpop_spl.assign( (rpop_dat.size()-1), 0 );//to give it length
    splinetable( rpop_yrz, rpop_dat, rpop_spl );
    lgf <<"\tdone!" << endl;
  }



//   //DEBUG!

  // lgf << "SPLhiv:" << endl;
  // for(unsigned int ti = 0; ti < h1549_yrz.size(); ++ti){lgf <<h1549_yrz.at(ti)<<",";} lgf << endl;
  // for(unsigned int ti = 0; ti < h1549_dat.size(); ++ti){lgf <<h1549_dat.at(ti)<<",";} lgf << endl;
  // for(unsigned int ti = 0; ti < h1549_spl.size(); ++ti){lgf <<h1549_spl.at(ti)<<",";} lgf << endl;
  // cout << h1549_yrz.back()<< "|" << h1549_dat.back() <<"|" <<h1549_spl.back()<<endl;
  // exit(-1);

  //migrations
  // fname = wkdir + migrateFN;
  //fname = "data/" + migrateFN;
  fname = migrateFN;
  popfile.open ( fname.c_str() );
  // popfile.open ("data/migrate.dat");
  if ( popfile.is_open() ){
    getline( popfile, line );//yrz
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){mr_yrz.push_back( tempa );}
    getline( popfile, line );//data
    ss.clear();
    ss.str(line);
    while ( ss >> tempa ){mr_dat.push_back( tempa );}
    popfile.close();
    lgf << "migration rates read..." << endl;
    if ( mr_dat.size() != mr_yrz.size() ){ lgf << "**migration-rate misread!**" << endl;} else {j += 1;}
  } else {
    lgf << "**failed to open migrate.dat!**" << endl;
    cerr << "warning: failed to find migration rate data!" << endl;
    exit(-1);
  }
  //precomputing splines...
  lgf << "precomputing mr splines..." << endl;
  mr_spl.assign( (mr_dat.size()-1), 0 );//to give it length
  splinetable( mr_yrz, mr_dat, mr_spl );
  lgf <<"\tdone!" << endl;

  //switched to salomon murray

  // reading in lifetable
  vector< double > tempvec;//this is used for reading hh, age, and LT data
  string temps;
  fname = LTbFN;
  popfile.open ( fname.c_str() );
  if ( popfile.is_open() ){
    while( getline( popfile, line ) ){//life
      ss.clear();
      tempvec.clear();
      ss.str(line);
      while ( ss >> tempa ){tempvec.push_back( tempa );}
      UNLTb.push_back( tempvec );
    }
    popfile.close();
  } else {
    lgf << "**failed to open lifetable!**" << endl;
    cerr << "warning: failed to find life table data!" << endl;
  }



  // cout<<"LT:"<<UNLTb.size()<<","<<UNLTb.at(0).size()<<endl;
  // for(unsigned int li=0; li < UNLTb.size();++li){
  //   for(unsigned int lj=0; lj < UNLTb.at(li).size();++lj){
  //     cout << UNLTb.at(li).at(lj)<<",";
  //   }
  //   cout<<endl;
  // }

  
  //pretty sure this not needed when reading in a population
//   if( hhFLAG ){
//   //READING IN HOUSEHOLD SIZES...
//   lgf << "reading household sizes..." << endl;
//   popfile.open ("data/nmcounts.dat");
//   if ( popfile.is_open() ){
//     while( getline( popfile, line ) ){//hhz
//       ss.clear();
//       tempvec.clear();
//       ss.str(line);
//       while ( ss >> tempa ){tempvec.push_back( tempa );}
//       hhsizes.push_back( tempvec );
//     }
//     popfile.close();
//   } else {
//     lgf << "**failed to open nmcounts.dat!**" << endl;
//     cerr << "warning: failed to find household size data!" << endl;
//      tempvec.assign( 10, 1 ); 
//      for ( unsigned int ii = 0; ii < 10; ++ ii){//even distribution
//        hhsizes.push_back( tempvec );
//     }
//   }

//   hhtotsizes.push_back(0);
//   //renormalize
//   nohh0 = 0;
//   for( unsigned int ii = 0; ii < hhsizes.size(); ++ii ){
//     for (unsigned int jj = 0; jj < hhsizes.at(ii).size(); ++jj){
//       //nohh0 += hhsizes.at(ii).at(jj);
//       if( (ii+jj) != 0 ){//ignore 0	
// 	nohh0 += hhsizes.at(ii).at(jj);
// 	if ( ii+jj >= hhtotsizes.size() ){
// 	  hhtotsizes.push_back( hhsizes.at(ii).at(jj) );
// 	} else{
// 	  hhtotsizes.at( ii + jj ) += hhsizes.at(ii).at(jj);
// 	}
//       }//ignore the zeros
//     }
//   }
//   //lgf << nohh0 << " initial households..." << endl;
//   for( unsigned int ii = 0; ii < hhsizes.size(); ++ii ){
//     for (unsigned int jj = 0; jj < hhsizes.at(ii).size(); ++jj){
//      hhsizes.at(ii).at(jj) /= nohh0 ;
//     }
//   }
//   hhsizes.at(0).at(0) = 0;
//   for( unsigned int ii = 0; ii < hhtotsizes.size(); ++ii ){hhtotsizes.at(ii) /= nohh0;}
//   hhtotsizes.at(0) = 0;

//   //compute weights for alternative method...
//   //going to imagine buffer of 10% empty households
//   //NB the above are proportions...
//   f_weight.assign( hhtotsizes.size(), 0.0 );
//   a_weight.assign( hhtotsizes.size(), 0.0 );
//   f_weight.at(0) = hhtotsizes.at(1);// / 0.1 ;//buffer
//   double temp = 0.0; //buffer
//   a_weight.at(0) = 1 - temp;//buffer
//   for( unsigned int  ii = 1; ii < hhtotsizes.size()  - 1; ++ii ){
//     f_weight.at(ii) =  (ii+1) * hhtotsizes.at(ii+1) / ( hhtotsizes.at(ii) + 1e-10 );
//     if ( hhtotsizes.at(ii) < 1e-10 ){ f_weight.at(ii) = 0; }
//     temp += (1 - 0.0) * hhtotsizes.at(ii);//buffer
//     a_weight.at(ii) = ( 1 - temp )  /  ( hhtotsizes.at(ii) + 1e-10 );
//     if ( hhtotsizes.at(ii) < 1e-10 ){ a_weight.at(ii) = 0; }
//   }
//   f_weight.at( hhtotsizes.size()-1 ) = 0;
//   a_weight.at( hhtotsizes.size()-1 ) = 0;

//   lgf << "N = ";
//   for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
//     lgf << hhtotsizes.at(ii) <<",";
//   }
//   lgf << endl;

//   lgf << "F = ";
//   for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
//     lgf << f_weight.at(ii) <<",";
//   }
//   lgf << endl;

//   lgf << "A = ";
//   for( unsigned int  ii = 0; ii < hhtotsizes.size(); ++ii ){
//     lgf << a_weight.at(ii) <<",";
//   }
//   lgf << endl;

//   lgf <<"\tdone!" << endl;

//   }//end hhFLAG

  //READING IN THE AGE-MIXING DATA
  if ( !AGEflag ){
    lgf << "age-mixing off, agedata not read..." << endl;
  } else {

    lgf << "reading in the age-mixing data..." << endl;
    //reading age bins
    // fname = wkdir + "Abins.dat";
    //fname = "shared/Abins.dat";
    fname = abFN;
    popfile.open ( fname.c_str() );
    if ( popfile.is_open() ){
      	getline( popfile, line );
	ss.clear();
	tempvec.clear();
	ss.str(line);
	while ( ss >> tempa ){agebins.push_back( tempa );}
	popfile.close();
    } else {
      lgf << "**failed to open Abins.dat!**" << endl;
      cerr << "warning: failed to find age bin data!" << endl;
    }

    //reading age data
    // fname = wkdir + amFN;
    //fname = "data/" + amFN;
    fname = amFN;
    popfile.open ( fname.c_str() );
    //popfile.open ("data/AM.dat");
    if ( popfile.is_open() ){
      while( getline( popfile, line ) ){//agez
	ss.clear();
	tempvec.clear();
	ss.str(line);
	while ( ss >> tempa ){tempvec.push_back( tempa );}
	agedata.push_back( tempvec );
      }
      popfile.close();
    } else {
      lgf << "**failed to open AM.dat!**" << endl;
      cerr << "warning: failed to find age mixing data!" << endl;
    }

    //check matching...
    if( !( (agedata.size()==agedata.at(0).size()) && (agedata.size()==agebins.size()) )){
      lgf <<"problem with bins and age data sizes not matching!"<< endl;
      lgf << agedata.size() << "," << agedata.at(0).size() <<","<< agebins.size() << endl;
      exit(-1);
    } else {
      lgf << "\tage dimensions match ok..." << endl;
    }

  }//end age-data else

  // WARN about degenerate seeding due to too many runs!
  if( noruns >= 10000 ){
    cout << "Err...you need to rethink the seeding then! Bailing..." <<endl;
    lgf << "Err...you need to rethink the seeding then! Bailing..." <<endl;
    exit(-1);
  }


  return j;
}
////////////////////////////////////////////////////////////////////////////////////////////
int parameters::readpopdata( const gsl_rng *r, ofstream& lgf ){//read the data only

  int j = 0;
  string line, fname;
  istringstream ss;
  ifstream parfile;
  //fname = "data/" + synpop0;
  fname = synpop0;
  parfile.open( fname.c_str() );
  if(!parfile){
    lgf << "**failed to open population file...(bail)!**" << endl;
    exit(-1);
  } else {
    lgf << "reading in the baseline population state..." << endl;
    string name;
    int comit(1), hhbar(0), gend(0);//household bar code and gender
    double ag(0),agi(0);//age
    int chh(0);//current community, current household
    unsigned int totpopcount(0), tothhcount(0);//total counts
    unsigned int ui(0), chhcount(0);
    ss.clear(); 
    while( getline( parfile, line ) ){
      ss.str(line);
      
      if( ui > 0  ){//ignore the title line, && ui < 50000
	ss >> comit >> hhbar >> gend >> agi;// assign this person's variables

	//data cleaning
	//gender
	if( gend == 2 ){ gend = 0;}//convert gender coding
	if( !( gend == 0 || gend == 1) ){
	  lgf << "warning: fishy gender at " << ui << endl;
	  gend = 0;
	}
	//age
	if( !(ag >= 0 && ag < 100) ){
	  ag = - log(gsl_rng_uniform(r) / (1.0/50.0));
	  agi = (int)auxFloor(ag);
	  lgf << "warning: fishy age at " << ui << endl;
	}

	//store data
	vector< int > tmp; tmp.push_back( comit ); tmp.push_back( hhbar ); tmp.push_back( gend ); tmp.push_back( agi );
	popdata.push_back( tmp );

	// 	if( ui % 1000 == 0 ){
	// 	  lgf << "\t...ui  = "<< ui << endl;
	// 	  lgf << "\t...comit,hhbar,gend,ag  = "<<comit<<","<<hhbar<<","<<gend<<","<<ag<< endl;
	// 	}
	++totpopcount;


	if( chh != hhbar ){//a new household!
	  chh = hhbar;//adjust current household
	  ++chhcount;
	  ++tothhcount;
	  //loc = chhcount - 1;
	}
	
      }//end ignore title line
      
      ++ui;
      ss.clear();
    }
    parfile.close();

    j = (int) ui;

    //log progress
    lgf << " total population count = "<< totpopcount << endl;
    lgf << " total household count = "<< tothhcount << endl;

  }//end open file?

  lgf << "...population data read-in completed!" << endl;

  return j;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//birthrate interpolation wrapper
double parameters::br( double yr ){
  if ( yr < br_yrz.front() ){
     return br_dat.front() / 1000;
  } else if ( yr > br_yrz.back() ){
    return br_dat.back() / 1000;
  } else {
    return splinefun( yr, br_yrz, br_dat, br_spl ) / 1000;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//deathrate interpolation wrapper
double parameters::dr( double yr ){
  double ans(0);
  if ( yr < dr_yrz.front() ){
    ans = dr_dat.front();
  } else if ( yr > dr_yrz.back() ){
    ans = dr_dat.back();
  } else {
    ans = splinefun( yr, dr_yrz, dr_dat, dr_spl );
  }
  return ans + 10;		
  // this is a bodge to cope with disease-induced undershoot
  // will also introduce random death to hit population data
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//deathrate interpolation wrapper
double parameters::rpop( void ){
  // this is an interpolation of the relative population data read in
  double ans(0);
  double yr = time;		// copied into parz during age update
  if ( yr < rpop_yrz.front() ){
    ans = rpop_dat.front();
  } else if ( yr > rpop_yrz.back() ){
    ans = rpop_dat.back();
  } else {
    ans = splinefun( yr, rpop_yrz, rpop_dat, rpop_spl );
  }
  return ans;		
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//migration rate interpolation wrapper
double parameters::mr( double yr ){
  if ( yr < mr_yrz.front() ){
    return mr_dat.front() / 1000;
  } else if ( yr > mr_yrz.back() ){
    return mr_dat.back() / 1000;
  } else {
    return splinefun( yr, mr_yrz, mr_dat, mr_spl ) / 1000;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HIVrate interpolation wrapper
double parameters::h1549( double yr ){
  double ans = 0;
  if( hivincF ){		      // using file-read spline-fit data
    if ( yr < h1549_yrz.front() || yr > h1549_yrz.back() ){
      ans = 0;
    } else {
      ans = splinefun(  yr, h1549_yrz, h1549_dat, h1549_spl ) / 100;
    }
  } else{//this one is via Salamon/Murray - 5 parameters

    double e = 2.71828183;
    double decline = ( yr > 2010 ? exp(-decline2010*(yr-2010)/100 ) : 1.0 );
    yr = yr - h_t0;//peaks at h_beta*(h_alpha-1)
    if( yr > 0 ){
      if( yr <= h_beta * (h_alpha - 1) ){
	ans = pow( e*yr, h_alpha-1 ) * exp( - yr/h_beta ) / pow( h_beta*(h_alpha-1), h_alpha-1 );
      } else {
	ans = h_theta + (1-h_theta) * pow( e*yr, h_alpha-1 ) * exp( - yr/h_beta ) / pow( h_beta*(h_alpha-1), h_alpha-1 );
      }
      ans *= h_peak / 100;
      ans *= decline;
    }
  } // end else
  return ans;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double parameters::fastrisk( double x ){
  double A(0.138), B(0.0406);	// no idea where these numbers come from
  double a(10), b(20);
  double ans(0);

  // real EV version has A and B swapped!
  double tmp(0);
  tmp = A; A = B; B = tmp;	// swap

  if( x < a ){			
    ans = A;
  }
  if ( x >= a && x < b ){
    ans = A + ( B - A ) * ( x - a ) / ( b - a );
  }
  if ( x >= b ){
    ans = B;
  }

  ans /= B;			// new! rr - normalised to 1

  return ans * fscale;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//age-mixing function
double parameters::A( double agee, double agor ){
  //this may be slower than giving people states...
  double ans(1);
  if( AGEflag ){  //failsafe
    unsigned int ae(0), ar(0);//age indices to be looked up
    //nudging, given the end bin age is meaningless
    if( agee > agebins.back() ){agee = agebins.back() - 0.1;}
    if( agor > agebins.back() ){agor = agebins.back() - 0.1;}
    //finding the bin
    ae = bisect_find( agee, agebins );
    ar = bisect_find( agor, agebins );
    //looking up the answer
    ans = agedata.at(ar).at(ae);//NB the order
  }
  return ans;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//what is the ART coverage at a given point in time
double parameters::ARTcov( double t ){
  double y = ( t > ARTst ? (t-ARTst)/ARTrt : 0.0 );
  //  y = pow( y, 2 );
  y = y/(1+y);
  y *= ARTcovmax;
  return y;
}
//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////
////////////////////////////////
//what is the ART coverage at a given point in time
double parameters::extratime( double lifeleft, double cd4 ){
    // double LE = 17.10 - 0.15 * age;//extra time based on mills 
    // // data from Uganda
    // if( cd4 > 100 ){LE += 11.9*2;} else{
    //   if( cd4 > 50 ){ LE += 11.9;}
    // }
  // double Lmax = 60.0;
  // double LE = ( (Lmax-age) > 0 ? (Lmax-age) : 0 );
  double LE = ( lifeleft > 0 ? lifeleft : 0 );
  LE *= pow( cd4 / cd40, expET );
  return LE;
}
//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///////////////// utility inverse
double parameters::invHexp( double a, double b, double x ){
  double ans(0);
  if ( b < 1e-6*exp(-a)*x ){
    ans = exp(-a)*x;
  } else {
    ans = log( 1 + b*exp(-a)*x ) / b;
  }
  return ans;
}

/////////////////////////////////////////////////////////////////////////////
///////////////// utility non-inverse
double parameters::Hexp( double a, double b, double x ){
  double ans = exp(a) * (exp( b * x ) - 1 ) / b;
  return ans;
}



//////////////////////////////////////////////////////////////////////////////////////////////
// lifeexpectancy for age-dependent mortality
double parameters::lifeexpectancy( const gsl_rng *r, double E0, double age ){
  if( E0 > 99 ){cout << "Warning! E0 very high - may cause trouble" << endl;}
  unsigned int top = bisect_find( E0, UNLTb.at(0) ); // the top col for ages
  double p  = (UNLTb.at(0).at(top) - E0) / (UNLTb.at(0).at(top) - UNLTb.at(0).at(top-1) );
  // how close to bottom age ^
  vector< double > surv( UNLTb.size(), 0 );
  vector< double > ages( UNLTb.size(), 0 );
  for( unsigned int ui = 1; ui < UNLTb.size(); ++ui ){
    surv.at(ui) = p * UNLTb.at(ui).at(top-1) + (1-p) * UNLTb.at(ui).at(top);
    surv.at(ui) *=  -1;		// make negative to deal with bisect_find for increasing
    ages.at(ui) = UNLTb.at(ui).at(0);
  }
  surv.at(0) = -(1 + 1e-9); ages.at(0) = 0.0;
  double fact( 1.0 );
  if( age > dt ){		// instead use conditional lifeexpectancy
    top = bisect_find( age, ages ); // which age bin are they in
    p = (ages.at(top) - age) / (ages.at(top) - ages.at(top-1) ); // how close to bottom
    fact = p * surv.at(top-1) + (1-p) * surv.at(top);
    fact *= -1;			// this is S(age); use S(l|a)=S(l)/S(a)
  }
  double u = - gsl_rng_uniform(r) * fact;
  // cout << -u << endl;
  // for( unsigned int ui = 0; ui < surv.size(); ++ui ){
  //   cout << -surv.at(ui) << "," << UNLTb.at(ui).at(0)<<";";
  // }
  // cout << endl;
  double ans = 5.0;
  if( u > surv.back() ){
    ans = 99;
  } else {
    top = bisect_find( u, surv ); // which age cat
    p = (surv.at(top) - u) / (surv.at(top)-surv.at(top-1)); // how close to bottom
    // cout << top << "|" << p << "|"<<UNLTb.at(top).at(0)<<endl;
    //ans = p * UNLTb.at(top-1).at(0) + (1-p) * UNLTb.at(top).at(0); // interpolate ages
    ans = p * ages.at(top-1) + (1-p) * ages.at(top); // interpolate ages
  }

  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
// lifeexpectancy for age-dependent mortality -- overloaded wrapper
// also deals with LE given age...
double parameters::lifeexpectancy( const gsl_rng *r, double age ){
  double e0 = dr( time );
  double ans = lifeexpectancy( r, e0, age );
  // cout << time <<","<< age <<"," << ans << endl;
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::ARTcd4asymptote( double cd4A ){
  double ans = cd42a / (1 + cd42b * exp( -cd42c * cd4A ));
  if( ans > 0.75 * cd40 ){ ans = 0.75 * cd40;}
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::Rxpfn( ){	// done like this because of multi-run worry
  double ans = Rxp;
  if( Rxflag && time > Rxtime ){ ans = Rxp2;}
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::Rxppfn( ){
  double ans = Rxpp;
  if( Rxflag && time > Rxtime ){ ans = Rxpp2;}
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::GL_g( ){
  // the function that determines GL initially (at age 15)
  return gg;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::GL_p( ){
  // the function that determines rates of becoming GL while HIV-
  double ans = ( time > glrst ? (time - glrst) : 0.0 );
  ans = pow( ans / pscale, pn );
  ans /= (1.0 + ans);
  ans *= pmaxval;
  if( time > 2011.5 ){ 
    double frac = GLfracnoteli; // guideline ready fraction amongst no ART
    ans = 5.0 * (frac + .05) * ( SQGLprop - frac);
    if( !SQ && time > 2014 ){ans = 15.0 * (frac + .05) * ( 0.8 - frac);}
  } 
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::GL_newp( double frac ){
  // the function that determines rates of becoming GL while HIV-
  double pm  = pmaxval;
  if(time > 2014){pm = 2*pm;}
  double ans = pscale * (frac + .05) * (pm - frac);
  if( time < glrst ){ans = 0;}
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::GL_q( ){
  // the function that determines the rate of becoming GL while HIV+
  return GL_p()*kq;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::GL_newq( double frac ){
  // the function that determines the rate of becoming GL while HIV+
  double ans = pn * (frac + .05) * (kq - frac);
  if( time < glrst ){ans = 0;}
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double parameters::IPTcovPf( void ){
  double ans = 0;
  if( time > ARTst ){		// NB the ART starttime
    ans = IPTcovP1;
    if( time > IPTcovPt1 ){
      if( time < IPTcovPt2 ){	// in between
	ans = IPTcovP1 * (IPTcovPt2 - time) + IPTcovP2 * (time - IPTcovPt1);
	ans /= (IPTcovPt2 - IPTcovPt1);
      } else {			// ramped up
	ans = IPTcovP2;
      }
    } // started to ramp
  }
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//end parameter member funs
