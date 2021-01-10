#include <TSystem.h>

R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libanamodule)

int run(const int nEvents = 1)
{
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;
  const bool cosmic = true;

  recoConsts* rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  if(cosmic)
  {
    rc->init("cosmic");
    rc->set_BoolFlag("COARSE_MODE", true);
    rc->set_DoubleFlag("KMAGSTR", 0.);
    rc->set_DoubleFlag("FMAGSTR", 0.);
  }
  rc->Print();

  GeomSvc::UseDbSvc(true);  
  GeomSvc* geom_svc = GeomSvc::instance();

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(0);

  AnaModule* ana = new AnaModule();
  ana->set_output_filename("ana.root");



  ana->registerDetector("D2X");     //register detector by its name, all detectors that do not directly partipate the tracking can be used
  ana->registerDetector("D2Xp");
  ana->registerDetector("D3pX");
  ana->registerDetector("D3pXp");     
  ana->registerDetector("D3mX");     
  ana->registerDetector("D3mXp");
     
  ana->registerDetector("DP1TR");
  ana->registerDetector("DP1BL");
  ana->registerDetector("DP1BR");
  ana->registerDetector("DP2TL");
  ana->registerDetector("DP2TR");
  ana->registerDetector("DP2BL");

  ana->registerDetector("H1B");
  ana->registerDetector("H1T");
  ana->registerDetector("H2B");
  ana->registerDetector("H2T");
  ana->registerDetector("H2L");
  ana->registerDetector("H2R");
  ana->registerDetector("H3B");
  ana->registerDetector("H3T");
  ana->registerDetector("H4B");
  ana->registerDetector("H4T");
/*
  ana->registerDetector("H4Tu");
  ana->registerDetector("H4Td");
  ana->registerDetector("H4Bu");
  ana->registerDetector("H4Bd");
*/
 ana->registerDetector("H4Y1L");
  ana->registerDetector("H4Y1R");
  ana->registerDetector("H4Y2L");
  ana->registerDetector("H4Y2R");

  /*  ana->registerDetector("P1X1");
  ana->registerDetector("P1X2");
  ana->registerDetector("P2X1");
  ana->registerDetector("P2X2");
  ana->registerDetector("P1Y1");
  ana->registerDetector("P1Y2");
  ana->registerDetector("P2Y1");
  ana->registerDetector("P2Y2");
*/  
se->registerSubsystem(ana);

  Fun4AllInputManager* in = new Fun4AllDstInputManager("DSTIN");
  in->Verbosity(0);
 //in->fileopen("result.root");
 se->registerInputManager(in);

  //we do not really need an output manager
  //Fun4AllDstOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", "result.root");
  //se->registerOutputManager(out);

//multpleDST run
const char* fn_list_dst="list_dst.txt";
const int   n_dst_ana=0;
vector<string> list_dst;
  string fn_dst;
  ifstream ifs(fn_list_dst);
  while (ifs >> fn_dst) list_dst.push_back(fn_dst);
  ifs.close();

  int n_dst = list_dst.size();
  cout << "N of DSTs: all = " << n_dst;
  if (n_dst_ana > 0 && n_dst > n_dst_ana) n_dst = n_dst_ana;
  cout << ", to be analyzed = " << n_dst << endl;
  for (int i_dst = 0; i_dst < n_dst; i_dst++) {
    string fn_dst = list_dst[i_dst];
    cout << "DST: " << i_dst << "/" << n_dst << ": " << fn_dst << endl;
    in->fileopen(fn_dst);
    se->run();
    //se->run(nEvents);
  }
//se->run(nEvents);


  //se->run(nEvents);

  // finish job - close and save output files
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;

  delete se;
  gSystem->Exit(0);
  return 0;
}
