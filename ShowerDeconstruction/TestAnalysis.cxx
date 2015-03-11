#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "Exception.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/internal/BasicRandom.hh"

#include "Pythia8/Pythia.h" //rizki
#include "Pythia8Plugins/HepMC2.h" //rizki

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "AnalysisParameters.h"
#include "Exception.h"

#include "Message.h"

#include "HBBModel.h"
#include "BackgroundModel.h"
#include "ISRModel.h"
#include "Deconstruct.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h" //added by rizki
#include "TH1D.h"
#include "TCanvas.h"

#include "ParseUtils.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace Deconstruction;
using namespace fastjet;
using namespace Pythia8;

int main(int argc, char *argv[]) {

  LOGLEVEL(INFO); // set the level of messages you want to see here //commented by rizki
  // this will be ignored if compiled with -DNOLOG, which would be faster ...
  //LOGLEVEL(DEBUG); //added by rizki

  std::string inputcard = "input_card.dat";
  std::string outFileName = "out_Top:qqbar2ttbar_2btagcondition_dR015.root";
  std::string LHAname = "lha.lhe";
  std::string analysisFile = "analysis_test.cmnd";
  int help = 0;
  static struct extendedOption extOpt[] = {

        {"help",          no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"inputcard",     required_argument,     0, 'c', "The input card text file with one parameter per line in the format 'name  value'.", &inputcard, extendedOption::eOTString},
        {"out",           required_argument,     0, 'o', "Output ROOT file with analysis' results.", &outFileName, extendedOption::eOTString},
        {"lhaOut",        required_argument,     0, 'l', "LHA output filename", &LHAname, extendedOption::eOTString},
        {"analysisFile",  required_argument,     0, 'a', "Pythia 8 analysis input file.", &analysisFile, extendedOption::eOTString},
        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

  if (!parseArguments(argc, argv, extOpt) || help) {
    dumpHelp("TestWAnalysis", extOpt, "TestWAnalysis:\nGenerates Pythia events and runs Shower Deconstruction on them looking for W signal.\n\n");
    return 0;
  } else {
    std::cout << "Command line: ";
    for (int k = 0; k < argc; ++k) std::cout << argv[k] << " ";
    std::cout << endl;
    std::cout << "Input options:" << std::endl;
    dumpOptions(extOpt);
  }

  // Create class to store parameters
  AnalysisParameters param(inputcard);

  // SET UP PYTHIA:
  ///// Initialize HEPMC:
  // write out to HEPMC File
  // Interface for conversion from Pythia8::Event to HepMC one. 
  //HepMC::I_Pythia8 ToHepMC; // commented by rizki
  HepMC::Pythia8ToHepMC ToHepMC; //added by rizki
  ToHepMC.set_crash_on_problem();

  // Specify file where HepMC events will be stored.
  //  HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);
  // Following two lines are deprecated alternative.
  // HepMC::IO_Ascii ascii_io(argv[2], std::ios::out);
  // HepMC::IO_AsciiParticles ascii_io(argv[2], std::ios::out);


  // for SD
  HBBModel *signal = 0;
  BackgroundModel *background = 0;
  ISRModel *isr = 0;
  Deconstruct *deconstruct = 0;

  signal = new HBBModel(param);
  background = new BackgroundModel(param);
  isr = new ISRModel(param);
  deconstruct = new Deconstruct(param, *signal, *background, *isr);

  // for ROOT analysis
  TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
  TTree *tree = new TTree("analysis", "");

  double h_time_sd = 0;
  double h_sig = 0;
  double h_bkg = 0;
  double h_sb = 0;

  tree->Branch("h_time_sd", &h_time_sd);
  tree->Branch("h_sig", &h_sig);
  tree->Branch("h_bkg", &h_bkg);
  tree->Branch("h_sb", &h_sb);

  double xsec = 0;

  // Generator.
  Pythia pythia;

  // Create an LHAup object that can access relevant information in pythia.
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  myLHA.openLHEF(LHAname.c_str());
  myLHA.setInit();
  myLHA.initLHEF();

  // Read in commands from external file.
  pythia.readFile(analysisFile);

  int nEvents = pythia.mode("Main:numberOfEvents");

  // initialize randomnumber for b-tag:
  srand(pythia.mode("Random:seed"));

  // Initialization. Beam parameters set in .cmnd file.
  pythia.init();

  for (int i = 0; i < nEvents; ++i) {
    //if (i%1000 == 0) //commented by rizki
    //if (i<4300) continue; //added by rizki
    if (i%1 == 0) //added by rizki
      std::cout << "Event " << i << "/" << nEvents << endl;
    if (!pythia.next())
      continue;


    // Fill also LHE-Output File
    myLHA.setEvent();
    myLHA.eventLHEF();

    vector<PseudoJet> hadrons;

    for (int p = 0; p < pythia.event.size(); ++p) {
      int pid = (int) fabs(pythia.event[p].id());

      // particles from final state
      if (pythia.event[p].isFinal()) {
       PseudoJet mom(pythia.event[p].px(),
                     pythia.event[p].py(),
                     pythia.event[p].pz(),
                     pythia.event[p].e());

        // neutrinos
        if (pid == 12 || pid == 14|| pid == 16) {
          continue;
        }
        hadrons.push_back(mom);

      } // end if final-states
    }
    hadrons = sorted_by_pt(hadrons);

    // hadrons are in the hadrons vector
    for(unsigned int ii = 0; ii < hadrons.size(); ii++) {
      if(std::fabs(hadrons[ii].rap()) > 5.0) {
        hadrons.erase(hadrons.begin()+ii);
        ii--;
      }
    }

    // build jets
    //JetDefinition jet_def(fastjet::cambridge_algorithm, 1.0); // commented by rizki
    JetDefinition jet_def(fastjet::cambridge_algorithm, param[p_R]); // added by rizki

    ClusterSequence clust_seq(hadrons, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(350));

    //std::cout << "number of Fat Jets = " << jets.size() << endl;
    LOG(DEBUG) << " jets: " << endl; // commented by rizki
    for(unsigned ii=0; ii<jets.size(); ii++) {
      LOG(DEBUG) << "user_index = " << jets[ii].user_index()<<endl; // commented by rizki
      printjet(jets[ii]);
    }
    LOG(DEBUG) << endl; // commented by rizki

    // large-R jets are done
    // get constituents
    vector<vector<PseudoJet> > constjets;
    for(unsigned k = 0; k < jets.size(); k++) {
      constjets.push_back(clust_seq.constituents(jets[k]));
    }

    LOG(DEBUG) << "fat jets: " << endl; // commented by rizki
    for(unsigned jj=0; jj<jets.size(); jj++) {
      printjet(jets[jj]);
    }

    if (jets.size() < 1)
      continue;
    //std::cout << "leading Fat Jet pT = "<< jets[0].pt() << endl; //rizki   
    bool passM = true;
    bool passEta = true;
    for (int k = 0; k < 2; ++k) {
      if (jets[k].m() < 100) {
        passM = false;
      }
      if(std::fabs(jets[k].rap()) > 2.0) {
        passEta = false;
      }
    }
    if (!passM) continue;
    if (!passEta) continue;

    vector<PseudoJet> constits1 = clust_seq.constituents(jets[0]);

    LOG(DEBUG) << "constiuents of fat jet1: " << constits1.size() << endl; // commented by rizki
    //std::cout << "constiuents of fat jet1: " << constits1.size() << endl; // added by rizki

    // collect constits in small cones:
    // start value for small jets

    //JetDefinition jet_def_small(fastjet::kt_algorithm, 0.2); // commented by rizki
    JetDefinition jet_def_small(fastjet::kt_algorithm, param[p_Rsmall]); // added by rizki
    ClusterSequence clust_seq_small1(constits1, jet_def_small);
    vector<PseudoJet> jets_small1 = sorted_by_pt(clust_seq_small1.inclusive_jets(15));

    LOG(DEBUG) << "first fat jet,  microjetconesize: " << param[p_Rsmall] << ", has " << jets_small1.size() << " microjet(s). " << endl; // commented by rizki
    //std::cout << "first fat jet,  microjetconesize: " << param[p_Rsmall] << ", has " << jets_small1.size() << " microjet(s). " << endl; // addded by rizki

    LOG(DEBUG) << "jets_small before erase: " << jets_small1.size() << endl; // commented by rizki
    for (unsigned int ii = 0; ii < jets_small1.size(); ii++) {
      LOG(DEBUG) << "user_index = " << jets_small1[ii].user_index() << std::endl; // commented by rizki
      printjet(jets_small1[ii]);
    }

    if (jets_small1.size() > 9) {
      jets_small1.erase(jets_small1.begin() + (int) 9,
                        jets_small1.begin() + jets_small1.size());
    }
    LOG(DEBUG) << "jets_small after erase: " << jets_small1.size() << endl; // commented by rizki
    for (unsigned ii=0; ii< jets_small1.size(); ii++) {
      printjet(jets_small1[ii]);
      //std::cout << "microjet no." << ii << ": pt = "<<jets_small1[ii].pt()<< ", user_inder = "<< jets_small1[ii].user_index()<<endl; //added by rizki
    }

    //testing btagging rizki - start
    for (unsigned mji=0; mji< jets_small1.size(); mji++) { //loop microjet
      //std::cout <<"mji = "<< mji <<std::endl;
      double mjeta = jets_small1[mji].eta();
      double mjphi = jets_small1[mji].phi();
	for (int p = 0; p < pythia.event.size(); ++p) { //loop pythia hadrons
	  //if(std::fabs(mjeta-pythia.event[p].eta())<1e-2 && std::fabs(mjphi-pythia.event[p].phi()<1e-2)){
	  double peta = pythia.event[p].eta();
	  double pphi = pythia.event[p].phi();

	  double dEta = mjeta - peta;
	  double dPhi = mjphi - pphi;
	  while (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
	  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
	  double dR = sqrt(dEta*dEta + dPhi*dPhi);

	  if(std::fabs(dR < param[p_Rsmall])){
	    int pid = (int) fabs(pythia.event[p].id());
	    //std::cout<<"we have a microjet-hadron match! pdgid = "<< pid<< endl;
	    std::div_t divresult;
	    divresult = div(pid, 500);
	    if( pid ==5 ||  divresult.quot==1){
	      //std::cout << "Btagged! -----> b quark/B hadron : pdgid = " <<pid << std::endl;
	      jets_small1[mji].set_user_index(1);
	      break;
	    }
	  }
	}
    }
    
    int btagi = 0;
    for (unsigned mji=0; mji< jets_small1.size(); mji++) { //check how many btagged microjets
      if(jets_small1[mji].user_index()==1) btagi++;
      //only let events with 2 btags among two high pt mirojets
      //if(mji==1 && btagi<2) break;
    }

    if(btagi<2){
      std::cout <<"btag microjets less than 2, discarding event "<< i << std::endl;
      continue; //if less than two discard event
    }
    //testing btagging rizki - end 

    vector<PseudoJet> input_deconstruction1;
    for(unsigned int k = 0; k < jets_small1.size(); ++k) {
      PseudoJet momsd(jets_small1[k].px(),
                      jets_small1[k].py(),
                      jets_small1[k].pz(),
                      jets_small1[k].e()+0.1);
      // we shift the energy of the microjets such that it is not exactly on-shell (otherwise the code breaks) this value should be checked under experimental considerations...
      momsd.set_user_index(jets_small1[k].user_index());
      input_deconstruction1.push_back(momsd);
    }

    input_deconstruction1 = sorted_by_pt(input_deconstruction1);

    // SD or ED follows
    // calculate the number of events used:
    xsec += 1;

    // first do the Shower Deconstruction
    double wSignal;
    double wBackground;
    double chi;
    wSignal = 0;
    wBackground = 0;
    chi = 0;

    std::multimap<double, std::vector<fastjet::PseudoJet> > signalWeight;

    clock_t clo_start_sd = clock();

    try {
      chi = deconstruct->deconstruct(input_deconstruction1, wSignal, wBackground);
      signalWeight = deconstruct->signalWeight();
    } catch(Deconstruction::Exception &e) {
      std::cout << "Exception while running SD: " << e.what() << std::endl;
    }

    clock_t clo_end_sd = clock();

    //fastjet::PseudoJet top1(0,0,0,0);
    //top1 = sum(signalWeight.rbegin()->second);

    h_bkg = (wBackground);
    //std::cout<<"wBackground = "<< wBackground << std::endl; //added by rizki
    h_sig = (wSignal);
    //std::cout<<"wSignal = "<< wSignal << std::endl; //added by rizki
    h_sb = (chi);
    cout << scientific;
    cout.precision(10);
    std::cout << "Chi = " << chi << std::endl;
    std::cout << "" << std::endl; //added by rizki

    h_time_sd = (clo_end_sd - clo_start_sd);

    tree->Fill();
  }
  std::cout << "Number of events that passed selection = " << xsec << " of a total of " << nEvents << std::endl;
  std::cout << "Generated cross section = " << pythia.info.sigmaGen() << ", err = " << pythia.info.sigmaErr() << " [mb]" << std::endl;
  std::cout << "Luminosity = " << nEvents/pythia.info.sigmaGen() << " [mb^-1]" << std::endl;
  xsec = xsec/(nEvents/(pythia.info.sigmaGen()));
  std::cout << "Fiducial cross section = " << xsec << "[mb]" << std::endl;

  //  pythia.statistics(); //commented by rizki
  pythia.stat();
  myLHA.updateSigma();
  myLHA.closeLHEF(true);

  outFile->Write();
  outFile->Close();

  return 0;
}

