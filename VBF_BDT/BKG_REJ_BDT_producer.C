#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


string bdt_type = "VBF";

int BKG_REJ_BDT_producer(TString rootIn=bdt_type + "_BKG_REJ_nopTt",TString myMethodList = "" )
{
  TMVA::Tools::Instance();
  std::map<std::string,int> Use;
  // Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  // ---------------------------------------------------------------
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;


  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return 1;
      }
      Use[regMethod] = 1;
    }
  }
  
  // --------------------------------------------------------------------------------------------------
  TString fname_sig;
  if(bdt_type == "ggF"){
   fname_sig = "ntuples/photon_bdt_ntuples_sig_ggF_ALL.root";
  } else { 
   fname_sig = "ntuples/photon_bdt_ntuples_sig_VBF_ALL.root";
  }
  TString fname_bkg = "ntuples/photon_bdt_ntuples_bkg_ALL.root"; //photon_bdt_ntuples_ALL_nott.root";//

  if (gSystem->AccessPathName( fname_bkg ) || gSystem->AccessPathName( fname_sig ))  // file does not exist in local directory
    gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");


  TFile *input_bkg = TFile::Open( fname_bkg );
  TFile *input_sig = TFile::Open( fname_sig );
  std::cout << "--- TMVAClassification       : Using sig input file: " << input_sig->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using bkg input file: " << input_bkg->GetName() << std::endl;

  // Register the training and test trees
  TTree *signalTree     = (TTree*)input_sig->Get("tree");
  TTree *backgroundTree = (TTree*)input_bkg->Get("tree");


  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "MVA_output/TMVA_" + rootIn + ".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( ("TMVAClassification_" + rootIn) , outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

 if(bdt_type=="ggF"){
  dataloader->AddVariable("min_dR", 'F');
//  dataloader->AddVariable("pTt", 'F');
  dataloader->AddVariable("pt_mass", 'F');
  dataloader->AddVariable("max_dR", 'F');
  dataloader->AddVariable("cosTheta", 'F');
  dataloader->AddVariable("costheta", 'F');
  dataloader->AddVariable("photon_res", 'F');
  dataloader->AddVariable("l1_rapidity", 'F');
//  dataloader->AddVariable("phi", 'F');
  dataloader->AddVariable("photon_mva", 'F');
//  dataloader->AddVariable("photon_prap", 'F');
  dataloader->AddVariable("l2_rapidity", 'F');
 } else{
  dataloader->AddVariable("cosTheta", 'F');
//  dataloader->AddVariable("pTt", 'F');
  dataloader->AddVariable("photon_prap", 'F');
  dataloader->AddVariable("max_dR", 'F');
  dataloader->AddVariable("l2_rapidity", 'F');
  dataloader->AddVariable("min_dR", 'F');
  dataloader->AddVariable("pt_mass", 'F');
  dataloader->AddVariable("costheta", 'F');
  dataloader->AddVariable("photon_mva", 'F');
 }

/*
  //Variables to choose from
  dataloader->AddVariable("min_dR", 'F');
  dataloader->AddVariable("photon_mva", 'F');
  dataloader->AddVariable("max_dR", 'F');
  dataloader->AddVariable("pt_mass", 'F');
  dataloader->AddVariable("cosTheta", 'F');
  dataloader->AddVariable("costheta", 'F');
  dataloader->AddVariable("photon_res", 'F');
  dataloader->AddVariable("photon_prap", 'F');
  dataloader->AddVariable("l1_rapidity", 'F');
  dataloader->AddVariable("l2_rapidity", 'F');
  dataloader->AddVariable("phi", 'F');
  dataloader->AddVariable("pTt", 'F');
*/



  //Applying the weight to the MC samples
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  // Apply additional cuts on the signal and background samples (can be different)
  // Set up an easy way to tell the testing and training sample apart
  TCut niceCut;
  if(bdt_type == "ggF"){
    niceCut= "binning_bdt < 0.7 && drmin_yj > 0.4";
  } else {
    niceCut= "binning_bdt > 0.7 && drmin_yj > 0.4";
  }
  TCut training = niceCut + "event_number%2==0";
  TCut testing  = niceCut + "event_number%2==1";


  dataloader->AddTree( signalTree, "Signal", 1.0, training, TMVA::Types::kTraining);
  dataloader->AddTree( signalTree, "Signal", 1.0, testing,  TMVA::Types::kTesting);

  dataloader->AddTree( backgroundTree, "Background", 1.0, training, TMVA::Types::kTraining);
  dataloader->AddTree( backgroundTree, "Background", 1.0, testing,  TMVA::Types::kTesting);

  dataloader->SetSignalWeightExpression( "weight" );
  dataloader->SetBackgroundWeightExpression( "weight" );



  //This block has all the settings for the BDT.
  //nTrain_Signal is number of signal training events
  //nTrain_Background is number of bkg. training events
  //SplitSeed is seed that randomly splits events (from SplitMode==Random) 
  //dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
  //          "SplitMode=Random:SplitSeed=1000:NormMode=NumEvents:!V" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
      "!H:!V:NTrees=100:MinNodeSize=7.0%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
//      "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=200:MinNodeSize=3.0%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//			 "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  if (Use["BDTB"]) // Bagging
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
			 "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
			 "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
  if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
			 "!H:!V:NTrees=100:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------
  // Save the output
  outputFile->Close();
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  delete factory;
  delete dataloader;
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
  return 0;
}
int main( int argc, char** argv )
{
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod;
  }
  return BKG_REJ_BDT_producer(methodList);
}
