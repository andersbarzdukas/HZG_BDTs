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


int VBF_BDT_producer(TString rootIn="VBF_test_noypt",TString myMethodList = "" )
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
  TString fname_sig = "ntuples/ntuples_sig_VBF_all.root";
  TString fname_bkg = "ntuples/ntuples_sig_ggF_all.root";

  if (gSystem->AccessPathName( fname_bkg ) || gSystem->AccessPathName( fname_sig ))  // file does not exist in local directory
    gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");


  TFile *input_bkg = TFile::Open( fname_bkg );
  TFile *input_sig = TFile::Open( fname_sig );
  std::cout << "--- TMVAClassification       : Using sig input file: " << input_sig->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using bkg input file: " << input_bkg->GetName() << std::endl;

  // Register the training and test trees
  TTree *signalTree     = (TTree*)input_sig->Get("tree");
  TTree *background     = (TTree*)input_bkg->Get("tree");


  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "MVA_output/TMVA_" + rootIn + ".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( ("TMVAClassification_" + rootIn) , outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  //Standard variables used in Run 2 Dijet BDT
  dataloader->AddVariable("jj_deta", 'F');
  dataloader->AddVariable("jj_dphi", 'F');
  dataloader->AddVariable("Zyjj_dr", 'F');
  dataloader->AddVariable("pt_bal", 'F');
  dataloader->AddVariable("zep_var", 'F');
  //dataloader->AddVariable("y_pt", 'F');
  dataloader->AddVariable("mjj", 'F');
  dataloader->AddVariable("drmin_yj", 'F');
  dataloader->AddVariable("j1_pt", 'F');
  dataloader->AddVariable("j2_pt", 'F');

  //dataloader->AddVariable("pTt", 'F');

  //Potential variables (mjj is of particular interest)
   //dataloader->AddVariable("kinMVA", 'F');
  //dataloader->AddVariable("Zyjj_dphi", 'F');
//  dataloader->AddVariable("cent_var", 'F');
  //dataloader->AddVariable("drmax_yj", 'F');

//  dataloader->AddVariable("pTt_perp", 'F');
//  dataloader->AddVariable("pTt_wrong", 'F');
//  dataloader->AddVariable("mjj", 'F');


  //Applying the weight to the MC samples
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  dataloader->AddSignalTree    ( signalTree,     signalWeight );
  dataloader->AddBackgroundTree( background, backgroundWeight );

  dataloader->SetSignalWeightExpression( "weight" );
  dataloader->SetBackgroundWeightExpression( "weight" );

  // Apply additional cuts on the signal and background samples (can be different)
  // These cuts should be fairly self explanatory
  TCut niceCut = "dr_jj > 0.4 && drmin_yj > 0.4";
  TCut mycuts = niceCut;
  TCut mycutb = niceCut;


  //This block has all the settings for the BDT.
  //nTrain_Signal is number of signal training events
  //nTrain_Background is number of bkg. training events
  //SplitSeed is seed that randomly splits events (from SplitMode==Random) 
  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
            "SplitMode=Random:SplitSeed=1000:NormMode=NumEvents:!V" );

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
  return VBF_BDT_producer(methodList);
}
