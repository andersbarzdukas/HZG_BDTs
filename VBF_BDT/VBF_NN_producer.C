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
#include "TMVA/TMVAMultiClassGui.h"
     
using namespace TMVA;

int VBF_NN_producer(TString rootIn="Dijet_DNN_3sig_3bkg_withllyvars",TString myMethodList = "" )
//int VBF_NN_producer(TString rootIn="Test",TString myMethodList = "" )
{
  TMVA::Tools::Instance();
  std::map<std::string,int> Use;
  // Boosted Decision Trees
  Use["NN"]            = 0; // uses Neural Net Approach
  Use["BDTG"]          = 1; // uses Gradient Boost
  Use["DL_CPU"]        = 0;
  Use["DNN"]           = 0;

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
  TString fname_sig_ggF = "ntuples/ntuples_sig_ggF_all.root";
  TString fname_sig_VBF = "ntuples/ntuples_sig_VBF_all.root";
  TString fname_sig_VH  = "ntuples/ntuples_sig_VH_all.root";

  TString fname_bkg_DY = "ntuples/ntuples_bkg_DY_all.root";
  TString fname_bkg_tt = "ntuples/ntuples_bkg_tt_all.root";
  TString fname_bkg_EWK = "ntuples/ntuples_bkg_EWKZG_all.root";


  if (gSystem->AccessPathName( fname_bkg_DY ) || gSystem->AccessPathName( fname_bkg_tt ) ||
      gSystem->AccessPathName( fname_sig_ggF ) || gSystem->AccessPathName( fname_sig_VBF ) || gSystem->AccessPathName( fname_sig_VH )){  // file does not exist in local directory
    std::cout << "File is not present" << endl;
    gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");
  }

  TFile *input_sig_ggF = TFile::Open( fname_sig_ggF );
  TFile *input_sig_VBF = TFile::Open( fname_sig_VBF );
  TFile *input_sig_VH  = TFile::Open( fname_sig_VH );

  TFile *input_bkg_DY = TFile::Open( fname_bkg_DY );
  TFile *input_bkg_tt = TFile::Open( fname_bkg_tt );
  TFile *input_bkg_EWK = TFile::Open( fname_bkg_EWK );


  std::cout << "--- TMVAClassification       : Using ggF sig input file: " << input_sig_ggF->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using VBF sig input file: " << input_sig_VBF->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using VH sig input file: "  << input_sig_VH->GetName() << std::endl;

  std::cout << "--- TMVAClassification       : Using DY bkg input file: "  << input_bkg_DY->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using tt bkg input file: "  << input_bkg_tt->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using EWK bkg input file: "  << input_bkg_EWK->GetName() << std::endl;


  // Register the training and test trees
  TTree *signalTree_ggF = (TTree*)input_sig_ggF->Get("tree");
  TTree *signalTree_VBF = (TTree*)input_sig_VBF->Get("tree");
  TTree *signalTree_VH  = (TTree*)input_sig_VH->Get("tree");

  TTree *background_DY  = (TTree*)input_bkg_DY->Get("tree");
  TTree *background_tt  = (TTree*)input_bkg_tt->Get("tree");
  TTree *background_EWK = (TTree*)input_bkg_EWK->Get("tree");


  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "MVA_output/TMVA_NN_" + rootIn + ".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( ("MultiClass_NN_" + rootIn) , outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  //Standard variables used in Run 2 Dijet NN

  dataloader->AddVariable("photon_mva", 'F');
  dataloader->AddVariable("min_dR", 'F');
  dataloader->AddVariable("max_dR", 'F');
  dataloader->AddVariable("pt_mass", 'F');
  dataloader->AddVariable("cosTheta", 'F');
  dataloader->AddVariable("costheta", 'F');
  dataloader->AddVariable("photon_res", 'F');
  dataloader->AddVariable("photon_prap", 'F');
  dataloader->AddVariable("l1_rapidity", 'F');
  dataloader->AddVariable("l2_rapidity", 'F');
//  dataloader->AddVariable("phi", 'F');


  dataloader->AddVariable("jj_deta", 'F');
  dataloader->AddVariable("jj_dphi", 'F');
  dataloader->AddVariable("Zyjj_dr", 'F');
  dataloader->AddVariable("pt_bal", 'F');
  dataloader->AddVariable("zep_var", 'F');
  dataloader->AddVariable("mjj", 'F');
  dataloader->AddVariable("drmin_yj", 'F');
  dataloader->AddVariable("j1_pt", 'F');
  dataloader->AddVariable("j2_pt", 'F');



  
  //Applying the weight to the MC samples
  //Double_t signalWeight     = 1.0;
  //Double_t backgroundWeight = 1.0;

  // Apply additional cuts on the signal and background samples (can be different)
  // These cuts should be fairly self explanatory
  TCut niceCut = "drmin_yj > 0.4";//dr_jj > 0.4 && 
  //TCut mycuts = niceCut;
  //TCut mycutb = niceCut;

  dataloader->AddTree( signalTree_ggF, "Signal_ggF");//, signalWeight  );
  dataloader->AddTree( signalTree_VBF, "Signal_VBF");//, signalWeight );
  dataloader->AddTree( signalTree_VH,  "Signal_VH" );//, signalWeight );

  dataloader->AddTree( background_DY,  "Background_DY" );//, backgroundWeight );
  dataloader->AddTree( background_tt,  "Background_tt" );//, backgroundWeight );
  dataloader->AddTree( background_EWK, "Background_EWKZG" );//, backgroundWeight );

//  dataloader->SetSignalWeightExpression( "weight" );
//  dataloader->SetBackgroundWeightExpression( "weight" );
//  dataloader -> SetWeightExpression("weight");


  //This block has all the settings for the BDT.
  //nTrain_Signal is number of signal training events
  //nTrain_Background is number of bkg. training events
  //SplitSeed is seed that randomly splits events (from SplitMode==Random) 
  dataloader->PrepareTrainingAndTestTree( niceCut, "SplitMode=Random:NormMode=NumEvents:!V" );


  // Boosted Decision Trees
  if (Use["MLP"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
  if (Use["BDTG"]) // gradient boosted decision trees
      factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
  if (Use["NN"]) // Neural Network 
    factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP","!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
  if (Use["DL_CPU"]) {
    TString layoutString("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
    TString trainingStrategyString("TrainingStrategy=Optimizer=ADAM,LearningRate=1e-3,"
    "TestRepetitions=1,ConvergenceSteps=10,BatchSize=100,MaxEpochs=20");
    TString nnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
    "WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
    nnOptions.Append(":");
    nnOptions.Append(layoutString);
    nnOptions.Append(":");
    nnOptions.Append(trainingStrategyString);
    factory->BookMethod(dataloader,TMVA::Types::kDNN,"DNN CPU",nnOptions); //kDL vs DNN?
  }
  if(Use["DNN"]){
      TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");
      // Training strategies.
      TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
      TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1 + "|" + training2;
      // General Options.
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
      // Standard implementation, no dependencies.
      // Multi-core CPU implementation.
      TString cpuOptions = dnnOptions + ":Architecture=CPU";
      factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN CPU", cpuOptions);
   }
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
//  if (!gROOT->IsBatch()) TMVA::TMVAMultiClassGui( outfileName );
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
  return VBF_NN_producer(methodList);
}
