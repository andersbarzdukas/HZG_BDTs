#from ROOT import TMVA, TFile, TTree, TCut, gROOT
import ROOT
from os.path import isfile

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense,Activation
from tensorflow.keras.optimizers import SGD

#SetupTMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output=TFile.Open('TMVA_Dijet_NN_Multiclassifier.root','RECREATE')
factory=TMVA.Factory('TMVAClassification',output,
  '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=multiclass')

#Loaddata
if not isfile('tmva_example_multiple_background.root'):
  createDataMacro=str(gROOT.GetTutorialDir())+'/tmva/createData.C'
  print(createDataMacro)
  gROOT.ProcessLine('.L{}'.format(createDataMacro))
  gROOT.ProcessLine('create_MultipleBackground(4000)')


ggF_signal=TFile.Open('/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/ntuples/ntuples_sig_ggF_all.root')
VBF_signal=TFile.Open('/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/ntuples/ntuples_sig_VBF_all.root')
VH_signal=TFile.Open('/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/ntuples/ntuples_sig_VH_all.root')
tt_background=TFile.Open('/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/ntuples/ntuples_bkg_tt_all.root')
DY_background=TFile.Open('/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/ntuples/ntuples_bkg_DY_all.root')

ggF_tree=ggF_signal.Get('tree')
VBF_tree=VBF_signal.Get('tree')
VH_tree=VH_signal.Get('tree')
tt_tree=tt_background.Get('tree')
DY_tree=DY_background.Get('tree')

dataloader=TMVA.DataLoader('dataset')
dataloader.AddVariable("jj_deta",'F');
dataloader.AddVariable("jj_dphi",'F');
dataloader.AddVariable("Zyjj_dr",'F');
dataloader.AddVariable("pt_bal",'F');
dataloader.AddVariable("zep_var",'F');
dataloader.AddVariable("mjj",'F');
dataloader.AddVariable("drmin_yj",'F');
dataloader.AddVariable("j1_pt",'F');
dataloader.AddVariable("j2_pt",'F');


dataloader.AddTree(ggF_tree,'Signal_ggF')
dataloader.AddTree(VBF_tree,'Signal_VBF')
dataloader.AddTree(VH_tree,'Signal_VH')
dataloader.AddTree(tt_tree,'Background_tt')
dataloader.AddTree(DY_tree,'Background_DY')

dataloader.PrepareTrainingAndTestTree(TCut(''),
  'SplitMode=Random:NormMode=NumEvents:!V')

#Generatemodel

#Definemodel
model=Sequential()
model.add(Dense(32,activation='relu',input_dim=4))
model.add(Dense(5,activation='softmax'))

#Setlossandoptimizer
model.compile(loss='categorical_crossentropy',optimizer=SGD(learning_rate=0.01),weighted_metrics=['accuracy',])

#Storemodeltofile
model.save('modelMultiClass.h5')
model.summary()

#Bookmethods
factory.BookMethod(dataloader,TMVA.Types.kFisher,'Fisher',
  '!H:!V:Fisher:VarTransform=D,G')
factory.BookMethod(dataloader,TMVA.Types.kPyKeras,'PyKeras',
  'H:!V:VarTransform=D,G:FilenameModel=modelMultiClass.h5:FilenameTrainedModel=trainedModelMultiClass.h5:NumEpochs=20:BatchSize=32')

#RunTMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
