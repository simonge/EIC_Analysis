import ROOT
from ROOT import TMVA, TFile, TTree, TCut
from subprocess import call
from os.path import isfile
 
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.optimizers import Adam
 
# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
 
ROOT.EnableImplicitMT();
   
output = TFile.Open('tempnew2.root', 'RECREATE')
factory = TMVA.Factory('TMVARegression', output,
        '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=Regression')
 
data = TFile.Open('/scratch/EIC/Analysis/temp.root')
tree = data.Get('temp')
tree.SetEntries(500000)
 
dataloader = TMVA.DataLoader('dataset')

dataloader.AddVariable("real_cut[0].fCoordinates.fY")
dataloader.AddVariable("real_cut[0].fCoordinates.fZ")
dataloader.AddVariable("real_vector[0].fCoordinates.fX")
dataloader.AddVariable("real_vector[0].fCoordinates.fY")
        
dataloader.AddTarget('eE')
dataloader.AddTarget('thetaE')
dataloader.AddTarget('cos(phiV)')
dataloader.AddTarget('sin(phiV)')
 
dataloader.AddRegressionTree(tree, 1.0)
dataloader.PrepareTrainingAndTestTree(TCut('(Tag1_4||Tag2_4)&&iFilter'),
        'nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V')
 
# Generate model
 
# Define model
model = Sequential()
model.add(Dense(64, activation='tanh', input_dim=4))
model.add(Dense(4, activation='linear'))
 
# Set loss and optimizer
model.compile(loss='mean_squared_error', optimizer=Adam())
 
# Store model to file
model.save('modelRegression.h5')
model.summary()
 
# Book methods
factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras',
        'H:!V:VarTransform=D,G:FilenameModel=modelRegression.h5:FilenameTrainedModel=trainedModelRegression.h5:NumEpochs=200:BatchSize=1024')
#factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDTG',
#        '!H:!V:VarTransform=D,G:NTrees=1000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4')
 
# Run TMVA
factory.TrainAllMethods()
#factory.TestAllMethods()
#factory.EvaluateAllMethods()
