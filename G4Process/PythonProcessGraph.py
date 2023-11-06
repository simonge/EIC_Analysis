import uproot
import numpy as np
import awkward as ak
import tensorflow as tf

np.set_printoptions(precision=10,threshold=1000)

inDir  = '/scratch/EIC/G4out/'

for fNb in range(0,1):
    #inName        = 'qr_10x100_ab'+str(fNb)+'.edm4hep.root'
    #inName        = 'qr_18x275_beam_Small4_'+str(fNb)+'.edm4hep.root'
    inName        = 'qr_bx_18x275_Small_'+str(fNb)+'.edm4hep.root'
    #inName        = 'qr_10x100_ab'+str(fNb)+'.edm4hep.root'
    fileNameInput = inDir+inName
    fileInput     = uproot.open(fileNameInput)
    print(fileNameInput)
    print(fileInput.keys())

    events = fileInput['events']

    # Decode LowQ2 id fields
    moduleID     = events["TaggerTrackerHits.cellID"].array()>>8&0b11
    layerID      = events["TaggerTrackerHits.cellID"].array()>>10&0b1111
    xID          = ak.values_astype(events["TaggerTrackerHits.cellID"].array()>>32&0b1111111111111111,np.int16)
    yID          = ak.values_astype(events["TaggerTrackerHits.cellID"].array()>>48&0b1111111111111111,np.int16)
    hitT         = events["TaggerTrackerHits.time"].array()
    hitE         = events["TaggerTrackerHits.EDep"].array()

    # Take truth values
    MCParticlesTruth = events.arrays(["MCParticles.momentum.x","MCParticles.momentum.y","MCParticles.momentum.z","MCParticles.generatorStatus","MCParticles.PDG"])

    # Separate truth values
    lowQ2ref     = events["TaggerTrackerHits#0.index"].array()
    hitparts     = MCParticlesTruth[lowQ2ref]
    xMom         = hitparts["MCParticles.momentum.x"]
    yMom         = hitparts["MCParticles.momentum.y"]
    zMom         = hitparts["MCParticles.momentum.z"]
    PID          = hitparts["MCParticles.PDG"]
    genStat      = hitparts["MCParticles.generatorStatus"]

    # Make awkward array
    zipped = ak.zip({"modID":moduleID,"layID":layerID,"xID":xID,"yID":yID,"hitT":hitT,"hitE":hitE,"xMom":xMom,"yMom":yMom,"zMom":zMom,"pid":PID,"genStat":genStat,"refID":lowQ2ref})
    #zipped.show(limit_rows=(1000), limit_cols=(10000),type=True)

    # Convert awkward to pandas
    df = ak.to_dataframe(zipped)

    # Convert pandas to tensorflow
    rowIDs = df.index.get_level_values(0)
    ragA = tf.RaggedTensor.from_value_rowids(df,rowIDs)

    print(ragA)
    print(ragA.shape)
    
