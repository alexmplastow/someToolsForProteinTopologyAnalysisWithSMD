import sys
sys.path.append('/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone8/someToolsForProteinTopologyAnalysisWithSMD')

import objects
import functions

import MDAnalysis as mda
import glob
import tqdm

import matplotlib.pyplot as plt
import seaborn as sns


#TODO: finish adding the new directory for storing PDB files

#######################################################################################################
#Parameters
#######################################################################################################

pdbFile = '/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone8/someToolsForProteinTopologyAnalysisWithSMD/exampleFiles/trajData/dryChimera10_autopsf.pdb'

trajectoryFile = '/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone8/someToolsForProteinTopologyAnalysisWithSMD/exampleFiles/trajData/WT700ns.chimera.10.RepNumber.5.nowat.dcd'

chainName = 'A'

pdbRegen = False

pearsonTensorRegen = True

#NOTE: I would not recomend setting this to true unless you are prepared to eat up the RAM on your
    #NOTE: machine
visualizePearsonTensor = False

#NOTE: you should make the analytical output resemble the old VMD approach closely
    #NOTE: you will need to distinguish between short paths and long paths with low correlations
vmdOptPathFileName = 'VMDoptPath.txt'

contactDataFileName = 'contact.dat'

#######################################################################################################
#Main
#######################################################################################################

#static parameters
startTime = 5
frameSpeed = 0.2

initialStructure = objects.frameStructure(pdbFile)
contactDF = initialStructure.getContacts(cutoff=9)

if pdbRegen:
        #Collecting the frame to be used in the analysis
        u = mda.Universe(pdbFile, trajectoryFile)
        u.frameMax = len(u.trajectory)

        for t in range(0, u.frameMax):

            u.trajectory[t]

            with mda.Writer(f"../tmp/frame_{t:04d}.pdb", u.atoms.n_atoms) as writer:
                writer.write(u)

if pearsonTensorRegen:
    sim = objects.simulationAsPDBs(f"../tmp", chainName = chainName)
    sim.savePearsonTensor(startTime,f"../tmp")

if visualizePearsonTensor:

    pearsonCorrelationTensor = objects.pearsonCorrelationTensor("../tmp")
    fileIndices = range(0, len(glob.glob(f"../tmp/*npy")))

    for i in tqdm.tqdm(fileIndices[startTime:]):
        correlationMatrix = pearsonCorrelationTensor.getMovingAverage(i, 5)

        plt.style.use('dark_background')
        plt.figure(figsize=(10, 8))
        sns.heatmap(correlationMatrix, vmin=-1, vmax=1, cmap="nipy_spectral")
        plt.title(f'Protealysin Dynamical Cross Correlation Heatmap, t={round(i*frameSpeed, 3)} ns')
        plt.xlabel('Residue Index')
        plt.ylabel('Residue Index')
        plt.savefig(f"../tmp/DNAheatMap{i}.png")

pearsonCorrelationTensor = objects.pearsonCorrelationTensor("../tmp")
fileIndices = range(0, len(glob.glob(f"../tmp/*npy")))

for i in tqdm.tqdm(fileIndices[startTime:]):
    frameEdges = pearsonCorrelationTensor.getFrameEdges(i)   
    network = objects.residueNetwork(contactDF, frameEdges)
    network.generateContactDataFile(contactDataFileName)
    #TODO: the graph is completely connected, we should have many suboptimal paths
    network.generateVMDoptimalPathsFile(vmdOptPathFileName) 
