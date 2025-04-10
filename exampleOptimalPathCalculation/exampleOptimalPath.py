import sys
sys.path.append('/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone6/someToolsForProteinTopologyAnalysisWithSMD')

import objects
import functions

#######################################################################################################
#Parameters
#######################################################################################################

pdbFile = '/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone6/someToolsForProteinTopologyAnalysisWithSMD/exampleFiles/trajData/dryChimera10_autopsf.pdb'

trajectoryFile = '/home/alpal/projects/SMD/gaussianIntegralPaperCode/githubClone6/someToolsForProteinTopologyAnalysisWithSMD/exampleFiles/trajData/WT700ns.chimera.10.RepNumber.5.nowat.dcd'

#######################################################################################################
#Main
#######################################################################################################

initialStructure = objects.frameStructure(pdbFile)
contactDF = initialStructure.getContacts(cutoff=9)


#pearsonTensor = objects.pearsonCorrelationTensor(numpyDir)
#frameEdges = pearsonTensor.getFrameEdges(i)     
#network = objects.residueNetwork(contactDF, frameEdges)
#minResidueIndex, maxResidueIndex = network.getLargestConnectedComponent()
#network.generateOptimalPathsFile(minResdiueIndex, maxResidueIndex, /path/to/saved/optimal/path/file)
