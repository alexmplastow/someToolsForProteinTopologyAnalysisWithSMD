import matplotlib.pyplot as plt
import Bio.PDB.PDBParser as parser
import warnings
from Bio import BiopythonWarning
import numpy as np
import pandas as pd
import statistics
from scipy.signal import find_peaks
from scipy.stats import linregress
import glob
from tqdm import tqdm
import functions
import pickle
import math
import os


#I found my code became slightly more concise this way
vmdColors = {
    0: 'blue',
    1: 'red',
    2: 'gray',
    3: 'orange',
    4: 'yellow',
    5: 'tan',
    6: 'silver',
    7: 'green',
    8: 'white',
    9: 'pink',
    10: 'cyan',
    11: 'purple',
    12: 'greenyellow',
    13: 'orchid',
    14: 'darkgoldenrod',
    15: 'slateblue',
    16: 'black',
    17: 'y',
    18: 'yellowgreen',
    19: 'forestgreen',
    20: 'springgreen',
    21: 'darkturquoise',
    22: 'dodgerblue',
    23: 'mediumblue',
    24: 'darkblue',
    25: 'darkslateblue',
    26: 'blueviolet',
    27: 'magenta',
    28: 'fuchsia',
    29: 'mediumvioletred',
    30: 'crimson',
    31: 'peru',
    32: 'goldenrod'
}
colors = {
    'T': 'mediumseagreen',
    'E': 'white',
    'B': 'yellow',
    'H': 'cornflowerblue',
    'G': 'tomato',
    'I': 'violet',
    'C': 'black'
}


class forceExtension:
    def __init__(self, dataFile, pullingVelocity, pullingAxis = 'z', numpyFile=False, manualInputs=False, **kargs):
        if numpyFile:
            self.dataFile = dataFile
            forceTensor = np.load(dataFile)
            self.pullingVelocity = pullingVelocity
            self.t = forceTensor[:, 0]
            self.z = forceTensor[:, 1]
            self.F = forceTensor[:, 2]
            self.frame_t = self.t / (10**5)
            self.tInNs = self.t * 2 * (10**-6)

        #Note: this routine is only accessable if a dictionary is passed in
            #E.g. SMDCpd_i = objects.forceExtension(**SMDEDictC_i, dataFile = None, manualInputs=True)
        elif manualInputs:
            self.pullingVelocity = pullingVelocity
            self.t = kargs.get('t')
            self.z = kargs.get('z')
            self.F = kargs.get('F')
            self.frame_t = kargs.get('frame_t')
            self.tInNs = kargs.get('tInNs')


        else:
            self.dataFile = dataFile
            self.pullingVelocity = pullingVelocity
            self.t = functions.SMDLogReader(dataFile, 'timestep')
            self.F = functions.SMDLogReader(dataFile, 'force', pullingAxis = pullingAxis)
            self.z = self.t * self.pullingVelocity
            self.frame_t = self.t / (10**5)
            self.tInNs = self.t * 2 * (10**-6)

    def convolve(self, kernelSize=100):
        # smoothing data a bit
        self.F = np.convolve(self.F, np.ones(kernelSize) / kernelSize, mode='same')

    def returnDomainSpecificForceData(self, lowerBound, upperBound):
        domainSpecificForceData = self.F[(lowerBound < self.z) & (self.z < upperBound)]

        return domainSpecificForceData

    def temporalDataDebugRoutine(self):
        fig, ax = plt.subplots()
        ax.hist(self.t, bins=10000)
        plt.title("Your hopefully uniform distribution")
        plt.show()

    def dataRepeatFix(self):
        forceTensor = zip(self.t, self.z, self.F)

        fixedForceTensor = np.array([0, 0, 0])
        fixedForceTensor = fixedForceTensor[:, np.newaxis]

        timeStamps = set()
        with tqdm(total=len(self.t)) as progressBar:
            for row in forceTensor:
                if row[0] in timeStamps:
                    progressBar.update(1)
                    continue
                else:
                    timeStamps.add(row[0])
                    row = np.array(row)[:, np.newaxis]
                    fixedForceTensor = np.append(fixedForceTensor, np.array(row), axis=1)
                    progressBar.update(1)
        self.t = fixedForceTensor[0, :]
        self.z = fixedForceTensor[1, :]
        self.F = fixedForceTensor[2, :]

    def forceTensorToNumpyArrays(self, filename):
        t = self.t[:, np.newaxis]
        z = self.z[:, np.newaxis]
        F = self.F[:, np.newaxis]
        forceTensor = np.append(t, z, axis=1)
        forceTensor = np.append(forceTensor, F, axis=1)
        np.save(filename, forceTensor)

    def generateCleanNumpyFileForForceExtension(self, numpyFileName):
        self.dataRepeatFix()
        self.forceTensorToNumpyArrays(numpyFileName)

    def readForceTensorNumpyFile(self, filename):
        forceTensor = np.load(filename)
        self.t = forceTensor[:, 0]
        self.z = forceTensor[:, 1]
        self.F = forceTensor[:, 2]

    def AUC(self, limitsOfIntegration):

        if limitsOfIntegration[0] > limitsOfIntegration[-1]:
            raise Exception("The order of your limits of integration is reversed")

        lowerLimitOfIntegration = limitsOfIntegration[0]
        upperLimitOfIntegration = limitsOfIntegration[1]

        lowerLimitIndex = functions.findClosestValueIndex(self.z, lowerLimitOfIntegration)
        upperLimitIndex = functions.findClosestValueIndex(self.z, upperLimitOfIntegration)

        changeInZ = ((self.z[upperLimitIndex] - self.z[lowerLimitIndex]) /
                     (upperLimitIndex - lowerLimitIndex))

        riemannSum = np.sum(self.F[lowerLimitIndex:upperLimitIndex] * changeInZ)

        return riemannSum

    def elapsedDataDecorator(method):
        def wrapper(self, frameNum):
            zMaximum = ((100000) * self.pullingVelocity * frameNum)

            z, F = method(self, frameNum)

            if frameNum > 2:
                zProgress = z[z <= zMaximum]
                forceDataProgress = F[0:len(zProgress)]

                if np.sum(forceDataProgress) > 0:
                    totalEnergy = \
                        round(np.sum(forceDataProgress * (1000) * self.pullingVelocity), 3)
                    totalEnergy = format(totalEnergy, ".3f")

            else:
                zProgress = z[z <= zMaximum]
                forceDataProgress = F[0:len(zProgress)]

            # Returning the elapsed times
            return zProgress, forceDataProgress

        return wrapper

    def dF(self):

        step_length = 1000  
        dFs = []          
        z_values = []  
        for n in range(0, len(self.z) - step_length, step_length):
            differential = (self.F[n + step_length] - self.F[n]) / (self.z[n + step_length] - self.z[n])
            dFs.append(differential)
            z_values.append(self.z[n])

        z = np.array(z_values)

        return z, dFs

    @elapsedDataDecorator
    def elapsedForceExtension(self, frameNum):
        return self.z, self.F

    @elapsedDataDecorator
    def elapsed_dF(self, frameNum):
        z, dFs = self.dF()
        return z, dFs

    def shortenMaxDistance(self, newDistance, redefineAttributes = True):
        frameNum = newDistance/((100000) * self.pullingVelocity)
        z, F = self.elapsedForceExtension(frameNum)
        t = self.t[0:len(z)]
        if redefineAttributes:
            self.z = z
            self.t = t
            self.F = F
        return z, t, F

    def findPeaks(self, prominance=25, width=10, plot=False):

        widthInSamples = len(self.shortenMaxDistance(width, redefineAttributes=False)[0])
        peaks, _ = find_peaks(self.F, prominence=prominance, width = widthInSamples)

        self.peakAbcissa = self.z[peaks]
        self.peakOrdinate = self.F[peaks]

        if plot:
            plt.plot(self.z, self.F)
            plt.scatter(self.peakAbcissa, self.peakOrdinate, color='red', marker='*')
            plt.show()

        return self.peakAbcissa, self.peakOrdinate

     
    def removeSnippets(self, user_time, time_lag):
        # Calculate the start and end time for the range to be removed
        start_time = user_time - time_lag
        end_time = user_time

        # Find the indices corresponding to the start and end time
        indices_to_remove = np.where((self.tInNs >= start_time) & (self.tInNs <= end_time))[0]

        # Extract the snippets
        self.tSnippet = self.t[indices_to_remove]
        self.zSnippet = self.z[indices_to_remove]
        self.fSnippet = self.F[indices_to_remove]
        self.tInNsSnippet = self.tInNs[indices_to_remove] 
        self.frame_tSnippet = self.frame_t[indices_to_remove]        

        #Keep old values to avoid debacles
        self.tOriginal = self.t.copy()
        self.zOriginal = self.z.copy()
        self.FOriginal = self.F.copy()
        self.tInNsOriginal = self.tInNs.copy()
        self.frame_tOriginal = self.frame_t.copy()

        # Remove the values from the original attributes
        self.t = np.delete(self.t, indices_to_remove)
        self.z = np.delete(self.z, indices_to_remove)
        self.F = np.delete(self.F, indices_to_remove)
        self.tInNs = np.delete(self.tInNs, indices_to_remove)
        self.frame_t = np.delete(self.frame_t, indices_to_remove)

    def resetSnippets(self):
        self.t = self.tOriginal
        self.z = self.zOriginal
        self.F =  self.FOriginal
        self.tInNs = self.tInNsOriginal
        self.frame_t = self.frame_tOriginal


    def findLinearSnippetTrendLine(self):

        try:
            m, b, _, _, _ = linregress(self.zSnippet, self.fSnippet)
            self.linearF = m*self.zSnippet + b
            self.m = m
            return self.zSnippet, self.linearF, self.m

        except:
            print("subsection is empty, returning last values in the register")
            return self.zSnippet, self.linearF, self.m

    def findPiecewiseSnippetTrendline(self, lag=20, lowRes=False, lagType = 'Temporal', N = None):

        pieceWiseF = np.array([])
        lowResPieceWiseF = np.array([])
        lowResPieceWiseZ = np.array([])

        dF = np.array([])
        dZ = np.array([])

        if lagType == 'Temporal':

            t = 0
            for i in range(0, math.floor(max(self.tInNs/lag))):
                t+=lag
                self.removeSnippets(t, lag)
            
                if len(self.zSnippet) == 0 and len(self.fSnippet) == 0:
                    print("Assuming F(z) end has nearly reached")
                    continue
                zSnippet, F, m = self.findLinearSnippetTrendLine()
                dF = np.append(dF, np.array([m]))
                dZ = np.append(dZ, np.array([self.zSnippet[-1]]))
                pieceWiseF = np.append(pieceWiseF, F)
                lowResPieceWiseF = np.append(lowResPieceWiseF, F)
                lowResPieceWiseZ = np.append(lowResPieceWiseZ, np.copy(zSnippet))           
                self.resetSnippets()

            if lowRes:
            
                self.pieceWiseZ = lowResPieceWiseZ
                self.pieceWiseF = lowResPieceWiseF
                self.dZ = dZ
                self.dF = dF
        
        
                return self.pieceWiseZ, self.pieceWiseF

            else:

                self.pieceWiseF = pieceWiseF
                self.pieceWiseZ = self.z[0:self.pieceWiseF.shape[0]]
                self.dZ = dZ
                self.dF = dF
        
        
                return self.pieceWiseZ, self.pieceWiseF

        elif lagType == 'Spatial':
            
            raise Exception("This code is not configured to handle spatial cases")

        elif lagType == 'Resolution':

            dz = (self.z[-1] - self.z[0]) / N
            z = self.z[0]
            started = False
            self.pieceWiseZ = np.array([])
            self.pieceWiseF = np.array([])
            
            for i_z in range(0, N - 1):
                
                Z_I = np.where((self.z >= z) & (self.z <= z + dz))[0]

                z = z + dz               

                self.zSnippet = self.z[Z_I]
                self.fSnippet = self.F[Z_I]

                if not started:
                    zSnippet, linearF, _ = self.findLinearSnippetTrendLine()
                    self.pieceWiseZ = np.append(self.pieceWiseZ, zSnippet[0])
                    self.pieceWiseF = np.append(self.pieceWiseF, linearF[0])
                    self.pieceWiseZ = np.append(self.pieceWiseZ, zSnippet[-1])
                    self.pieceWiseF = np.append(self.pieceWiseF, linearF[-1])
                    started = True

                else:
                    try:
                        oldPieceWiseZ = self.pieceWiseZ
                        oldPieceWiseF = self.pieceWiseF

                        zSnippet, linearF, _ = self.findLinearSnippetTrendLine()
                        self.pieceWiseZ = np.append(self.pieceWiseZ, zSnippet[-1])
                        self.pieceWiseF = np.append(self.pieceWiseF, linearF[-1])

                    except:
                        print("Assuming we've snagged an empty vector by accident")
                        self.pieceWiseZ = np.append(self.pieceWiseZ, self.pieceWiseZ[-1])
                        self.pieceWiseF = np.append(self.pieceWiseF, self.pieceWiseF[-1])


                        

            return self.pieceWiseZ, self.pieceWiseF
       
    def find_dFWithPiecewise(self, lag=20):

        pieceWiseF = np.array([])
        dF = np.array([])
        dZ = np.array([])
        

        t = 0
        for i in range(0, math.floor(max(self.tInNs/lag))):
            t+=lag
            self.removeSnippets(t, lag)
            if len(self.zSnippet) == 0 and len(self.fSnippet) == 0:
                print("Assuming F(z) end has nearly reached")
                continue
            _, F, m = self.findLinearSnippetTrendLine()
            dF = np.append(dF, np.array([m]))
            dZ = np.append(dZ, np.array([self.zSnippet[-1]]))
            self.resetSnippets()

        self.dZ = dZ
        self.dF = dF
        
        return self.dZ, self.dF

    #Note: this function will misbehave if data is not convolved
    def shrinkOrdinate(self, maxAFMpeak, SMDprominance=25, SMDwidth=10, plot=False):
        peakAbcissa, peakOrdinate = self.findPeaks(prominance=SMDprominance, 
                                                   width=SMDwidth, plot=plot)
 
        maxLength = self.z[round(len(self.z) * 0.25)]
        peaks = zip(peakAbcissa, peakOrdinate)
        
        filteredPeaks = [peak for peak in peaks if peak[0] < maxLength]
        peakOrdinate = [peak[-1] for peak in filteredPeaks]
       
        maxSMDpeak = max(peakOrdinate)
        scalingFactor = maxSMDpeak/maxAFMpeak
        self.F_original = self.F
        self.F = self.F / scalingFactor

    def orderAlongTempDomain(self):
        
        combined = np.column_stack((self.t, self.z, self.F))
        sorted_combined = combined[np.argsort(combined[:, 0])]
        sorted_t, sorted_z, sorted_F = np.hsplit(sorted_combined, 3)

        self.t = sorted_t
        self.z = sorted_z
        self.F = sorted_F
    
        return sorted_t, sorted_z, sorted_F

    def checkZdataTypes(self):

        zDataTypes = set()
        zDataTypes.update(np.unique(self.z.dtype))
        return zDataTypes

class zDistribution:
    def __init__(self, residue_z_dict):
        self.residue_z_dict = residue_z_dict
        self.residues = list(residue_z_dict)
        self.zDistribution = np.array(
            [residue_z_dict[residue] for residue in self.residues])

    def histogramDecorator(method):
        def wrapper(self, normalize=True, *args, **kwargs):
            # I don't want my bins any  larger than an amino acid
            # residue
            binWidth = 4

            # Defining my histogram parameter
            n_bins = int(((max(self.zDistribution)
                           - min(self.zDistribution))
                          / binWidth)) + 1

            hist, binEdges = method(self, n_bins, *args, **kwargs)

            binEdges = binEdges[:-1]

            hist, binEdges = np.array(hist), np.array(binEdges)

            if normalize:
                norm = np.linalg.norm(hist)
                hist = np.array(hist) / norm

            return hist, binEdges
        return wrapper

    @histogramDecorator
    def histogramTransformation(self, n_bins):
        hist, binEdges = np.histogram(self.zDistribution,
                     n_bins,
                     (int(min(self.zDistribution)),
                      int(min(self.zDistribution))
                      + 4 * n_bins))
        return hist, binEdges

    @histogramDecorator
    def histogramTransformationWithSubtraction(self, n_bins, subtraction=2):
        #This will output faulty values in the event that there is no cluster
            #To be found
        hist, binEdges = np.histogram(self.zDistribution,
                     n_bins,
                     (int(min(self.zDistribution)),
                      int(min(self.zDistribution))
                      + 4 * n_bins))

        hist = np.array([max(discreteValueCount - subtraction, 0) for discreteValueCount in hist])
        return hist, binEdges

    def histogramFit(self, subtraction = 2):
        hist, bins = self.histogramTransformationWithSubtraction(
            normalize=True, subtraction=subtraction)

        gaussianX, gaussianY = functions.histogramFit(hist, bins)

        return gaussianX, gaussianY

    def criticalValues(self, p_value=0.05):
        gaussianX, gaussianY = self.histogramFit()

        #I'm making the somewhat heavy-handed supposition that
            # all values of x are qually spaced apart
        dX = gaussianX[1] - gaussianX[0]

        def criticalValueFinder(gaussianX, gaussianY, p_value = p_value):
            p = 0
            i = 0
            while p <= p_value:
                p+=dX*gaussianY[i]
                i+=1
            critical_X = gaussianX[i]
            return critical_X

        leftCriticalValue = criticalValueFinder(gaussianX, gaussianY)
        rightCriticalValue = criticalValueFinder(np.flip(gaussianX, axis=0),
                                                 np.flip(gaussianY, axis=0))

        return leftCriticalValue, rightCriticalValue

#TODO: this code would be slightly better if I could generate this
    #TODO: object directly from the appropriate file
    #TODO: but generating an instance of the structure
    #TODO: object will work for now
class secondaryStructure:
    def __init__(self, index, ssIndex, resIDs,
                 residueIndices, type, vmdColors = vmdColors):

        self.index = index
        self.descriptor = f"{type}{ssIndex}"
        self.color = vmdColors[index]
        self.resIDs = resIDs
        self.residueIndices = residueIndices
        self.ssType = type
        self.ssDetachmentPoint = None
        self.ssDetachmentFrame = None

    def isDetached(self, zAxisProjection, p = 0.05):
        zTransformation = self.resIDs
        residueLocationDict = zAxisProjection.residue_z_dict
        residueKeys = []
        for resID in self.resIDs:
            for residueKey in residueLocationDict.keys():
                if resID == residueKey._id[1]:
                    residueKeys.append(residueKey)
                    break
        residueLocations = [residueLocationDict[residueKey]
                            for residueKey in residueKeys]

        leftCriticalValue, rightCriticalValue = zAxisProjection.\
            criticalValues(p_value= p )

        def isOutsideCriticalValues(residueLocation,
                                   leftCriticalValue,
                                   rightCriticalValue):
            if residueLocation >= leftCriticalValue and\
                residueLocation <= rightCriticalValue:
                return False
            else:
                return True

        criticalValueBoolean = lambda residueLocation: \
            isOutsideCriticalValues(residueLocation,
            leftCriticalValue=leftCriticalValue,
            rightCriticalValue=rightCriticalValue)

        isDetachedList = map(criticalValueBoolean, residueLocations)

        #I somewhat arbitrarily define detachment as the point
            #at which all values fall outside the critical value
        if False not in isDetachedList:
            isDetached_finalSS_bool = True
        else:
            isDetached_finalSS_bool = False

        self.isDetached_SS_bool = isDetached_finalSS_bool
        return isDetached_finalSS_bool    

    def addUnfoldingInfo(self):

        residueUnfoldingCounts = [0 for index in self.residueIndices]
        self.residueUnfoldingCounts = residueUnfoldingCounts

        residueUnfoldingProcesses = [False for index in self.residueIndices]
        self.residueIsUnfolding = residueUnfoldingProcesses

        residueCompleteUnfoldings = [False for index in self.residueIndices]
        self.residueCompleteUnfoldings = residueCompleteUnfoldings

        residueStartedUnfolding = [False for index in self.residueIndices]
        self.residueStartedUnfolding = residueStartedUnfolding


    def addLabel(self):
        self.label = secondaryStructureLabel(self)

    def estimate_r(self, zAxisProjection):
        zTransformation = self.resIDs
        residueLocationDict = zAxisProjection.residue_z_dict
        residueKeys = []
        for resID in self.resIDs:
            for residueKey in residueLocationDict.keys():
                if resID == residueKey._id[1]:
                    residueKeys.append(residueKey)
                    break
        residueLocations = [residueLocationDict[residueKey]
                            for residueKey in residueKeys]
        r_z = np.mean(np.array(residueLocations))
        self.r_z = r_z      
        return r_z

class secondaryStructureLabel(secondaryStructure):
    def __init__(self, x, fontsize, yLim, ssCount,
                 secondary_structure_instance, shape = '.', size = 1):

        for key, value in vars(secondary_structure_instance).items():
            setattr(self, key, value)

        self.text = f"{self.descriptor}"
        self.fontsize = fontsize
        self.shape = shape
        self.yLim = yLim
        self.ssCount = ssCount
        self.size = size
        self.x = x

        plottingWidth = (yLim[-1] - yLim[0])/2

        plottingIncrement = plottingWidth/ssCount
        self.y = 0.2*((yLim[0] + plottingWidth) + plottingIncrement*(self.index - 1))
 
class ssDecomposition:
    def __init__(self, pathTo_ssDecomposition,
                 structureObject,
                 pullingVelocity,
                 vmdColors=vmdColors,
                 colors=colors):

        # Routine for collecting secondary structure data
        with open(pathTo_ssDecomposition, 'r') as f:
            data = np.array([x.split() for x in f.readlines() if "#" not in x])

        # Separate ss data into categories
        categories = ['T', 'E', 'B', 'H', 'G', 'I', 'C']
        category_indices = []
        for i in range(len(data)):
            if data[i][4] in categories:
                category_indices.append(data[i][4])

        self.z = data[:, 3].astype(int) * 100000 * pullingVelocity
        self.y = data[:, 0].astype(int)
        self.colorCategories = \
            [colors[color] for color in category_indices]
        self.secondaryStructures = structureObject.\
            getSecondaryStructures()
        self.pullingVelocity = pullingVelocity
        self.structure = structureObject
        self.legendLabels = ["Turn", "Extended Configuration", "Isolated Bridge", "Alpha Helix",
                             "3-10 Helix", "Pi-helix", "Coil (other)"]

    def getTickLabels(self, fontsize = 10):

        structureTickQuants = [round(statistics.median(ssIndices.resIDs))
                               for ssIndices in
                               self.secondaryStructures]
        tickColors = [secondaryStructure.color
                      for secondaryStructure
                      in self.secondaryStructures]

        self.structureTickQuants = structureTickQuants
        self.tickColors = tickColors

        self.fontSize = 10
        self.decompositionLabels = [structure.descriptor for
                       structure in self.secondaryStructures]


        return structureTickQuants, \
            tickColors, self.decompositionLabels

    def unfoldingDescriptor(self, tau=50):
        resnum = len(list(self.structure.getResidues_shortcut()))
        frameNum = int(len(self.z)/resnum)

        #C is just the name of my matrix with color labels
        C_reshape = np.array(self.colorCategories).reshape((frameNum,resnum))

        ssStructureIndices_dict = {ss.descriptor:ss.resIDs
                                   for ss in self.secondaryStructures}


        list_of_rows = []
        for row in C_reshape:
            ss_unfoldingDict = {}
            for ss in ssStructureIndices_dict.keys():
                ss_unfoldingDict[ss] = row[
                                       ssStructureIndices_dict[ss][0]:
                                       ssStructureIndices_dict[ss][-1]]
            list_of_rows.append(ss_unfoldingDict)

        colorDF = pd.DataFrame(list_of_rows)

        for secondaryStructure in self.secondaryStructures:
            secondaryStructure.addUnfoldingInfo()

        def findUnfoldings(column, tau=tau):

            ssStructure = [ssStructure for ssStructure in
                           self.secondaryStructures if
                           ssStructure.descriptor == column.name][0]

            unfoldingPerResidueColumn = []


            for i in range(len(column) - 1):
                comparingArrays = list(zip(column[i], column[i + 1]))
                unfoldingPerResidueColumn_j = []
                for j in range(0, len(comparingArrays)):
                    column_n = comparingArrays[j][0]
                    column_n_plus_1 = comparingArrays[j][1]

                    if ssStructure.residueCompleteUnfoldings[j]:
                        unfoldingPerResidueColumn_j.append(0)
                        continue
                    else:
                        pass
                    if ssStructure.residueUnfoldingCounts[j] == tau:
                        ssStructure.residueCompleteUnfoldings[j] = True
                        unfoldingPerResidueColumn_j.append(1)
                        continue

                    else:
                        unfoldingPerResidueColumn_j.append(0)


                    if column_n != 'black' and column_n_plus_1 == 'black':
                        ssStructure.residueStartedUnfolding[j] = True
                        ssStructure.residueUnfoldingCounts[j] += 1
                        continue

                    elif ssStructure.residueStartedUnfolding[j]:

                        if column_n == 'black' and column_n_plus_1 == 'black':
                            ssStructure.residueUnfoldingCounts[j] += 1
                            continue

                        else:
                            ssStructure.residueStartedUnfolding[j] = False
                            ssStructure.residueUnfoldingCounts[j] == 0
                            continue
                    else:
                        continue

                unfoldingPerResidueColumn.append(unfoldingPerResidueColumn_j)
            return unfoldingPerResidueColumn

        finalUnfoldingDF = pd.DataFrame()

        for column in colorDF.columns:
            unfoldingPerResidueColumn = findUnfoldings(colorDF[column])
            finalUnfoldingDF[column] = unfoldingPerResidueColumn

        self.unfoldingDescription = finalUnfoldingDF
        return finalUnfoldingDF

    def checkUnfoldingEvent(self, frame, unfoldingDescription):

        frameStructureUnfoldingData = unfoldingDescription.iloc[frame]

        unfoldingSSs = []
        for column in frameStructureUnfoldingData.index:

            if 1 in frameStructureUnfoldingData[column]:
                unfoldingSSs.append(column)

        return unfoldingSSs


    def ssSoFar(self, frameNum):

        proteinLength = len(list(self.structure.getResidues_shortcut()))
        z_elapsed = self.z[(frameNum - 1)*proteinLength:frameNum*proteinLength]
        Folding_elapsed = self.y[(frameNum - 1)*proteinLength:frameNum*proteinLength]
        colorCateogiresForMatplotlib = self.colorCategories[
            (frameNum - 1)*proteinLength:
            frameNum*proteinLength]

        return z_elapsed, Folding_elapsed, colorCateogiresForMatplotlib



class frameStructure:
    def __init__(self, pathToPDB, pathTo_ssFile_ID=None,
                 pathTo_ssFile_residues=None):
        self.pathToPDB = pathToPDB
        
        #Attributes usefull in ss analysis
        self.pathTo_ssFile_ID = pathTo_ssFile_ID
        self.pathTo_ssFile_residues = pathTo_ssFile_residues
        

    def getResidues_shortcut(self, chainName = "X"):
        p = parser()
        warnings.simplefilter("ignore", BiopythonWarning)

        #The identifier is technically different from the chain name
	#But because my use is synonymous  in this case, I see
	#no harm in proceeding with the chainName as the identifier

        structure = p.get_structure(chainName, self.pathToPDB)
        chains = structure.get_chains()
        
        for chain in chains:
            if chain.id == chainName:
                proteinChain = chain
                break

        residues = proteinChain.get_residues()
        return residues

    def getPullingLength(self):
        return 4 * \
            len(list(
                self.getResidues_shortcut()))

    def zAxisProjection(self, chainName = "X"):

        residues = self.getResidues_shortcut(chainName)

        resCoords = {}

        for res in residues:
            for atom in res.get_atoms():
                if "CA" in atom.fullname:
                    resCoords[res] = atom.coord[-1]

        frame_n_zDistribution = zDistribution(resCoords)

        return frame_n_zDistribution

    #This method is philistine, philistine I tell you
    def getSecondaryStructures(self):

        def alpha_beta_partitioning(ssFile):
            data_dict = {}

            with open(ssFile, 'r') as file:
                for line in file:
                    if '{' not in line:
                        parts = line.split(' ')
                    else:
                        parts = line.split('{')

                    key = int(parts[0].strip())
                    values = [int(x) for x in parts[1].split() if x.isdigit()]

                    data_dict[key] = values

            keys = list(data_dict.keys())
            partition_index = None

            # Find the index where the drop in values occurs
            for i in range(len(keys) - 1):
                current_values = data_dict[keys[i]]
                next_values = data_dict[keys[i + 1]]

                if current_values[-1] > next_values[0]:
                    partition_index = i
                    break

            if partition_index is not None:
                alphaHelices = {k: v for k, v in data_dict.items() if k <= keys[partition_index]}
                betaSheets = {k: v for k, v in data_dict.items() if k > keys[partition_index]}

            else:
                raise Exception(f"It looks like something is wrong with your ssStructure file, \n"
                                f"check {ssFile}")

            return alphaHelices, betaSheets

        alphaHelices_resID, betaSheets_resID = \
            alpha_beta_partitioning(self.pathTo_ssFile_ID)

        alphaHelices_residues, betaSheets_residues = \
            alpha_beta_partitioning(self.pathTo_ssFile_residues)

        def ssObjects(ssStructures_ID, ssStructures_residues, type):
            ssIndex = 0
            ssStructureObjects = []

            if ssStructures_ID.keys() != ssStructures_residues.keys():
                raise Exception(f"Check {self.pathTo_ssFile_ID} and \n"
                                f"{self.pathTo_ssFile_residues} \n"
                                f"Your structure keys have a mismatch, they are \n"
                                f"{ssStructures_ID.keys()} and \n"
                                f"{ssStructures_residues.keys()}")

            for key in ssStructures_ID.keys():
                ssIndex+=1

                ssStructureObjects.append(
                    secondaryStructure(index = key, ssIndex=ssIndex,
                                       resIDs=ssStructures_ID[key],
                                       residueIndices=
                                       ssStructures_residues[key],
                                       type = type))
            return ssStructureObjects

        alphaHelixObjects = ssObjects(alphaHelices_resID,
                                      alphaHelices_residues, 'α')
        betaSheetObjects = ssObjects(betaSheets_resID,
                                     betaSheets_residues,
                                     'β')

        self.ssStructures = alphaHelixObjects + betaSheetObjects
        self.alphaHelices = alphaHelixObjects
        self.betaSheets = betaSheetObjects

        return alphaHelixObjects + betaSheetObjects

    def frame_SSr_DebugRoutine(self, chainName = "X", i = 0, p_value = 0.05,
                               saveDir = "./tmp"):
        
        #Generating the usual nuts and bolts of p(z)
        zDistribution = self.zAxisProjection(chainName = chainName)
        hist, binEdges = zDistribution.histogramTransformationWithSubtraction()
        gaussianX, gaussianY = zDistribution.histogramFit()

        plt.clf()

        ssObjects = self.getSecondaryStructures()
        for ssObject in ssObjects:
            ssObject.estimate_r(zDistribution)
            r_z = ssObject.r_z
            color = ssObject.color
            label = ssObject.descriptor
            
            plt.axvline(x=r_z, color=color, linewidth=2, label=label)
            plt.text(r_z, plt.gca().get_ylim()[1] * 0.9, 
                     label, color=color, fontsize=10, horizontalalignment='left')

        leftCrit, rightCrit = zDistribution.criticalValues(p_value = p_value)
        plt.axvline(x=leftCrit, color='black', linewidth=4)
        plt.axvline(x=rightCrit, color='black', linewidth=4)
        
        plt.title(f'Debug frame {i}')
        plt.bar(binEdges, hist)
        plt.plot(gaussianX, gaussianY)
        plt.xlabel(r"$r_z$ (Å)")
        plt.ylabel("p($r_z$)")
        plt.xlim(0, 500)
        plt.ylim(0, 1)
        plt.savefig(f'{saveDir}/debugFrame{i}.png')


class movingResidue:
    def __init__(self, biopythonResidue, t, expectedPosition = np.array([0.00, 0.00, 0.00]),
                 tPosition = np.array([0.00, 0.00, 0.00]), position_d = np.array([0.00, 0.00, 0.00])):
        self.biopythonResidue = biopythonResidue
        self.expectedPosition = expectedPosition
        self.tPosition = tPosition
        self.position_d = position_d
        self.t = t

    def differentiateResidue(self):
        self.position_d = np.subtract(self.tPosition, self.expectedPosition)

class residuePair:

    def __init__(self, residue_i, residue_j):
        self.residue_i = residue_i
        self.residue_j = residue_j

    def getCovarianceMatrix(self):
        d_ri = self.residue_i.position_d[:, np.newaxis]
        d_rj = self.residue_j.position_d[:, np.newaxis]
        #I need to transpose the j residue
        d_rj = np.transpose(d_rj)
        M_ij = np.outer(d_ri,d_rj)
        self.M_ij = M_ij
        return M_ij

    def getNormalizedCorrelationMatrix(self):
        self.getCovarianceMatrix()
        M_mm = np.diagonal(self.M_ij)[:, np.newaxis]
        M_nn = np.diagonal(self.M_ij)[:, np.newaxis]
        M_nn = np.transpose(M_nn)

        M_diagonal_square = np.outer(M_mm, M_nn)
        M_diagonal_square_sqrt = np.sqrt(M_diagonal_square)

        D_ij = np.divide(self.M_ij, M_diagonal_square_sqrt)
        return D_ij

#I rather dislike this analytical approach, oh well, consider it a cautionary tale
	#for future projects, I'll need to stick with PDBs until further notice
class simulationAsPDBs:
    def __init__(self, allSimStructureDir):
        self.allSimStructureDir = allSimStructureDir

    def getDynamicalAnalysisFrame(self, t, tau = 5 ):
        if t < tau:
            raise Exception("You goofed, you can't probe the expected value before"\
                            "the simulation started")
        else:           
            PDBs = sorted(glob.glob(f"{self.allSimStructureDir}/*pdb"))
            PDBsToAnalyze = PDBs[t-tau:t]
            structuresToAnalyze = [frameStructure(PDB) for PDB in PDBsToAnalyze]

            #I figured it be handy to keep the time lag
                #Associated with the object
            self.tau = tau

            return structuresToAnalyze
    def getExpectedPositions(self, t, getFinalPositions = True):
        structuresToAnalyze = self.getDynamicalAnalysisFrame(t)
        myResidues = [movingResidue(biopythonResidue, t) for biopythonResidue in
                        structuresToAnalyze[0].getResidues_shortcut()]
        for structure in structuresToAnalyze:
            
            frameResidues = structure.getResidues_shortcut()
            
            for myResidue, biopythonResidue in zip(myResidues, frameResidues):
                #TODO: it's a little sloppy to have repetative code, but not especially problematic
                    #TODO: here, fix it if you can
                CA = next((atom for atom in biopythonResidue.get_atoms() if atom.name == "CA"))
                myResidue.expectedPosition = np.add(np.array(CA.coord), myResidue.expectedPosition)

        for myResidue in myResidues:
            myResidue.expectedPosition = myResidue.expectedPosition / self.tau

        if getFinalPositions:
            tResidues = structuresToAnalyze[-1].getResidues_shortcut()

            for myResidue, biopythonResidue in zip(myResidues, tResidues):
                CA = next((atom for atom in biopythonResidue.get_atoms() if atom.name == "CA"))
                myResidue.tPosition = np.array(CA.coord)


        return myResidues

    def get_dPositions(self, t):
        myResidues = self.getExpectedPositions(t)
        for residue in myResidues:
            residue.differentiateResidue()

        return myResidues

    def pairResidues(self, t):
        myResidues = self.get_dPositions(t)

        residuePairs_ij = []

        for residue_i in myResidues:
            residuePairs_i = []

            for residue_j in myResidues:
                residuePairs_i.append(residuePair(residue_i, residue_j))

            residuePairs_ij.append(residuePairs_i)

        return residuePairs_ij

    def getCovrianceMatrices(self, t):

        residuePairs_ij = self.pairResidues(t)

        for row in residuePairs_ij:
            for residue_ij in row:
                print(residue_ij.getCovarianceMatrix())

                raise Exception("SOTP")

    def getNormalizedCorrelationMatrix(self, t):

        residuePairs_ij = self.pairResidues(t)

        for row in residuePairs_ij:
            for residue_ij in row:
                print(residue_ij.getNormalizedCorrelationMatrix())

                raise Exception("SOTP")

class pearsonCorrelationTensor:
    def __init__(self, pathToNumpyFiles):
        self.dir = pathToNumpyFiles
        self.files = functions.sortFilePaths(glob(f"{self.dir}/*npy"))

    def getFrame(self, i):
        return np.load(self.files[i])

    def getFrameEdges(self, i):
        frame = np.load(self.files[i])
        edgify = lambda C: -np.log(abs(C))
        edgeVectorization = np.vectorize(edgify)
        edges = edgeVectorization(frame)
        return edges

    def getNeighboringMatrixFileNames(self, i, n):
        if i < n:
            raise Exception("You goofed, you cannot have a window size which is larger than your"
                            " initial vector")

        neighbors = [self.files[fileIndex] for fileIndex in range(i - n, i + n)]

        return neighbors

    def getNeighboringMatrices(self, i, n):
        neighbors = self.getNeighborsMatrixFileNames(i, n)
        neighboringMatrices = [np.load(neighbor) for neighbor in neighbors]
        return neighboringMatrices

    def getMovingAverage(self, i, n):
        neighboringMatrices = self.getNeighboringMatrices(i, n)
        movingAverage = np.mean(neighboringMatrices, axis=0)

        return movingAverage

#TODO: this class should probably be pulled
class trajectoryDirWithForceExtension:
    def __init__(self, t, realpath, forceExtensionObject):
        self.t = t
        self.realpath = realpath
        self.forceExtension = forceExtensionObject
        self.fileName = realpath.split("/")[-1]

    def plotPickeling(self):
        fig, ax = plt.subplots()
        ax.plot(self.forceExtension.tInNs, self.forceExtension.F, color='blueviolet')
        ax.plot(self.forceExtension.tInNsSnippet, self.forceExtension.fSnippet, color='black')
        ax.set_title(f"NLeC Rep 1 @ {self.t} ns")
        ax.set_xlabel("t (ns)")
        ax.set_ylabel("F (pN)")
        ax.legend(labels=["F(z)", "20ns frame"])

        # Save the figure
        with open(f'{self.realpath}/plot.fig.pickle', 'wb') as f:
            pickle.dump(fig, f)

    def plotSaving(self, fileName, pieceWiseZ=np.array([]), pieceWiseF=np.array([])):
        fig, ax = plt.subplots()
        ax.plot(self.forceExtension.z, self.forceExtension.F, color='blueviolet')
        ax.plot(self.forceExtension.zSnippet, self.forceExtension.fSnippet, color='black')
        ax.set_title(f"{fileName} @ {self.t} ns")
        ax.set_xlabel("pulling distance along z-axis (Å)")
        ax.set_ylabel("F (pN)")
        
        if pieceWiseZ != np.array([]) and pieceWiseF != np.array([]):
            ax.plot(pieceWiseZ, pieceWiseF, color = 'white')

        ax.legend(labels=["F(z)", "20ns frame", "F(z) approximation"])


        # Save the figure
        #plt.savefig(f'{self.realpath}/plot.png')
        plt.show()


    def plotSavingWithDF(self, fileName, color, dIndex,
                         pieceWiseZ=np.array([]), pieceWiseF=np.array([]),
                         dZ=np.array([]), dF=np.array([])):
       fig, axes = plt.subplots(2, 1)
       axes[0].plot(self.forceExtension.z, self.forceExtension.F, color=color)
       axes[0].plot(self.forceExtension.zSnippet, self.forceExtension.fSnippet, color='black')
       axes[0].set_title(f"{fileName} @ {self.t} ns")
       axes[0].set_xlabel("pulling distance along z-axis (Å)")
       axes[0].set_ylabel("F (pN)")
        
       if pieceWiseZ != np.array([]) and pieceWiseF != np.array([]):
           axes[0].plot(pieceWiseZ, pieceWiseF, color = 'white')

       axes[0].legend(labels=["F(z)", "20ns frame", "F(z) approximation"])
       
       indices_to_remove = dIndex
       
       dZsnippet = dZ[indices_to_remove]
       dFsnippet = dF[indices_to_remove]
            
       axes[1].scatter(dZ, dF, color=color)
       axes[1].plot(dZ, dF, color=color)
       axes[1].set_xlabel("pulling distance along z-axis (Å)")
       axes[1].set_ylabel("stiffness, ⟨∂F/∂z⟩ (pN/Å)")

       axes[1].scatter(dZsnippet, dFsnippet, color='black')
       axes[1].plot(dZsnippet, dFsnippet, color='black')


       
       # Save the figure
       plt.savefig(f'{self.realpath}/plot.png')
       #plt.show()     

class pearsonCorrelationTensor:
    def __init__(self, pathToNumpyFiles):
        self.dir = pathToNumpyFiles
        self.files = functions.sortFilePaths(glob(f"{self.dir}/*npy"))

    def getFrame(self, i):
        return np.load(self.files[i])

    def getFrameEdges(self, i):
        frame = np.load(self.files[i])
        edgify = lambda C: -np.log(abs(C))
        edgeVectorization = np.vectorize(edgify)
        edges = edgeVectorization(frame)
        return edges

    def getNeighboringMatrixFileNames(self, i, n):
        if i < n:
            raise Exception("You goofed, you cannot have a window size which is larger than your"
                            " initial vector")

        neighbors = [self.files[fileIndex] for fileIndex in range(i - n, i + n)]

        return neighbors

    def getNeighboringMatrices(self, i, n):
        neighbors = self.getNeighborsMatrixFileNames(i, n)
        neighboringMatrices = [np.load(neighbor) for neighbor in neighbors]
        return neighboringMatrices

    def getMovingAverage(self, i, n):
        neighboringMatrices = self.getNeighboringMatrices(i, n)
        movingAverage = np.mean(neighboringMatrices, axis=0)

        return movingAverage

class residueNetwork:
    def __init__(self, contactDF, pearsonFrameEdges):
        contactDF.to_numpy()
        self.networkMatrix = np.array(contactDF*pearsonFrameEdges)

        # Create an empty graph
        G = nx.Graph()

        # Add nodes
        for i in range(len(self.networkMatrix)):
            G.add_node(i)


        # Add weighted edges
        for i in range(len(self.networkMatrix)):
            for j in range(i + 1, len(self.networkMatrix[i])):
                weight = self.networkMatrix[i][j]
                if weight > 0:
                    G.add_edge(i, j, weight=weight)

        self.network = G

    def draw(self):
        pos = nx.spring_layout(self.network)
        nx.draw(self.network, pos, with_labels=True, node_color='purple', font_weight='bold')
        plt.show()

    def drawCircle(self):
        pos = nx.circular_layout(self.network)
        nx.draw(self.network, pos, with_labels=True, node_color='purple')
        plt.show()

    def dijkstrasPath(self, startPoint, endPoint):
        path = nx.dijkstra_path(self.network,
                                startPoint, endPoint)
        return path

   def generateVMDoptimalPathsFile(self, startPoint, endPoint, savePath):
       raise Exception("This method has not been configured yet")

   def generateVMDdynamicNetworkFile(self, savePath):
       raise Exception("This method has not been configured yet")

   def getLargestConnectedComponent(self)
       if True:
           raise Exception("This method has not been configured yet")

       else:
           return 0, 100

    

#Looking a little incomplete

class VMDDNAoutput:
    def __init__(self, contactCorrDataFile):
        self.contactCorrDataFile = contactCorrDataFile
        self.contactCorrMatrix = functions.parse_matrix_file(self.contactCorrDataFile)
        with np.errstate(divide='ignore', invalid='ignore'):
            edgeMatrix = -np.log(np.abs(self.contactCorrMatrix))
        edgeMatrix[np.isinf(edgeMatrix)] = 0
        self.edgeMatrix = edgeMatrix

        # Create an empty graph
        G = nx.Graph()

        # Add nodes
        for i in range(len(self.edgeMatrix)):
            G.add_node(i)

            # Add weighted edges
        for i in range(len(self.edgeMatrix)):
            for j in range(i + 1, len(self.edgeMatrix[i])):
                weight = self.edgeMatrix[i][j]
                if weight > 0:
                    G.add_edge(i, j, weight=weight)

        self.network = G

    def draw(self):
        pos = nx.spring_layout(self.network)
        nx.draw(self.network, pos, with_labels=True, node_color='purple', font_weight='bold')
        plt.show()

    def drawCircle(self, show=True, saveDir=None, title=None):
        pos = nx.circular_layout(self.network)
        plt.clf()
        nx.draw(self.network, pos, with_labels=True, node_color='purple')

        # Enforce square plot
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis('off')  # Turn off the axis


        if show:
            plt.show()
        elif not show:
            print("The code presumes if you don't show the figure, it must save it")
            plt.title(title)
            plt.savefig(f"{saveDir}/{title}.png")


    def allSimplePaths(self, startingPoint, endingPoint):
        paths = nx.all_shortest_paths(self.network,
                                startingPoint, endingPoint)
        return paths

    def shortestPath(self, startingPoint, endingPoint):
        path = nx.shortest_path(self.network,
                            startingPoint, endingPoint)
        return path
    def generateNewContactDataFileWithScaling(self, scale=1):
        matrix = scale*self.contactCorrMatrix
        matrixAsString = functions.array_to_string(matrix)
        header = open(self.contactCorrDataFile).readlines()[0]
        file = header + matrixAsString
        file = file.replace("0.0", "0")
        with open(f"{self.contactCorrDataFile}.test.dat", "w") as corrFile:
            corrFile.write(file)

    def getConnectedComponents(self):
        connectedComponents = nx.connected_components(self.network)
        return connectedComponents

    def getLargestConnectedComponent(self):
        largest_cc = max(self.getConnectedComponents(), key=len)

        self.minFromLargestCC = min(largest_cc)
        self.maxFromLargestCC = max(largest_cc)

        return largest_cc
    
    def getDegreeDistribution(self):
        self.degreeDistribution = np.array([node[-1] for node in self.network.degree()])
        return self.degreeDistribution

    def getLargestConnectedComponentDegreeDistribution(self):
        largest_cc = list(self.getLargestConnectedComponent())
        self.getDegreeDistribution()
        largest_ccDegreeDistribution = self.degreeDistribution[largest_cc]
        return largest_ccDegreeDistribution

class suboptimalPathFile:
    def __init__(self, VMDDNA, leftResID, rightResID):
        self.shortestDistance = ""
        self.edgeLengthOffset = ""
        self.suboptimalPathMax = ""
        self.sources = leftResID
        self.targets = rightResID
        self.finalPaths = list(VMDDNA.allSimplePaths(leftResID, rightResID))
        self.numberOfPaths = len(self.finalPaths)
        self.residuesThatOccurMost = []
        self.edgesThatOccurMost = []

    def generate(self, saveDir):
        shortestDistance = [f"The shortest distance is here {self.shortestDistance}"]
        edgeLengthOffset = [f"The edge length offset is {self.edgeLengthOffset}"]
        suboptimalPathMax = [f"Suboptimal paths will be no longer than {self.suboptimalPathMax}"]
        sources = [f"The sources are: {self.sources}"]
        targets = [f"The targets are: {self.targets}"]
        finalPaths = ["The final paths are:"] + self.finalPaths
        numberOfPaths = [f"Number of paths is {self.numberOfPaths}"]

        residuesThatOccurMost = ([f"The residues that occur the most in the suboptimal paths are:"]
                                 + self.residuesThatOccurMost)

        edgesThatOccurMost = ([f"The edges that occur the most in the suboptimal paths are:"]
                              + self.edgesThatOccurMost)


        fileList = (shortestDistance + edgeLengthOffset + suboptimalPathMax + sources
                    + targets + finalPaths + numberOfPaths + residuesThatOccurMost + edgesThatOccurMost)

        with open(f"{saveDir}/network.out", "w") as file:
            for item in fileList:
                if isinstance(item, list):
                    file.write(", ".join(map(str, item)) + "\n")
                else:
                    file.write(str(item) + "\n")

        print("RANNNNNNNNNNNNNNNNNNNNNNNNNNN")


class contactNetwork:
    def __init__(self, contactDF):
        contactArray = contactDF.to_numpy()

        self.networkMatrix = contactArray

        # Create an empty graph
        G = nx.Graph()

        # Add nodes
        for i in range(len(self.networkMatrix)):
            G.add_node(i)


        # Add edges
        for i in range(len(self.networkMatrix)):
            for j in range(i + 1, len(self.networkMatrix[i])):
                weight = self.networkMatrix[i][j]
                if weight > 0:
                    G.add_edge(i, j, weight=weight)

        self.network = G

    def getDegreeDistribution(self):
        self.degreeDistribution = np.array([node[-1] for node in self.network.degree()])
        return self.degreeDistribution



class AFMData:
    def __init__(self, filename, replicate):
        self.filename = filename
        df = pd.read_csv(filename)
        
        self.replicate = replicate
        self.sep_vector = np.array(df['Sep'])
        self.F = np.array(df['Force'])
        self.d = self.sep_vector / (10 ** - 10)
        
                                        
    def plotSep_vs_Force(self, unconvertedVoltage=False):
        plt.figure(figsize=(10, 6))
        if unconvertedVoltage:
            plt.plot(self.d, self.V, label=f'Replicate {self.replicate}')
            plt.ylabel('Voltage (V)')
            plt.title('Distance vs. Voltage')

        else:
            plt.plot(self.d, self.F, label=f'Replicate {self.replicate}')
            plt.ylabel('Force (?)')
            plt.title('Distance vs. Force')

        plt.xlabel('Distance (Å)')
        plt.legend()
        plt.grid(True)
        plt.show()
    def plotVoltageFFT(self):
         fft_values = np.fft.fft(self.V)
    
         # Compute the corresponding frequencies
         n = len(self.V)
         timestep = np.mean(np.diff(self.sep_vector))  # Assuming uniform sampling
         freq = np.fft.fftfreq(n, d=timestep)
    
         # Plot the FFT values
         plt.figure(figsize=(10, 6))
         plt.plot(freq, np.abs(fft_values))
         plt.title('FFT of the Signal')
         plt.xlabel('Frequency (Hz)')
         plt.ylabel('Amplitude')
         plt.grid(True)
         plt.show()
    
    def findLowerQuartileData(self):
        
        lowerQuartileIndex = int(len(self.d) * 0.25)
        
        lowerQuartile_d = self.d[0:lowerQuartileIndex]
        lowerQuartile_F = self.F[0:lowerQuartileIndex]
        
        self.lowerQuartile_d = lowerQuartile_d
        self.lowerQuartile_F = lowerQuartile_F
        
        return lowerQuartile_d, lowerQuartile_F

    def findProteinPeaks(self, prominence = 0.000000000005, width = 10, plot=False):
        
        self.findLowerQuartileData()
        peaks, _ = find_peaks(self.lowerQuartile_F, prominence = prominence, width = width)      
        
        self.peakAbcissa = self.d[peaks]
        self.peakOrdinate = self.F[peaks]

        if plot:
            plt.style.use('dark_background')
            plt.figure(figsize=(10, 6))
            plt.plot(self.d, self.F, label=f'Replicate {self.replicate}', color='purple')
            plt.scatter(self.peakAbcissa, self.peakOrdinate, marker='X', color='red')
            plt.ylabel('Force (?)')
            plt.title('Distance vs. Force')
            plt.xlabel('Distance (Å)')
            plt.legend()
            plt.show()

        return self.peakAbcissa, self.peakOrdinate  

class crossCorrelationFit:
    #psi(z_i) is a cross correlation within the range of -1 to 1
    def __init__(self, AFMobject, SMDobject, fitted_SMDobject=None, ρ=None):      
        self.SMDobject = SMDobject
        self.AFMobject = AFMobject
        self.fitted_SMDobject = fitted_SMDobject
        self.ρ = ρ
        self.max_ρ = ρ

    def determineTranslationPolicy(self, splits, maxAFMLength=None):
        
        self.N = splits

        if maxAFMLength == None:
            d = (
                self.AMFobject.sep_vector[-1] - 
                self.AFMobject.sep_vector[0]
                )
        else:

            d = maxAFMLength

        Δ = d/splits
        
        self.Δ = Δ 
        return Δ

    def updateAbcissa(self, Δ):
        self.fitted_SMDobject.z+=Δ
      
    def plot(self, trendLine=False):
        plt.clf()
        plt.plot(self.fitted_SMDobject.z, self.fitted_SMDobject.F)
        plt.plot(self.AFMobject.d, self.AFMobject.F)
        if trendLine:
            plt.plot(self.pieceWiseZ, self.pieceWiseF)
        plt.show()
            
    def scatterSMD(self):
        plt.scatter(self.fitted_SMDobject.z, self.fitted_SMDobject.F)
        plt.plot(self.AFMobject.d, self.AFMobject.F)
        plt.show()

    def findSMDTrendLine(self):
        AFM_x, AFM_y = self.AFMobject.d, self.AFMobject.F
        SMD_x, SMD_y = self.fitted_SMDobject.z, self.fitted_SMDobject.F

        AFM_CCindices = np.where( (AFM_x >= SMD_x[0]) & (AFM_x <= SMD_x[-1]) )

        AFM_y_CC = AFM_y[AFM_CCindices]
        AFM_x_CC = AFM_x[AFM_CCindices]

        #I cordially dislike this approach, but I see no alternatives which keep my code clean
        self.AFMalongSMD_F = AFM_y_CC

        AFMresolution = functions.resolutionCalculation(AFM_x_CC, AFM_y_CC)
           
        self.pieceWiseZ, self.pieceWiseF = self.fitted_SMDobject.findPiecewiseSnippetTrendline(
                                                 lag=1/AFMresolution, lowRes=False)
        
        return self.pieceWiseZ, self.pieceWiseF 
 
    #Note: the 'findSMDTrendline' method must be executed first
    def findPearsonCorrelation(self):
        try:
            ρ = functions.generalizedCrossCorrelation(self.pieceWiseF, self.AFMalongSMD_F)
        except:
            print('most likely, you failed to execute the .findSMDTrendLine method')
            ρ = functions.generalizedCrossCorrelation(self.pieceWiseF, self.AFMalongSMD_F)

    def showPearsonCorrelationPlot(self, matplotlibSaveDir, N, maxTranslocationLength, show=False, 
        showOptimalCurvePairing = False, returnP=False, saveOptimalSMDcurve=False):
        self.determineTranslationPolicy(splits = N, maxAFMLength = maxTranslocationLength)
        
        P = np.array([])
        Z = np.array([])
        z = 0
 
        max_ρValue = 0
        max_ρIndex = 0
        
        for i in tqdm(range(0, N)):
            self.updateAbcissa(self.Δ)
            self.findSMDTrendLine()
            ρ = functions.crossCorrelateSMDandAFMdata(self.fitted_SMDobject, self.AFMobject)
            P = np.append(P, ρ)
            Z = np.append(Z, z)
            z+=self.Δ

            if ρ > max_ρValue:
                max_ρValue = ρ
                max_ρIndex = i           

            plt.clf()
            plt.figure(figsize=(8, 6))
            plt.subplot(1, 2, 1)
            SMD, = plt.plot(self.fitted_SMDobject.z, self.fitted_SMDobject.F)
            AFM, = plt.plot(self.AFMobject.d, self.AFMobject.F)
            plt.xlabel('Pulling Distance (Å)')
            plt.ylabel('Force')
            plt.legend(handles=[SMD, AFM], labels=['SMD', 'AFM'])
            plt.title('Force-Extension curves')

            plt.subplot(1, 2, 2)
            plt.plot(Z, P)
            plt.xlim(0, self.Δ * N)
            plt.ylim(0,1)

            SMDdataFilename = os.path.basename(self.SMDobject.dataFile)
            AFMdataFilename = os.path.basename(self.AFMobject.filename)
                       
            plt.xlabel('Translation Distance (Å)')
            plt.ylabel('ρ')
            plt.title('ρ vs translation distance')

            plt.savefig(f"{matplotlibSaveDir}/{SMDdataFilename}_{AFMdataFilename}_ρ_at{round(z):04.0f}.png")
            
            if show:
                plt.show()
            self.P = P

            plt.close()
         
        #Resets to the original position
        self.fitted_SMDobject.z -= (N*self.Δ)
        #Places in the correct posiiton
        self.fitted_SMDobject.z += ((max_ρIndex + 1)*self.Δ)


        if saveOptimalSMDcurve:
            
            self.fitted_SMDobject.forceTensorToNumpyArrays(f"{matplotlibSaveDir}/{SMDdataFilename}_"
                                                           f"{AFMdataFilename}_optimal_ρ.npy")
            
        if showOptimalCurvePairing:

            plt.clf()
            plt.figure(figsize=(8, 6))
            plt.subplot(1, 2, 1)

            SMD, = plt.plot(self.fitted_SMDobject.z, self.fitted_SMDobject.F)
            AFM, = plt.plot(self.AFMobject.d, self.AFMobject.F)
            plt.xlabel('Pulling Distance (Å)')
            plt.ylabel('Force')
            plt.legend(handles=[SMD, AFM], labels=['SMD', 'AFM'])
            plt.suptitle(f'Force-Extension curves for {SMDdataFilename} & {AFMdataFilename}')

            plt.subplot(1, 2, 2)
            plt.plot(Z, P)
            plt.xlim(0, self.Δ * N)
            plt.ylim(0,1)
            plt.show()
            plt.close()
   
        if returnP:
            P = self.P
            return P
    
    #That policy needs to be found first
    def findP(self, matplotlibSaveDir, N, maxTranslocationLength):
        for i in tqdm(range(0, N)):
            self.updateAbcissa(self.Δ)
            self.findSMDTrendLine()
            ρ = functions.crossCorrelateSMDandAFMdata(self.fitted_SMDobject, self.AFMobject)
            P = np.append(P, ρ)
        self.P = P
        return P
            
    
       

