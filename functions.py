import objects

from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.ndimage import zoom

from Bio.SeqUtils import seq1
from Bio import PDB
import pandas as pd
import re
import csv

def histogramFit(hist, binEdges):

    # Function to fit a Gaussian PDF
    def gaussian_pdf(x, mean, std_dev, amplitude):
        return amplitude * norm.pdf(x, mean, std_dev)

    # Fit the Gaussian PDF to the histogram data
    p0 = [np.mean(binEdges), np.std(binEdges), np.max(hist)]  # Initial guess for parameters
    params, covariance = curve_fit(gaussian_pdf, binEdges, hist, p0=p0, maxfev=100000)

    #Returning our vectors of interest
    X = np.linspace(min(binEdges), max(binEdges), 100)
    Y = gaussian_pdf(X, *params)

    return X, Y

def hitogramData_and_histogramFit(PDBfile, ssDataFile):
    frame_0_structure = objects.frameStructure(PDBfile, ssDataFile)
    frame_0_zDistribution = frame_0_structure.zAxisProjection()
    hist, bins = frame_0_zDistribution.histogramTransformationWithSubtraction(
        normalize=True, subtraction=2)
    gaussianX, gaussianY = histogramFit(hist, bins)

    return hist, bins, gaussianX, gaussianY

def SMDLogReader(filename, SMDdatumType, pullingAxis = 'z'):
    data = []
    with open(filename, 'r') as file:
        for line in file.readlines():
            parts = line.split(" ")
            if parts[0] == 'ENERGY:' or parts[0] == 'GPRESSURE:' \
                    or parts[0] == 'PRESSAVG:' or parts[0] == 'ENERGSMD'\
                    or parts[0] == 'TIMING:' or parts[0] == 'GPRESSAVG:':
                continue
            if SMDdatumType == 'force' and pullingAxis == 'z':
                datum = float(parts[-1])
            
            elif SMDdatumType == 'force' and pullingAxis == 'y':
                datum = float(parts[-2])

            elif SMDdatumType == 'force' and pullingAxis == 'x':
                datum = float(parts[-3])
          

            elif SMDdatumType == 'timestep':
                datum = float(parts[2])
            else:
                raise Exception("You need to specify which kinda of data you would like \n"
                                "Acceptable values are \'force\' or \'timestep\'")

            data.append(datum)
    return np.array(data)

def parse_initial_ss_text_file(filename):
    # Initialize an empty dictionary
    data_dict = {}

    # Open the text file for reading
    with open(filename, 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Split the line into two parts using '{' as the separator
            if '{' not in line:
                parts = line.split(' ')
            else:
                parts = line.split('{')

            # Extract the key (the first digit) and the values (the list in brackets)
            key = int(parts[0].strip())
            values = [int(x) for x in parts[1].split() if x.isdigit()]

            # Add the key-value pair to the dictionary
            data_dict[key] = values

    return

def sortFilePaths(file_paths):
    return sorted(file_paths, key=lambda x: int(x.split('_')[-1].split('.')[0]))

def calculate_residue_distance_and_build_csv(structureFileName):

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", structureFileName)

    # Create an empty dictionary to store residue coordinates
    residue_coordinates = {}
    # I'll need this for later
    residue_identifiers = []

    # Iterate through the structure and extract coordinates for each residue
    for model in structure:
        for chain in model:
            for residue in chain:
                residueIdentifier = seq1(residue.get_resname())
                # I'm going to remove the number later, but for debugging
                # Purposes, I'd rather keep it there for now
                residue_key = f"{residueIdentifier}_{residue.id[1]}"
                residue_coordinates[residue_key] = residue["CA"].coord
                residue_identifiers.append(residue_key)

    # Now I'm generating an empty data frame, which I'll populate with the data
    data = [[None] * len(residue_identifiers) for _ in range(len(residue_identifiers))]
    distanceDF = pd.DataFrame(data, columns=residue_identifiers, index=residue_identifiers)

    # Calculate distances between residues
    # distances = {}
    for key1, coord1 in residue_coordinates.items():
        for key2, coord2 in residue_coordinates.items():
            distance = np.linalg.norm(np.array(coord1) - np.array(coord2))
            if distance > 10:
                # distances[(key1, key2)] = int(10)
                distanceDF.loc[key1, key2] = np.uint8(10)

            else:
                # distances[(key1, key2)] = int(distance)
                distanceDF.loc[key1, key2] = np.uint8(round(distance))

    return distanceDF

def findZeros(array, index):
    # Search backward
    for i in range(index, -1, -1):
        if array[i] == 0:
            backward_index = i
            break
    else:
        backward_index = 0

    # Search forward
    for i in range(index, len(array)):
        if array[i] == 0:
            forward_index = i
            break
    else:
        forward_index = len(array) - 1

    return backward_index, forward_index


def parse_matrix_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Ensure there's more than just a header
    if len(lines) < 2:
        raise ValueError("The file does not contain enough data to form a matrix.")

    # Skip the header line
    data_lines = lines[1:]

    # Parse the matrix
    matrix = []
    for line in data_lines:
        row = list(map(float, line.split()))
        matrix.append(row)

    # Convert the list of lists to a NumPy array
    matrix_np = np.array(matrix)

    # Verify the matrix is square
    if matrix_np.shape[0] != matrix_np.shape[1]:
        raise ValueError("The matrix is not square.")

    return matrix_np

# Function to extract the numerical part from the file name
def extract_number(file_path):
    match = re.search(r'(\d+)nsTrajectory', file_path)
    return int(match.group(1)) if match else None

def array_to_string(arr):
    rows = []
    for row in arr:
        row_str = " ".join(map(str, row))
        rows.append(row_str)
    return "\n".join(rows)


#Note: this is highly specialized to Kate's specific format, I would need to trust she will keep giving
	#Note: me files in precisely this  format
def parse_csv_to_afm_data(file_path, NleCFormat = False):
    # Read the CSV file while skipping the first row
    df = pd.read_csv(file_path, skiprows=1)
        
    # Extract unique filenames
    if not NleCFormat:
        repNames = [col for col in df.columns if 'Unnamed' not in col]
        forceIdentifiers = [col for col in df.columns if 'Unnamed' in col]
    else:
        repNames = [col for col in df.columns if 'Sep' in col]
        forceIdentifiers = [col for col in df.columns if 'Force' in col]

    
    # Initialize list to hold AFMData objects
    afm_data_list = []

    # Iterate over the columns in pairs (Sep, Force)
    for repName, forceIdentifier in zip(repNames, forceIdentifiers):
                
        sep_vector = df[repName].dropna().values[1:].astype('float64')
        force_vector = df[forceIdentifier].dropna().values[1:].astype('float64')
                        
        # Create an instance of AFMData and add to list
        afm_data = objects.AFMData(file_path, sep_vector, repName, force_vector)
        afm_data_list.append(afm_data)

    return afm_data_list

def extractColumnsFromAFM_dataConverted2CSV(file_path):
    sep = []
    V = []

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            # Append the first and second columns to the respective lists
            sep.append(float(row[0]))
            V.append(float(row[1]))

    # Convert lists to numpy arrays
    sep = np.array(sep)
    V = np.array(V)

    return file_path, sep, V

def plotForceExtensionsFromSMDandAFM(SMDforceExtensionObject, AFMData, xlimMax = 2000):
    AFM_x, AFM_y = AFMData.d, AFMData.F
    SMD_x, SMD_y = SMDforceExtensionObject.z, SMDforceExtensionObject.F


    plt.style.use('dark_background')
    plt.xlabel('Pulling Distance (Å)')
    plt.ylabel('Force (units)')
    afmLine, = plt.plot(AFM_x, AFM_y, color = 'red')
    smdLine, = plt.plot(SMD_x, SMD_y, color = 'blue')
    plt.xlim(-50, xlimMax)
    plt.legend(handles=[afmLine, smdLine], labels=['AFM', 'SMD'])
    plt.show()

def resolutionCalculation(d, F):
    # Sort the distance vector (if not already sorted)
    sorted_indices = np.argsort(d)
    sorted_distance = d[sorted_indices]
    sorted_force = F[sorted_indices]

    # Calculate intervals between consecutive distance points
    intervals = np.diff(sorted_distance)
 
    # Calculate the average interval
    average_interval = np.mean(intervals)
    
    # Calculate the resolution (readings per unit distance)
    resolution = 1 / average_interval
    return resolution

def generalizedCrossCorrelation(X_1, X_2):
    X_12dot = np.dot(X_1, X_2)
   
    X_1_norm = np.sqrt(np.linalg.norm(X_1)**2)
    X_2_norm = np.sqrt(np.linalg.norm(X_2)**2)

    CC = X_12dot/(X_1_norm * X_2_norm)
    return CC

#Figure out the deal with the regression line
def crossCorrelateSMDandAFMdata(SMDforceExtensionObject, AFMData, xlimMax = 2000, plot=False):
    
    #In the AFM experiments, dV != 0, so I will need to determine sampling resolution
    #In roughly the same region I will be cross correlating

    AFM_x, AFM_y = AFMData.d, AFMData.F
    SMD_x, SMD_y = SMDforceExtensionObject.z, SMDforceExtensionObject.F

    AFM_CCindices = np.where( (AFM_x >= SMD_x[0]) & (AFM_x <= SMD_x[-1]) )

    AFM_y_CC = AFM_y[AFM_CCindices]
    AFM_x_CC = AFM_x[AFM_CCindices]
       
    pieceWiseSMD_x, pieceWiseSMD_y = SMDforceExtensionObject.findPiecewiseSnippetTrendline(
                                     lowRes=True, 
                                     lagType = 'Resolution', 
                                     N = len(AFM_x_CC))
    
    #calculatingTheCrossCorrelation
    ρ = generalizedCrossCorrelation(pieceWiseSMD_y, AFM_y_CC)

            
    if plot==True:
              
        plt.xlabel('Pulling Distance (Å)')
        plt.ylabel('Force (pN)')
        smdLine, = plt.plot(SMD_x, SMD_y, color = 'blue')
        afmLine, = plt.plot(AFM_x, AFM_y, color = 'red')
        plt.xlim(-50, xlimMax)
        #plt.legend(handles=[afmLine, smdTrendLine, smdLine], labels=['AFM', 'SMD trendline', 'SMD'])
        plt.title(f"AFM replicate 9 with SMD replicate 1 and cross correlation {round(ρ, 3)}")
        plt.show()
    
    return ρ

def detachmentDebugCSV(ssStructures, saveDir = '/home/alpal/projects/SMD/gaussianIntegralPaperCode/debugTmp'):
    #The attributes you will need are .descriptor .ssDetachmentFrame
    data = []
    for ssStructure in ssStructures:
        data.append([ssStructure.descriptor, ssStructure.ssDetachmentFrame])

    df = pd.DataFrame(data, columns = ['ssName', 'frame of detachment'])
    df.to_csv(f'{saveDir}/ssDebug.csv')
        
        
        

###########################################################################################
#TODOnes
###########################################################################################

#TODO: figure out why these appear to cover radically different ranges

#TODO: figure out why the plot resolution alters when you change the domain

#TODO: Write a short routine to identify the smaller of the two resulting vectors
        #TODO: resize the data accordingly

#TODO: there is something egregiously wrong with my SMD data, I should seek to fix it

