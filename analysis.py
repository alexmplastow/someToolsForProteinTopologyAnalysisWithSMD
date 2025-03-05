import argparse
import sys
import functions
import objects
import matplotlib.pyplot as plt
import MDAnalysis as mda
import glob
import tqdm

def main(args):
    
    # Parameters
    forceExtensionFile = args.force_extension_file
    trajectoryFile = args.trajectory_file
    topologyFile = args.topology_file
    smdVelocity = float(args.smd_velocity)
    detachmentSearchDepth = int(args.detachment_search_depth)
    DCDsamplingResolution = int(args.dcd_sampling_resolution)
    ssDataFileResidue = args.ss_data_file_residue
    ssDataFileResnum = args.ss_data_file_resnum
    AFMfile = args.afm_file
    AFMprominence = float(args.afm_prominence)
    AFMwidth = float(args.afm_width)
    SMDprominance = float(args.smd_prominence)
    SMDwidth = int(args.smd_width)
    pdbRegen = args.pdb_regen
    p_value = float(args.p_value)
    AFMconversion = args.disable_afm_conversion
    plot = args.plot
    bimodal = args.bimodal
    
    if pdbRegen:
        #Collecting the frame to be used in the analysis
        u = mda.Universe(topologyFile, trajectoryFile)
        u.frameMax = len(u.trajectory)

        for t in range(0, u.frameMax, detachmentSearchDepth):

            u.trajectory[t]

            with mda.Writer(f"./tmp/frame_{t:04d}.pdb", u.atoms.n_atoms) as writer:
                writer.write(u)

        
    pdbs = sorted(glob.glob("./tmp/*pdb"))

    frameStructure_0 = objects.frameStructure(pdbs[0], ssDataFileResnum, ssDataFileResidue)
    sStructureNames = [name.descriptor for name in frameStructure_0.getSecondaryStructures()]
    sStructures = frameStructure_0.getSecondaryStructures().copy()
        
    AFMReplicate = objects.AFMData(AFMfile, 0)

    if AFMconversion == True:

        AFMReplicate.F = AFMReplicate.F/10e-12

    yLim = [0, max(AFMReplicate.F)]    

    detachedSSlabels = []


    for i, pdb in tqdm.tqdm(enumerate(pdbs)):

        #A little oldschool
        if i % detachmentSearchDepth != 0:
            continue

        frameStructure = objects.frameStructure(pdb, ssDataFileResnum, ssDataFileResidue)
        zAxisProjection = frameStructure.zAxisProjection()

        
        try:

            for sStructure in sStructures:

                if sStructure.isDetached(zAxisProjection, p = p_value, bimodal = args.bimodal):
       
                    if sStructure.ssDetachmentPoint == None:
                        sStructure.ssDetachmentPoint = float(smdVelocity)*i*DCDsamplingResolution
                        sStructure.ssDetachmentFrame = i

                        ssLabel = objects.secondaryStructureLabel(sStructure.ssDetachmentPoint,
                                                                        fontsize=9,
                                                                        yLim=yLim,
                                                                        ssCount=len(sStructures),
                                                                        secondary_structure_instance
                                                                        =sStructure)
                        detachedSSlabels.append(ssLabel)

        except:
            for sStructure in sStructures:
                if sStructure.ssDetachmentPoint == None:
                    sStructure.ssDetachmentPoint = float(smdVelocity)*detachmentSearchDepth*i*DCDsamplingResolution

    #Collecting the force extension curve from the SMD
    SMDforceExtension = objects.forceExtension(forceExtensionFile, float(smdVelocity), numpyFile=True)

    peakAbcissa, peakOrdinate = AFMReplicate.findProteinPeaks(prominence=AFMprominence, width=AFMwidth)

    highestAFMpeak = max(peakOrdinate)

    SMDforceExtension.shrinkOrdinate(highestAFMpeak, SMDprominance, SMDwidth=SMDwidth)

    functions.detachmentDebugCSV(sStructures)

    Δx = SMDforceExtension.z[0]

    #Adding SS data

    for label in detachedSSlabels:

        label.x+=Δx

        plt.text(label.x, label.y,
                 label.text,
                 fontsize=label.fontsize,
                 color=label.color)
        plt.scatter(label.x, label.y, color=label.color)


    functions.crossCorrelateSMDandAFMdata(SMDforceExtension, AFMReplicate, plot=plot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for analyzing SMD detachment events and AFM correlations.")
    parser.add_argument("--force_extension_file", required=True, help="Path to the force extension file.")
    parser.add_argument("--trajectory_file", required=True, help="Path to the trajectory file.")
    parser.add_argument("--topology_file", required=True, help="Path to the topology file.")
    parser.add_argument("--smd_velocity", default="0.000002", help="SMD pulling velocity.")
    parser.add_argument("--detachment_search_depth", type=int, default=10, help="Depth for detachment search.")
    parser.add_argument("--dcd_sampling_resolution", type=int, default=1000000, help="DCD sampling resolution.")
    parser.add_argument("--ss_data_file_residue", required=True, help="Path to the secondary structure residue file.")
    parser.add_argument("--ss_data_file_resnum", required=True, help="Path to the secondary structure residue number file.")
    parser.add_argument("--afm_file", required=True, help="Path to the AFM data file.")
    parser.add_argument("--afm_prominence", type=float, default=1.5e-15, help="AFM peak prominence.")
    parser.add_argument("--afm_width", type=float, default=0.1, help="AFM peak width.")
    parser.add_argument("--smd_prominence", type=float, default=1e-11, help="SMD peak prominence.")
    parser.add_argument("--smd_width", type=int, default=20, help="SMD peak width.")
    parser.add_argument("--pdb_regen", action="store_true", help="Regenerate PDB frames.")
    parser.add_argument("--p_value", type=float, default=0.01, help="P-value for detachment analysis.")
    parser.add_argument("--disable_afm_conversion", action="store_false", help="Convert AFM data scaling.")
    parser.add_argument("--plot", action="store_true", help="Plot the heatmap of SMD-AFM correlations.")
    parser.add_argument("--bimodal", action="store_true", help="Assume the residue distribution is bimodal")

    args = parser.parse_args()
    main(args)
