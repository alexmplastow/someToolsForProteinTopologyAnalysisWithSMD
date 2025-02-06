import argparse
import functions
import objects
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main(force_extension_dir, afm_dir, smd_velocity, smd_prominence, smd_width, afm_prominence, afm_width, max_translocation_length, N, plot, save_dir, show_optimal_curve_pairing):
    force_extension_files = sorted(glob.glob(f"{force_extension_dir}/*"))
    afm_files_list = sorted(glob.glob(f"{afm_dir}/*"))
    
    afm_files = []        # Names of AFM files
    smd_files = []        # Names of SMD files
    max_rho_values = []   # Corresponding max_rho values
    
    for force_extension_file in force_extension_files:
        for i, afm_file in enumerate(afm_files_list):
            SMDforceExtension = objects.forceExtension(force_extension_file, float(smd_velocity), numpyFile=True)
            SMDforceExtension.convolve()
            
            AFMReplicate = objects.AFMData(afm_file, i)
            peakAbcissa, peakOrdinate = AFMReplicate.findProteinPeaks(prominence=afm_prominence, width=afm_width)
            highestAFMpeak = max(peakOrdinate)                       
            SMDforceExtension.shrinkOrdinate(highestAFMpeak, smd_prominence, SMDwidth=smd_width)
            
            ρ = functions.crossCorrelateSMDandAFMdata(SMDforceExtension, AFMReplicate, plot=False)
            
            SMDEDictC_i = SMDforceExtension.__dict__.copy()
            SMDCpd_i = objects.forceExtension(**SMDEDictC_i, manualInputs=True)                
            CCF = objects.crossCorrelationFit(AFMReplicate, SMDforceExtension, SMDCpd_i, ρ)
            
            P = CCF.showPearsonCorrelationPlot(save_dir, N, max_translocation_length, show=False, showOptimalCurvePairing=show_optimal_curve_pairing, returnP=True, saveOptimalSMDcurve=True)
            
            CCF.max_ρ = np.nanmax(P)
            
            afm_files.append(os.path.basename(afm_file))
            smd_files.append(os.path.basename(force_extension_file))
            max_rho_values.append(CCF.max_ρ)
    
    data = pd.DataFrame({
        'AFM_File': afm_files,
        'SMD_File': smd_files,
        'max_rho': max_rho_values
    })
    
    heatmap_data = data.pivot(index='AFM_File', columns='SMD_File', values='max_rho')
    
    if plot:
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, annot=True, cmap="YlGnBu", cbar_kws={'label': 'Values'})
        plt.title("Heatmap of SMD and AFM File Correlations")
        plt.xlabel("SMD Files")
        plt.ylabel("AFM Files")
        plt.show()
    else:
        print(heatmap_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Pearson correlation between SMD and AFM data.")
    parser.add_argument("--force_extension_dir", type=str, required=True, help="Directory containing force extension .npy files")
    parser.add_argument("--afm_dir", type=str, required=True, help="Directory containing AFM .csv files")
    parser.add_argument("--smd_velocity", type=float, default=0.000002, help="SMD velocity value")
    parser.add_argument("--smd_prominence", type=float, default=5, help="SMD prominence for peak detection")
    parser.add_argument("--smd_width", type=int, default=2, help="SMD peak width")
    parser.add_argument("--afm_prominence", type=float, default=0.0000000000025, help="AFM prominence for peak detection")
    parser.add_argument("--afm_width", type=int, default=5, help="AFM peak width")
    parser.add_argument("--max_translocation_length", type=int, default=750, help="Max translocation length")
    parser.add_argument("--N", type=int, default=10, help="Number of comparisons")
    parser.add_argument("--plot", action="store_true", help="Show the heatmap plot")
    parser.add_argument("--save_dir", type=str, default="./", help="Directory to save output figures")
    parser.add_argument("--show_optimal_curve_pairing", action="store_true", help="Show optimal curve pairing in Pearson correlation plots")
    
    args = parser.parse_args()
    
    main(
        force_extension_dir=args.force_extension_dir,
        afm_dir=args.afm_dir,
        smd_velocity=args.smd_velocity,
        smd_prominence=args.smd_prominence,
        smd_width=args.smd_width,
        afm_prominence=args.afm_prominence,
        afm_width=args.afm_width,
        max_translocation_length=args.max_translocation_length,
        N=args.N,
        plot=args.plot,
        save_dir=args.save_dir,
        show_optimal_curve_pairing=args.show_optimal_curve_pairing
    )
           
        
                  
#TODO: make your heatmap print

#TODO: write an additional routine to retool the force extension curve
        
