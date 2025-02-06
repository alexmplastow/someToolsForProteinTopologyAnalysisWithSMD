<h2 align="center">Fitting SMD data with AFM data</h2>

<div align="center">
  <img src="https://github.com/user-attachments/assets/9f3e511f-5e83-4cc8-ab86-42cb6ad15838" width="500">
</div>

Optimizing these SMD-AFM curve fits is substantially harder than it looks. I'd recommend sticking to the CLI interface, even if you have solid programming knowledge. The output is just a numpy file that combines both AFM and SMD curves.

```
python pearsonCorrelation.py --force_extension_dir "/path/to/force_extension_files" --afm_dir "/path/to/afm_files" --N 15 --plot --save_dir "/path/to/save_dir" --show_optimal_curve_pairing
```

<div align="center">
  <video src="https://github.com/user-attachments/assets/83b8e621-1bd8-4425-82af-8109169472f0" controls></video>
</div>


This routine will sample the best fit for the SMD curve on the AFM curve.

<div align="center">
  <img src="https://github.com/user-attachments/assets/1622bb35-0854-4dd2-a957-add0bbe4c886" width="500">
</div>

For each AFM dataset and each SMD dataset, some numpy file will be produced, but you will need to infer which is the best pairing. In this case, my fourth AFM replicate matches well with my third SMD replicate. After the script finishes, it will generate the numpy file containing the best fit for all pairings displayed in the heat map.

#TODO:
I'll add additional information for file format soon





