<div align="center">
  <img src="https://github.com/user-attachments/assets/9f3e511f-5e83-4cc8-ab86-42cb6ad15838" width="500">
</div>

Optimizing these SMD-AFM curve fits is substantially harder than it looks. I'd recommend sticking to the CLI interface, even if you have solid programming knowledge. The output is just a numpy file that combines both AFM and SMD curves.

<div align="center">
  <img src="https://github.com/user-attachments/assets/1622bb35-0854-4dd2-a957-add0bbe4c886" width="500">
</div>

```
python pearsonCorrelation.py --force_extension_dir "/path/to/force_extension_files" --afm_dir "/path/to/afm_files" --N 15 --plot --save_dir "/path/to/save_dir" --show_optimal_curve_pairing
```
