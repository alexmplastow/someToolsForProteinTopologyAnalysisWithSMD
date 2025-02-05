# Steered Molecular Dynamics (SMD) for Topology Analysis

This github provides a set of analytical routines for analyzing steered molecular dynamics (SMD) simulation. It contains CLIs and a library of objects which are useful for analyzing detachment events of secondary structures and a pure-python implementation for Dynamical Network Analysis (DNA). I'd recommend relying on the library of objects for regular python users who have an intuition for how objects work in python and the CLI interface for everyone else.

## Dependencies
- Python 3.12+
  - NumPy
  - Pandas
  - Matplotlib
  - Scipy
  - MDAnalysis
  - Biopython
  - Glob
  - Tqdm
- VMD

<h2 align="center">Detachment Modeling</h2>

<div align="center">
  <video src="https://github.com/user-attachments/assets/96f636e1-e06d-456f-8585-ca0ca591b162" controls></video>
</div>

```
cd VMDroutines
vmd -e ssFinderResID.tcl /path/to/*pdb /path/to/*psf /path/to/*dcd /path/to/secondaryStructure_resID.txt
```
or
```
vmd -e ssFinderResidue.tcl /path/to/*pdb /path/to/*psf /path/to/*dcd /path/to/secondaryStructure_residue.txt
```

<div align="center">
  <video src="https://github.com/user-attachments/assets/ba0f12cc-8e1c-47e9-8cea-70e6ae93aa1b" controls></video>
</div>

The thinking for modeling detachments is analogous to statistical inference, anytime a set of residues associated with a secondary structure is about μ +/- 2σ, it is detached. This parameterization is arbitrary, so I encourage tinkering with the "p_value" parameter. The default critical values are displayed as black lines in this visualization. The assumption of normally distributed residues along a certain axis breaks down as a protein denatures. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/42d59615-e2e6-4c26-a563-911cf5512be0" width="500">
</div>

If a residue acquires a p_value <= 0.05, it may be considered detached

To reproduce the visualization shown above:

Access the secondary structures,

```
frameStructure_0 = objects.frameStructure(path/to/initial*pdb, ssDataFileResnum, ssDataFileResidue)
sStructures = frameStructure_0.getSecondaryStructures().copy()
```
Then loop through using either biopython or MDAnalysis generated PDBs

```
for pdb in listOfPDBfilePaths:
    frameStructure = objects.frameStructure(pdb, ssDataFileResnum, ssDataFileResidue)
    frameStructure.frame_SSr_DebugRoutine(i = i, p_value = p_value)
```

This method can be used to identify detachments in SMD curves, if fitted to an atomic force microscopy (AFM) curve, this infers the detachment events which take place in experimental data. 

<h2 align="center">Dynamical Network Analysis (DNA)</h2>

DNA is handy for determining the path along which force most easily propagates.

<div align="center">
  <img src="https://github.com/user-attachments/assets/617e7802-0984-42a4-aa8a-56438de133af" width="200">
  <br>
  <small>[1]</small>
</div>

The force will easily propagate if the cross-correlation between two particles is small. A less mathematical description is that those portions of a protein that are highly correlated with each other are being pulled about by their neighbors, but those particles that are weakly correlated with each other are the path along which force is propagating and determine the dynamics of the macromolecule.

That's why DNA is useful in SMD.

<div align="center">
  <video src="https://github.com/user-attachments/assets/b5b90ecd-8686-4bb7-8d27-652fe73412cb" controls></video>
</div>

<div align="center">
  <img src="https://github.com/user-attachments/assets/6c8867dd-16a1-43dc-80be-ed2661525704" width="200">
  <br>
  <small>[2]</small>
</div>

ρ is the Pearson correlation of the gradient for any two particles. Particle pairs with a high ρ are more likely to move as a consequence of each other than particle pairs with a low ρ. To rule out correlations that result from sheer dumb luck, contacts can be identified.

<div align="center">
  <img src="https://github.com/user-attachments/assets/6c56fabc-1a2b-4f81-846f-6381da1efe0c" width="900">
</div>

<div align="center">
  <img src="https://github.com/user-attachments/assets/a3fcbd41-d9a9-4bd4-8b87-33f8b282be3c" width="300">
</div>

In the VMD implementation, contacts are defined as particles that are within 4.5 Å of each other for 75% of the trajectory which is analyzed [3]. Matrix V is defined as 1 if there is a contact and 0 if any two residue are not contacting.

<div align="center">
  <img src="https://github.com/user-attachments/assets/64faa79e-5e98-40bc-bb90-87cd9eb70b5b" width="600">
</div>

The resulting network, G, can be visualized in VMD.

<div align="center">
  <img src="https://github.com/user-attachments/assets/8574cfa2-e13f-4575-9f34-0c23d8fca0b3" width="600">
</div>

As can optimal and subpotimal paths, this makes much of networkx's analytical reliability compatable with VMD, particularly in cases where the clusters become bimodal.

## Equation References
[1] C. Bernardi, R. Bernardi and K. Schulten, "Mapping mechanical force propagation through biomolecular complexes," Nano Letters, vol. 15, no. 11, pp. 7370-7376, 2015. 

[2] J. Berman, "Taming Information with Open Source Tools," in Data Simplification, ScienceDirect, 2016, pp. 135-187.

[3] J. Eargle, Z. Luthey-Schulten and A. Sethi, "Dynamical Network Analysis," NIH Resource for Macromolecular Modeling and Bioinformatics, Urbana, Illinois, 2012.




<!--
<p align="center">
  <img src="https://github.com/user-attachments/assets/6c56fabc-1a2b-4f81-846f-6381da1efe0c" width="600">
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/64faa79e-5e98-40bc-bb90-87cd9eb70b5b" width="600">
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/8574cfa2-e13f-4575-9f34-0c23d8fca0b3" width="600">
</p>
-->


