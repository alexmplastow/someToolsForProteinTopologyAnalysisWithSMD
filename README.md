# Steered Molecular Dynamics (SMD) for Topology Analysis

## Dependencies
- Python 3.9+
- NumPy (`pip install numpy`)
- MDAnalysis (`pip install MDAnalysis`)


TODO: add references

<h2 align="center">Detachment Modeling</h2>

<div align="center">
  <video src="https://github.com/user-attachments/assets/96f636e1-e06d-456f-8585-ca0ca591b162" controls></video>
</div>

This github provides a set of analytical routines for analyzing steered molecular dynamics (SMD) simulation. It contains CLIs and a library of objects which are useful for analyzing detachment events of secondary structures and a pure-python implementation for Dynamical Network Analysis (DNA). I'd recommend relying on the library of objects for regular python users who have an intuition for how objects work in python and the CLI interface for everyone else.
<div align="center">
  <video src="https://github.com/user-attachments/assets/ba0f12cc-8e1c-47e9-8cea-70e6ae93aa1b" controls></video>
</div>

The thinking for modeling detachments is analogous to statistical inference, anytime a set of residues associated with a secondary structure is about μ +/- 2σ. This parameterization is arbitrary, so I encourage tinkering with the "p_value" parameter. The default critical values are displayed as the black lines in this visualization. The assumption of normally distributed residues along a certain axis breaks down as a protein denatures. 

<div align="center">
  <img src="https://github.com/user-attachments/assets/42d59615-e2e6-4c26-a563-911cf5512be0" width="500">
</div>

If a residue acquires a p_value <= 0.05, it may be considered detached

<h2 align="center">Dynamical Network Analysis (DNA)</h2>

DNA is handy for determining the path along which force most easily propagates.

<div align="center">
  <img src="https://github.com/user-attachments/assets/617e7802-0984-42a4-aa8a-56438de133af" width="200">
</div>

The force will easily propagate if the cross-correlation between two particles is small. A less mathematical description is that those portions of a protein that are highly correlated with each other are being pulled about by their neighbors, but those particles that are weakly correlated with each other are the path along which force is propagating and determine the dynamics of the macromolecule.

That's why DNA is useful in SMD.

<div align="center">
  <video src="https://github.com/user-attachments/assets/b5b90ecd-8686-4bb7-8d27-652fe73412cb" controls></video>
</div>

<div align="center">
  <img src="https://github.com/user-attachments/assets/6c8867dd-16a1-43dc-80be-ed2661525704" width="300">
</div>

ρ is the pearson correlation of the gradient for any two particles. Particle pairs with a high ρ are more likely to moving as a consequence of eachother than particle pairs with a low ρ. To rule out correlations which result from sheer dumb luck, contacts can be identified.

<div align="center">
  <img src="https://github.com/user-attachments/assets/6c56fabc-1a2b-4f81-846f-6381da1efe0c" width="900">
</div>

<div align="center">
  <img src="https://github.com/user-attachments/assets/a3fcbd41-d9a9-4bd4-8b87-33f8b282be3c" width="300">
</div>

In the VMD implementation, contacts are defined as particles that are within 4.5 Å of each other for 75% of the trajectory which is analyzed. Matrix V is defined as 1 if there is a contact and 0 if any two residue are not contacting.

<div align="center">
  <img src="https://github.com/user-attachments/assets/64faa79e-5e98-40bc-bb90-87cd9eb70b5b" width="600">
</div>

The resulting network, G, can be visualized in VMD.

<div align="center">
  <img src="https://github.com/user-attachments/assets/8574cfa2-e13f-4575-9f34-0c23d8fca0b3" width="600">
</div>

As can optimal and subpotimal paths, this makes much of networkx's analytical reliability compatable with VMD, particularly in cases where the clusters become bimodal.

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


