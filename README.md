# Gc_protein_constrained
Optimal Protein Allocation Controls the Inhibition of GltA and AcnB in Neisseria gonorrhoeae
mw.txt and kcat.txt are the input files to tun Monte Carlo simulation (MCS_Kcat_MW.py)
To run phi_low.py/phi_high.py, the input files are NGO_557_irreversible_model.xml and kcat_mw.xlsx
To determine which simulation to pick for further analysis run minimum_phi.py (input file Kcat_MW_100simulation_input.xlsx) and take the lowest phi value simulation results. After that, incorporate the Kcat and MW values of this simulation result to build the final model.
