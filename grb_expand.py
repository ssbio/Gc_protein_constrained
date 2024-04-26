import cobra
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import matplotlib.pyplot as plt


def read_excel_data(xlsx_path):
    xls = pd.ExcelFile(xlsx_path)
    kcat_df = pd.read_excel(xls, 'Kcat', index_col=0)
    mw_df = pd.read_excel(xls, 'MW', index_col=0)
    return kcat_df, mw_df


def update_reaction_bounds(model, reactions_bound_10, reactions_bound_0):
    """Updates reaction bounds in the model based on specified lists."""
    for rxn_id in reactions_bound_10:
        model.reactions.get_by_id(rxn_id).bounds = (0, 10)
    for rxn_id in reactions_bound_0:
        model.reactions.get_by_id(rxn_id).bounds = (0, 0)
    model.reactions.get_by_id("EX_o2_e__b").bounds = (0, 20)
    model.reactions.get_by_id("EX_h_e__b").bounds = (0, 1000)
    model.reactions.get_by_id("biomass").bounds = (2.26261, 2.26261)

def optimize_phi(model, pi_values):
    m = gp.Model("minimize_phi")
    fluxes = {rxn.id: m.addVar(lb=0, ub=1000, name=rxn.id) for rxn in model.reactions}
    m.update()
    
    # Add mass balance constraints for each metabolite
    for met in model.metabolites:
        m.addConstr(
            gp.quicksum(fluxes[rxn.id] * rxn.metabolites[met] for rxn in met.reactions if rxn.id in fluxes) == 0, 
            name=f"mass_balance_{met.id}"
        )
    # Add a constraint that forces non-zero flux through an essential reaction
    m.addConstr(fluxes["biomass"] >= 2.26261, "demand_biomass")

    # Define the objective to minimize Φ with a small weight on the sum of fluxes to minimize them as a secondary objective
    primary_objective = gp.quicksum((fluxes[rxn_id] / 1000) * pi for rxn_id, pi in pi_values.items())
    secondary_objective = gp.quicksum((fluxes[rxn_id] / 1000) for rxn_id in fluxes) * 0.001  # Adjust the weight as necessary
    m.setObjective(primary_objective + secondary_objective, GRB.MINIMIZE)

    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        return m.objVal, {rxn: fluxes[rxn].x for rxn in fluxes}
    else:
        return float('inf'), {}

# Load the model
model = cobra.io.read_sbml_model("NGO_557_irreversible_model.xml")

# Update model bounds
reactions_bound_10 = [

    "EX_ala_L_e__b", "EX_arg_L_e__b", "EX_asn_L_e__b", "EX_asp_L_e__b",
    "EX_ca2_e__b", "EX_cl_e__b", "EX_co2_e__b", "EX_cobalt2_e__b",
    "EX_cu2_e__b", "EX_cys_L_e__b", "EX_cystine_L__b", "EX_fe3_e__b",
    "EX_glc_D_e__b", "EX_gln_L_e__b", "EX_glu_L_e__b", "EX_gly_e__b",
    "EX_gthrd_e__b", "EX_h_e__b", "EX_h2o_e__b", "EX_his_L_e__b",
    "EX_hxan_e__b", "EX_ile_L_e__b", "EX_k_e__b", "EX_leu_L_e__b",
    "EX_lys_L_e__b", "EX_met_L_e__b", "EX_mg2_e__b", "EX_mn2_e__b",
    "EX_mobd_e__b", "EX_na1_e__b", "EX_nh4_e__b", "EX_no3_e__b",
    "EX_phe_L_e__b", "EX_pi_e__b", "EX_pnto_R_e__b",
    "EX_pro_L_e__b", "EX_ser_L_e__b", "EX_so4_e__b", "EX_thm_e__b",
    "EX_thr_L_e__b", "EX_trp_L_e__b", "EX_tungs_e__b", "EX_tyr_L_e__b",
    "EX_ura_e__b", "EX_val_L_e__b", "EX_zn2_e__b", "EX_apoACP_c__b",
    "EX_trdrd_c__b",
]
# Subset of reactions with bounds set to 0
reactions_bound_0 = [

    "EX_ala_L_e__f", "EX_arg_L_e__f", "EX_asn_L_e__f", "EX_asp_L_e__f",
    "EX_ca2_e__f", "EX_cl_e__f", "EX_co2_e__f", "EX_cobalt2_e__f",
    "EX_cu2_e__f", "EX_cys_L_e__f", "EX_cystine_L__f", "EX_fe3_e__f",
    "EX_glc_D_e__f", "EX_gln_L_e__f", "EX_glu_L_e__f", "EX_gly_e__f",
    "EX_gthrd_e__f", "EX_h_e__f", "EX_h2o_e__f", "EX_his_L_e__f",
    "EX_hxan_e__f", "EX_ile_L_e__f", "EX_k_e__f", "EX_leu_L_e__f",
    "EX_lys_L_e__f", "EX_met_L_e__f", "EX_mg2_e__f", "EX_mn2_e__f",
    "EX_mobd_e__f", "EX_na1_e__f", "EX_nh4_e__f", "EX_no3_e__f",
    "EX_phe_L_e__f", "EX_pi_e__f", "EX_pnto_R_e__f",
    "EX_pro_L_e__f", "EX_ser_L_e__f", "EX_so4_e__f", "EX_thm_e__f",
    "EX_thr_L_e__f", "EX_trp_L_e__f", "EX_tungs_e__f", "EX_tyr_L_e__f",
    "EX_ura_e__f", "EX_val_L_e__f", "EX_zn2_e__f", "EX_apoACP_c__f",
    "EX_trdrd_c__f"
]
update_reaction_bounds(model, reactions_bound_10, reactions_bound_0)


kcat_df, mw_df = read_excel_data("output.xlsx")  # Update this path

phi_results = []
all_fluxes = {}
min_phi = float('inf')
optimal_flux_values = None
optimal_simulation = None

# Perform optimizations for each simulation set
for sim_num in range(1, 101):
    sim_str = f'Simulation {sim_num}:'  # Assuming the colon is part of the header based on the error
    if sim_str in kcat_df.columns and sim_str in mw_df.columns:
        sim_kcat = kcat_df[sim_str].dropna()
        sim_mw = mw_df[sim_str].dropna()
        pi_values = {rxn: sim_mw[rxn] / sim_kcat[rxn] for rxn in sim_kcat.index.intersection(sim_mw.index)}
        
        phi, fluxes = optimize_phi(model, pi_values)

    
    phi_results.append(phi)
    all_fluxes[sim_num] = fluxes
    
    if phi < min_phi:
        min_phi = phi
        optimal_flux_values = fluxes
        optimal_simulation = sim_num
# Write the results to an Excel file
with pd.ExcelWriter('optimization_results.xlsx', engine='openpyxl') as writer:  # Update this path
    pd.DataFrame({'Simulation': range(1, 101), 'Phi Value': phi_results}).to_excel(writer, sheet_name='Phi Values')
    if optimal_flux_values:
        pd.DataFrame({'Reaction ID': list(optimal_flux_values.keys()), 'Flux': list(optimal_flux_values.values())}).to_excel(writer, sheet_name=f'Optimal Fluxes Sim {optimal_simulation}')

print("Optimization completed. Results saved.")

# Now, create a histogram of the Phi values
plt.figure(figsize=(10, 6))
plt.hist(phi_results, bins=20, color='skyblue', edgecolor='black')
plt.title('Distribution of Phi Values Across 100 Simulations')
plt.xlabel('Phi Value')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# Save the histogram to a file
plt.savefig('phi_values_histogram.png', dpi=300)

# Optionally, show the histogram in a window (this line can be omitted if running in a non-interactive environment)
plt.show()

print("Optimization completed. Results and histogram saved.")