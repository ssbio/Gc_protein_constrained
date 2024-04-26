import pandas as pd
import cobra
import numpy as np
from cobra.flux_analysis import pfba
import matplotlib.pyplot as plt

# Load the COBRA model
model = cobra.io.read_sbml_model("NGO_557_irreversible_model.xml")

# Update reaction bounds based on the provided values
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
    "EX_ca2_e__f", "EX_cl_e__f", "EX_cobalt2_e__f", "EX_no_e__b",
    "EX_cu2_e__f", "EX_cys_L_e__f", "EX_cystine_L__f", "EX_fe3_e__f",
    "EX_glc_D_e__f", "EX_gln_L_e__f", "EX_glu_L_e__f", "EX_gly_e__f",
    "EX_gthrd_e__f", "EX_his_L_e__f", "EX_ade_e__b", "EX_ala_D_e__b",
    "EX_hxan_e__f", "EX_ile_L_e__f", "EX_k_e__f", "EX_leu_L_e__f",
    "EX_lys_L_e__f", "EX_met_L_e__f", "EX_mg2_e__f", "EX_mn2_e__f",
    "EX_mobd_e__f", "EX_na1_e__f", "EX_nh4_e__f", "EX_no3_e__f",
    "EX_phe_L_e__f", "EX_pi_e__f", "EX_pnto_R_e__f", "EX_gua_e__b",
    "EX_pro_L_e__f", "EX_ser_L_e__f", "EX_so4_e__f", "EX_thm_e__f",
    "EX_thr_L_e__f", "EX_trp_L_e__f", "EX_tungs_e__f", "EX_tyr_L_e__f",
    "EX_ura_e__f", "EX_val_L_e__f", "EX_zn2_e__f", "EX_apoACP_c__f",
    "EX_trdrd_c__f", "EX_HgFe3_e__b", "EX_ile_L_e__b", "EX_lac_D_e__b",
    "EX_o2_e__f", "EX_pyr_e__b", "EX_xan_e__b", "EX_lac_L_e__b", "EX_no2_e__b",
    "EX_ac_e__b", "EX_cit_e__b", "EX_akg_e__b", "URIDK2r_b","PPKr_b"
]

# Apply the bounds
for reaction_id in reactions_bound_10:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 10

for reaction_id in reactions_bound_0:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0

model.reactions.get_by_id("EX_o2_e__b").lower_bound = 0
model.reactions.get_by_id("EX_o2_e__b").upper_bound = 20

# Specifically setting EX_h_e__b bounds
model.reactions.get_by_id("EX_h_e__b").lower_bound = 0
model.reactions.get_by_id("EX_h_e__b").upper_bound = 1000
# Specifically setting EX_h_e__b bounds
model.reactions.get_by_id("EX_h_e__f").lower_bound = 0
model.reactions.get_by_id("EX_h_e__f").upper_bound = 1000

# Pre-load reaction and protein costs data
reaction_data = pd.read_excel("kcat_mw.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])
reaction_protein_costs = {row["Reactions"]: row["Protein_cost"] for _, row in reaction_data.iterrows() if row["Reactions"] in model.reactions}

# Initialize the protein cost constraint outside the loop
total_protein_cost_expression = sum(
    model.reactions.get_by_id(reaction).flux_expression * cost
    for reaction, cost in reaction_protein_costs.items()
)
protein_cost_constraint = model.problem.Constraint(total_protein_cost_expression, lb=600, ub=600)
model.add_cons_vars(protein_cost_constraint)
# Identifying ATP-producing reactions optimized
atp_producing_reactions = [rxn.id for rxn in model.reactions if any(met.id == 'atp_c[cytosol]' and coeff > 0 for met, coeff in rxn.metabolites.items())]
# Initialize lists to store results
biomass_fluxes = []
acetate_fluxes = []
glucose_fluxes_recorded = []
SUCD1_fluxes = []
SUCOAS_fluxes = []
PGK_fluxes = []
CS_fluxes = []
ACONTb_fluxes = []
all_fluxes_dict = {}

# Main loop
results = []
glucose_flux_range = np.arange(1, 5, 0.5)
for glc_flux in glucose_flux_range:
    model.reactions.EX_glc_D_e__b.bounds = (glc_flux, glc_flux)
    print(glc_flux)
    # First, remove the constraint if it exists from a previous iteration
    if 'custom_constraint' in model.constraints:
        model.remove_cons_vars('custom_constraint')
    
    atps4rpp_flux_expression = model.reactions.ATPS4rpp_f.flux_expression
    sucd1_flux_expression = model.reactions.SUCD1.flux_expression
    custom_constraint = model.problem.Constraint(atps4rpp_flux_expression - 0.65 * sucd1_flux_expression, lb=0, ub=0, name='custom_constraint')
    model.add_cons_vars(custom_constraint)
    # Step 1: Maximize Biomass Production
    model.objective = 'biomass'
    biomass_solution = model.optimize()
    max_biomass = biomass_solution.objective_value

    # Step 2: Set Biomass Production to its Maximum as a Constraint
    model.reactions.get_by_id('biomass').lower_bound = max_biomass

    # Then Optimize for Acetate Production
    model.objective = 'EX_ac_e__f'
    solution = model.optimize()
    if solution.status == 'optimal':
        results.append((glc_flux, solution.fluxes['biomass'], solution.fluxes['EX_ac_e__f'], solution.fluxes.get('SUCD1', 0), solution.fluxes.get('SUCOAS_b', 0)))
    # Record the results
    biomass_fluxes.append(solution.fluxes['biomass'])
    acetate_fluxes.append(solution.fluxes['EX_ac_e__f'])
    SUCD1_fluxes.append(solution.fluxes['SUCD1'])
    SUCOAS_fluxes.append(solution.fluxes['SUCOAS_b'])
    PGK_fluxes.append(solution.fluxes['PGK_b'])
        # Record the CS_f and ACONTb_f fluxes
    CS_fluxes.append(solution.fluxes.get('CS', 0))
    ACONTb_fluxes.append(solution.fluxes.get('ACONTb_f', 0))          
    glucose_fluxes_recorded.append(glc_flux)

    # Store fluxes for all reactions at the current glucose uptake rate
    all_fluxes_dict[glc_flux] = solution.fluxes.to_dict()    

# Plot the results
# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 10)) # 2x3 grid of subplots, with custom figure size

# Subplot 1: Glucose Uptake Rate vs. Biomass Flux
axs[0, 0].plot(glucose_fluxes_recorded, PGK_fluxes, label='PGK Flux', marker='o', color='blue')
axs[0, 0].set_title('', fontsize=11, fontweight='bold')
axs[0, 0].set_xlabel('Glucose Uptake Rate (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[0, 0].set_ylabel('PGK Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[0, 0].grid(True)

# Subplot 2: Glucose Uptake Rate vs. Acetate Flux
axs[0, 1].plot(glucose_fluxes_recorded, acetate_fluxes, label='Acetate Flux', marker='x', color='red')
axs[0, 1].set_title('')
axs[0, 1].set_xlabel('Glucose Uptake Rate (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[0, 1].set_ylabel('Acetate Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[0, 1].grid(True)

# Subplot 3: Biomass Flux vs. Acetate Flux
axs[0, 2].plot(biomass_fluxes, acetate_fluxes, label='Acetate Flux', marker='>', color='green')
axs[0, 2].set_title('')
axs[0, 2].set_xlabel('Growth Rate (/h)', fontsize=11, fontweight='bold')
axs[0, 2].set_ylabel('Acetate Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[0, 2].grid(True)

# Subplot 4: Glucose Uptake Rate vs. SUCD1 Flux
axs[1, 0].plot(glucose_fluxes_recorded, SUCD1_fluxes, label='SUCD1 Flux', marker='^', color='purple')
axs[1, 0].set_title('')
axs[1, 0].set_xlabel('Glucose Uptake Rate (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 0].set_ylabel('SUCD1 Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 0].grid(True)

axs[1, 1].plot(glucose_fluxes_recorded, CS_fluxes, label='CS Flux', marker='o', color='magenta')
axs[1, 1].set_title('')
axs[1, 1].set_xlabel('Glucose Uptake Rate (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 1].set_ylabel('CS Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 1].grid(True)

# New Subplot 6: Glucose Uptake Rate vs. ACONTb_f Flux
axs[1, 2].plot(glucose_fluxes_recorded, ACONTb_fluxes, label='ACONTb Flux', marker='o', color='cyan')
axs[1, 2].set_title('')
axs[1, 2].set_xlabel('Glucose Uptake Rate (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 2].set_ylabel('ACONTb Flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
axs[1, 2].grid(True)
# After plotting, before calling plt.tight_layout()
plt.suptitle('', fontsize=16, fontweight='bold', verticalalignment='top')

# Call plt.tight_layout() after setting the title to ensure that the layout is adjusted to accommodate the title.
plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the rect parameter as needed to fit the title
# Fine-tune figure; make subplots close to each other and hide x ticks for top plots and y ticks for right plots.
plt.tight_layout()

# Save the figure to a file
plt.savefig("glucose_flux_variation_subplots_high.png")  # Specify your desired path and filename

# Show the plots
plt.show()

print(f"Subplots saved to glucose_flux_variation_subplots.png")


# Optionally, save the results to an Excel file
results_df = pd.DataFrame({
    'Glucose Uptake Rate': glucose_fluxes_recorded,
    'Biomass Flux': biomass_fluxes,
    'Acetate Flux': acetate_fluxes,
    'SUCD1 Flux': SUCD1_fluxes,
    'PGK Flux': PGK_fluxes
})
# Convert the all_fluxes_dict to a DataFrame for easy writing to Excel
all_fluxes_df = pd.DataFrame(all_fluxes_dict).T
all_fluxes_df.index.name = 'Glucose Uptake Rate'

with pd.ExcelWriter("results_with_fluxes_high.xlsx", engine='openpyxl') as writer:
    results_df.to_excel(writer, sheet_name='Summary', index=False)
    all_fluxes_df.to_excel(writer, sheet_name='All Fluxes')
    # Directly generating ATP fluxes DataFrame and saving
    pd.DataFrame({glc_flux: {rxn: all_fluxes_dict[glc_flux].get(rxn, 0) for rxn in atp_producing_reactions} for glc_flux in glucose_flux_range}).T.to_excel(writer, sheet_name='ATP Producing Reactions Fluxes')

print("Results, fluxes, and ATP-producing reactions' fluxes saved to results_with_fluxes.xlsx")