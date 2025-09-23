from basico import *
import matplotlib.pyplot as plt

# --- New Model ---
new_model(name='Ecoli_ROS_Response')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (molecule counts, SBML IDs)
# Oxidative stress (H2O2 / OxyR)
set_species(name='H2O2', sbml_id='H2O2', initial_concentration=50)
set_species(name='OxyR', sbml_id='OxyR', initial_concentration=100)
set_species(name='OxyR*', sbml_id='OxyRox', initial_concentration=0)
set_species(name='KatG', sbml_id='KatG', initial_concentration=0)
set_species(name='AhpCF', sbml_id='AhpCF', initial_concentration=0)

# Superoxide stress (O2- / SoxRS)
set_species(name='O2-', sbml_id='O2m', initial_concentration=30)   # superoxide anion
set_species(name='SoxR', sbml_id='SoxR', initial_concentration=50) # inactive
set_species(name='SoxR*', sbml_id='SoxRox', initial_concentration=0) # oxidized active
set_species(name='SoxS', sbml_id='SoxS', initial_concentration=0)   # transcriptional activator
set_species(name='SodA', sbml_id='SodA', initial_concentration=0)   # superoxide dismutase

# --- Reactions ---
## OxyR system
add_reaction(name='OxyR_activation', sbml_id='R_act',
             scheme='OxyR + H2O2 -> OxyRox', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_act).k1', value=0.01)

add_reaction(name='KatG_synthesis', sbml_id='R_katG',
             scheme='OxyRox -> OxyRox + KatG', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_katG).k1', value=0.05)

add_reaction(name='AhpCF_synthesis', sbml_id='R_ahpCF',
             scheme='OxyRox -> OxyRox + AhpCF', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_ahpCF).k1', value=0.02)

add_reaction(name='KatG_detox', sbml_id='R_detox1',
             scheme='H2O2 + KatG -> KatG', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_detox1).k1', value=0.01)

add_reaction(name='AhpCF_detox', sbml_id='R_detox2',
             scheme='H2O2 + AhpCF -> AhpCF', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_detox2).k1', value=0.005)

add_reaction(name='KatG_degradation', sbml_id='R_dKatG',
             scheme='KatG -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dKatG).k1', value=0.001)

add_reaction(name='AhpCF_degradation', sbml_id='R_dAhp',
             scheme='AhpCF -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dAhp).k1', value=0.001)

## SoxRS system
add_reaction(name='SoxR_activation', sbml_id='R_soxR',
             scheme='SoxR + O2m -> SoxRox', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_soxR).k1', value=0.02)

add_reaction(name='SoxS_synthesis', sbml_id='R_soxS',
             scheme='SoxRox -> SoxRox + SoxS', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_soxS).k1', value=0.05)

add_reaction(name='SodA_synthesis', sbml_id='R_sodA',
             scheme='SoxS -> SoxS + SodA', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_sodA).k1', value=0.03)

add_reaction(name='SodA_detox', sbml_id='R_sodA_detox',
             scheme='O2m + SodA -> SodA', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_sodA_detox).k1', value=0.02)

add_reaction(name='SoxS_degradation', sbml_id='R_dSoxS',
             scheme='SoxS -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dSoxS).k1', value=0.005)

add_reaction(name='SodA_degradation', sbml_id='R_dSodA',
             scheme='SodA -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dSodA).k1', value=0.001)

# --- Simulation ---
result = run_time_course(duration=200, step_number=400, method='deterministic')

# --- Plot ---
fig, ax = plt.subplots(figsize=(10,6))
result.plot(y='H2O2', ax=ax, color='black', label='H2O2')
result.plot(y='O2m', ax=ax, color='orange', label='O2- (superoxide)')
result.plot(y='OxyRox', ax=ax, color='red', label='OxyR*')
result.plot(y='KatG', ax=ax, color='blue', label='KatG')
result.plot(y='AhpCF', ax=ax, color='green', label='AhpCF')
result.plot(y='SoxS', ax=ax, color='purple', label='SoxS')
result.plot(y='SodA', ax=ax, color='brown', label='SodA')

ax.set_xlabel("Time (a.u.)")
ax.set_ylabel("Molecule count")
ax.set_title("E. coli ROS Stress Response: OxyR (H2O2) + SoxRS (O2-)")
ax.legend()
plt.show()

# --- Export ---
# save_model('Ecoli_ROS_model.cps')
# save_model('Ecoli_ROS_model.xml', sbml=True)
# print("Model exported as Ecoli_ROS_model.cps and Ecoli_ROS_model.xml")

fig, axes = plt.subplots(2, 2, figsize=(10,8), sharex=True)

result.plot(y='H2O2', ax=axes[0,0], color='black', legend=None)
axes[0,0].set_title("Hydrogen Peroxide")

result.plot(y='O2m', ax=axes[0,1], color='orange', legend=None)
axes[0,1].set_title("Superoxide")

result.plot(y=['KatG','AhpCF'], ax=axes[1,0])
axes[1,0].set_title("OxyR regulon (peroxide defense)")

result.plot(y=['SoxS','SodA'], ax=axes[1,1])
axes[1,1].set_title("SoxRS regulon (superoxide defense)")

# plt.tight_layout()
# plt.show()

print(result)
