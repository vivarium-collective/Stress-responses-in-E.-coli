
from basico import *
import matplotlib.pyplot as plt

# --- New Model ---
new_model(name='OxyR_oxidative_stress')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (molecule counts, SBML-ready IDs)
set_species(name='H2O2', sbml_id='H2O2', initial_concentration=50)   # oxidative stress input
set_species(name='OxyR', sbml_id='OxyR', initial_concentration=100)  # inactive regulator
set_species(name='OxyR*', sbml_id='OxyRox', initial_concentration=0) # oxidized active regulator
set_species(name='KatG', sbml_id='KatG', initial_concentration=0)
set_species(name='AhpCF', sbml_id='AhpCF', initial_concentration=0)

# --- Reactions ---
# 1. OxyR activation by H2O2
add_reaction(name='OxyR_activation', sbml_id='R_act',
             scheme='OxyR + H2O2 -> OxyRox', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_act).k1', value=0.01)

# 2. KatG synthesis (activated by OxyRox)
add_reaction(name='KatG_synthesis', sbml_id='R_katG',
             scheme='OxyRox -> OxyRox + KatG', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_katG).k1', value=0.05)

# 3. AhpCF synthesis (activated by OxyRox)
add_reaction(name='AhpCF_synthesis', sbml_id='R_ahpCF',
             scheme='OxyRox -> OxyRox + AhpCF', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_ahpCF).k1', value=0.02)

# 4. H2O2 detox by KatG
add_reaction(name='KatG_detox', sbml_id='R_detox1',
             scheme='H2O2 + KatG -> KatG', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_detox1).k1', value=0.01)

# 5. H2O2 detox by AhpCF
add_reaction(name='AhpCF_detox', sbml_id='R_detox2',
             scheme='H2O2 + AhpCF -> AhpCF', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_detox2).k1', value=0.005)

# 6. Protein degradation
add_reaction(name='KatG_degradation', sbml_id='R_dKatG',
             scheme='KatG -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dKatG).k1', value=0.001)

add_reaction(name='AhpCF_degradation', sbml_id='R_dAhp',
             scheme='AhpCF -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dAhp).k1', value=0.001)

# --- Simulation ---
result = run_time_course(duration=200, step_number=400, method='deterministic')

# --- Plot ---
fig, ax = plt.subplots(figsize=(8,5))
result.plot(y='H2O2', ax=ax, color='black', label='H2O2')
result.plot(y='KatG', ax=ax, color='blue', label='KatG')
result.plot(y='AhpCF', ax=ax, color='green', label='AhpCF')
result.plot(y='OxyRox', ax=ax, color='red', label='OxyR* (active)')
ax.set_xlabel("Time (a.u.)")
ax.set_ylabel("Molecule count")
ax.set_title("E. coli OxyR Oxidative Stress Response (Minimal Model)")
ax.legend()
plt.show()

# --- Export for COPASI & SBML ---
save_model('../OxyR_stress_model.cps')          # COPASI format
save_model('../OxyR_stress_model.xml', sbml=True)  # SBML format
