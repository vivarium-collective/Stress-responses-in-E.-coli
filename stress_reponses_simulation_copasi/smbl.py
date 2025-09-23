import sys
from basico import *
import numpy as np
import matplotlib.pyplot as plt

from basico import *
import matplotlib.pyplot as plt


import matplotlib.pyplot as plt
import pandas as pd

# --- Create new model ---
new_model(name='sRNA_model')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (use counts for stochastic compatibility)
set_species(name='mRNA', sbml_id='M', initial_concentration=0)
set_species(name='sRNA', sbml_id='S', initial_concentration=0)
set_species(name='duplex', sbml_id='D', initial_concentration=0)
set_species(name='protein', sbml_id='P', initial_concentration=0)
set_species(name='Source', sbml_id='SRC', initial_concentration=1000000)  # infinite pool

# --- Reactions (Mass Action only) ---
# Transcription of mRNA
add_reaction(name='transcription_mRNA', sbml_id='R_mRNA',
             scheme='SRC -> SRC + M', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_mRNA).k1', value=0.01)

# Transcription of sRNA
add_reaction(name='transcription_sRNA', sbml_id='R_sRNA',
             scheme='SRC -> SRC + S', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_sRNA).k1', value=0.005)

# Degradation of mRNA
add_reaction(name='degradation_mRNA', sbml_id='R_dm',
             scheme='M -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dm).k1', value=0.1)

# Degradation of sRNA
add_reaction(name='degradation_sRNA', sbml_id='R_ds',
             scheme='S -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_ds).k1', value=0.1)

# Duplex formation
add_reaction(name='binding', sbml_id='R_bind',
             scheme='M + S -> D', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_bind).k1', value=0.01)

# Duplex degradation
add_reaction(name='degradation_duplex', sbml_id='R_dd',
             scheme='D -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dd).k1', value=0.5)
from basico import *
import matplotlib.pyplot as plt
import pandas as pd

# --- Create new model ---
new_model(name='sRNA_model')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (molecule counts for stochastic)
set_species(name='mRNA', sbml_id='M', initial_concentration=0)
set_species(name='sRNA', sbml_id='S', initial_concentration=0)
set_species(name='duplex', sbml_id='D', initial_concentration=0)
set_species(name='protein', sbml_id='P', initial_concentration=0)
set_species(name='Source', sbml_id='SRC', initial_concentration=1000000)  # large pool

# --- Reactions (Mass Action only) ---
# Transcription
add_reaction(name='transcription_mRNA', sbml_id='R_mRNA',
             scheme='SRC -> SRC + M', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_mRNA).k1', value=0.01)

add_reaction(name='transcription_sRNA', sbml_id='R_sRNA',
             scheme='SRC -> SRC + S', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_sRNA).k1', value=0.005)

# Degradation
add_reaction(name='degradation_mRNA', sbml_id='R_dm',
             scheme='M -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dm).k1', value=0.1)

add_reaction(name='degradation_sRNA', sbml_id='R_ds',
             scheme='S -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_ds).k1', value=0.1)

# Binding
add_reaction(name='binding', sbml_id='R_bind',
             scheme='M + S -> D', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_bind).k1', value=0.01)

# Duplex degradation
add_reaction(name='degradation_duplex', sbml_id='R_dd',
             scheme='D -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dd).k1', value=0.5)

# Translation
add_reaction(name='translation', sbml_id='R_trans',
             scheme='M -> M + P', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_trans).k1', value=2.0)

# Protein degradation
add_reaction(name='degradation_protein', sbml_id='R_dp',
             scheme='P -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(R_dp).k1', value=0.05)

# --- Run simulations ---
fig, ax = plt.subplots(figsize=(8,5))

trajectories = []

# 20 stochastic runs (cloud)
for i in range(20):
    result = run_time_course(duration=100, step_number=200,
                             method='stochastic', use_numbers=True)
    trajectories.append(result[['M','S','P']])
    result.plot(y='M', ax=ax, color='blue', alpha=0.2, legend=None)
    result.plot(y='S', ax=ax, color='red', alpha=0.2, legend=None)
    result.plot(y='P', ax=ax, color='green', alpha=0.2, legend=None)

# Compute average trajectory
all_traj = pd.concat(trajectories, axis=0, keys=range(len(trajectories)))
avg_traj = all_traj.groupby(level=1).mean()

# Plot averages (dashed)
ax.plot(avg_traj.index, avg_traj['M'], color='blue', linestyle='--', linewidth=2, label='mRNA (stoch avg)')
ax.plot(avg_traj.index, avg_traj['S'], color='red', linestyle='--', linewidth=2, label='sRNA (stoch avg)')
ax.plot(avg_traj.index, avg_traj['P'], color='green', linestyle='--', linewidth=2, label='protein (stoch avg)')

# Deterministic overlay (solid)
det_result = run_time_course(duration=100, step_number=200, method='deterministic')
det_result.plot(y='M', ax=ax, color='blue', linewidth=2, label='mRNA (det)')
det_result.plot(y='S', ax=ax, color='red', linewidth=2, label='sRNA (det)')
det_result.plot(y='P', ax=ax, color='green', linewidth=2, label='protein (det)')

# Labels
ax.set_xlabel("Time (min)")
ax.set_ylabel("Molecule count")
ax.set_title("sRNA–mRNA Regulation:\nStochastic Cloud + Mean + Deterministic Overlay")
ax.legend()
plt.show()

# --- Save model ---
save_model('srna_model_fixed.cps')          # COPASI format
save_model('srna_model_fixed.xml', sbml=True)  # SBML format
print("Model exported as srna_model_fixed.cps and srna_model_fixed.xml")
