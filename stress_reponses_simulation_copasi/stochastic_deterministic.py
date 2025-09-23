from basico import *
import matplotlib.pyplot as plt

# Create a new model
new_model(name='sRNA_model')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (all in molecule counts)
set_species(name='mRNA', initial_concentration=0)
set_species(name='sRNA', initial_concentration=0)
set_species(name='duplex', initial_concentration=0)
set_species(name='protein', initial_concentration=0)
set_species(name='Source', initial_concentration=1000000)  # infinite pool for transcription

# --- Reactions ---

# Transcription of mRNA (Source -> Source + mRNA)
add_reaction(name='transcription_mRNA',
             scheme='Source -> Source + mRNA',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(transcription_mRNA).k1', value=0.01)

# Transcription of sRNA (Source -> Source + sRNA)
add_reaction(name='transcription_sRNA',
             scheme='Source -> Source + sRNA',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(transcription_sRNA).k1', value=0.005)

# Degradation of mRNA
add_reaction(name='degradation_mRNA',
             scheme='mRNA -> ',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_mRNA).k1', value=0.1)

# Degradation of sRNA
add_reaction(name='degradation_sRNA',
             scheme='sRNA -> ',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_sRNA).k1', value=0.1)

# Duplex formation
add_reaction(name='binding',
             scheme='mRNA + sRNA -> duplex',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(binding).k1', value=0.01)

# Duplex degradation
add_reaction(name='degradation_duplex',
             scheme='duplex -> ',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_duplex).k1', value=0.5)

# Translation (mRNA -> mRNA + protein)
add_reaction(name='translation',
             scheme='mRNA -> mRNA + protein',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(translation).k1', value=2.0)

# Protein degradation
add_reaction(name='degradation_protein',
             scheme='protein -> ',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_protein).k1', value=0.05)

# --- Simulation ---
fig, ax = plt.subplots(figsize=(8,5))

# 20 stochastic runs (cloud)
for i in range(20):
    result = run_time_course(duration=100, step_number=200,
                             method='stochastic', use_numbers=True)
    result.plot(y='mRNA', ax=ax, color='blue', alpha=0.2, legend=None)
    result.plot(y='sRNA', ax=ax, color='red', alpha=0.2, legend=None)
    result.plot(y='protein', ax=ax, color='green', alpha=0.2, legend=None)

# Deterministic overlay
det_result = run_time_course(duration=100, step_number=200, method='deterministic')
det_result.plot(y='mRNA', ax=ax, color='blue', linewidth=2, label='mRNA (det)')
det_result.plot(y='sRNA', ax=ax, color='red', linewidth=2, label='sRNA (det)')
det_result.plot(y='protein', ax=ax, color='green', linewidth=2, label='protein (det)')

# # Axis labels
# ax.set_xlabel("Time (min)")
# ax.set_ylabel("Molecule count")
# ax.set_title("sRNA–mRNA Regulation:\nStochastic Cloud + Deterministic Overlay")
# ax.legend()
# plt.show()


fig, ax = plt.subplots(figsize=(8,5))

# 20 stochastic runs (cloud) — no labels
for i in range(20):
    result = run_time_course(duration=100, step_number=200,
                             method='stochastic', use_numbers=True)
    result.plot(y='mRNA', ax=ax, color='blue', alpha=0.2, legend=None)
    result.plot(y='sRNA', ax=ax, color='red', alpha=0.2, legend=None)
    result.plot(y='protein', ax=ax, color='green', alpha=0.2, legend=None)

# Deterministic overlay — add labels once
det_result = run_time_course(duration=100, step_number=200, method='deterministic')
det_result.plot(y='mRNA', ax=ax, color='blue', linewidth=2, label='mRNA (det)')
det_result.plot(y='sRNA', ax=ax, color='red', linewidth=2, label='sRNA (det)')
det_result.plot(y='protein', ax=ax, color='green', linewidth=2, label='protein (det)')

# Format axes
ax.set_xlabel("Time (min)")
ax.set_ylabel("Molecule count")
ax.set_title("sRNA–mRNA Regulation:\nStochastic Cloud + Deterministic Overlay")
ax.legend()
# plt.show()

# print(det_result)


