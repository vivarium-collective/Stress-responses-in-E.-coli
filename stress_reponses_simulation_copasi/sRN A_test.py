from basico import *
import matplotlib.pyplot as plt

# Create a new model
new_model(name='sRNA_stochastic')

# Compartment
add_compartment('cytoplasm', 1.0)

# Species (all counts)
set_species(name='mRNA', initial_concentration=0)
set_species(name='sRNA', initial_concentration=0)
set_species(name='duplex', initial_concentration=0)
set_species(name='protein', initial_concentration=0)
set_species(name='Source', initial_concentration=1000000)  # infinite pool for transcription

# Transcription (mass action from Source)
add_reaction(name='transcription_mRNA',
             scheme='Source -> Source + mRNA',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(transcription_mRNA).k1', value=1.0)

add_reaction(name='transcription_sRNA',
             scheme='Source -> Source + sRNA',
             rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(transcription_sRNA).k1', value=0.5)

# Degradation
add_reaction(name='degradation_mRNA', scheme='mRNA -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_mRNA).k1', value=0.1)

add_reaction(name='degradation_sRNA', scheme='sRNA -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_sRNA).k1', value=0.1)

# Binding
add_reaction(name='binding', scheme='mRNA + sRNA -> duplex', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(binding).k1', value=0.01)

# Duplex degradation
add_reaction(name='degradation_duplex', scheme='duplex -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_duplex).k1', value=0.5)

# Translation
add_reaction(name='translation', scheme='mRNA -> mRNA + protein', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(translation).k1', value=2.0)

# Protein degradation
add_reaction(name='degradation_protein', scheme='protein -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_protein).k1', value=0.05)

# --- Plot stochastic simulations ---
fig, ax = plt.subplots(figsize=(8,5))

for i in range(30):  # stochastic cloud
    result = run_time_course(duration=100, step_number=200,
                             method='stochastic', use_numbers=True)
    result.plot(y='mRNA', ax=ax, color='blue', alpha=0.2, legend=None)
    result.plot(y='sRNA', ax=ax, color='red', alpha=0.2, legend=None)
    result.plot(y='protein', ax=ax, color='green', alpha=0.2, legend=None)

# --- Deterministic overlay ---
det_result = run_time_course(duration=100, step_number=200, method='deterministic')
det_result.plot(y='mRNA', ax=ax, color='blue', linewidth=2, label='mRNA (deterministic)')
det_result.plot(y='sRNA', ax=ax, color='red', linewidth=2, label='sRNA (deterministic)')
det_result.plot(y='protein', ax=ax, color='green', linewidth=2, label='protein (deterministic)')

ax.set_xlabel("Time (min)")
ax.set_ylabel("Molecule count")
ax.set_title("sRNA–mRNA Regulation\nStochastic Runs vs Deterministic Overlay")
ax.legend()
plt.show()
