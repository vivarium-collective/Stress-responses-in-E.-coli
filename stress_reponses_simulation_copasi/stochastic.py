from basico import *
import matplotlib.pyplot as plt

# Create a new model
new_model(name='sRNA_stochastic')

# Add compartment (volume = 1, so conc ~ counts)
add_compartment('cytoplasm', 1.0)

# Add species with integer counts
set_species(name='mRNA', initial_concentration=0)
set_species(name='sRNA', initial_concentration=0)
set_species(name='duplex', initial_concentration=0)
set_species(name='protein', initial_concentration=0)

# Add reactions

# 1. Transcription mRNA
add_reaction(name='transcription_mRNA', scheme=' -> mRNA', rate_law='Constant')
set_reaction_parameters(name='(transcription_mRNA).Flux', value=1)

# 2. Transcription sRNA
add_reaction(name='transcription_sRNA', scheme=' -> sRNA', rate_law='Constant')
set_reaction_parameters(name='(transcription_sRNA).Flux', value=1)

# 3. Degradation mRNA
add_reaction(name='degradation_mRNA', scheme='mRNA -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_mRNA).k1', value=0.1)

# 4. Degradation sRNA
add_reaction(name='degradation_sRNA', scheme='sRNA -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_sRNA).k1', value=0.1)

# 5. Duplex formation
add_reaction(name='binding', scheme='mRNA + sRNA -> duplex', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(binding).k1', value=0.01)

# 6. Duplex degradation
add_reaction(name='degradation_duplex', scheme='duplex -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_duplex).k1', value=0.5)

# 7. Translation
add_reaction(name='translation', scheme='mRNA -> mRNA + protein', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(translation).k1', value=2.0)

# 8. Protein degradation
add_reaction(name='degradation_protein', scheme='protein -> ', rate_law='Mass action (irreversible)')
set_reaction_parameters(name='(degradation_protein).k1', value=0.05)

# --- Run many stochastic simulations ---
fig, ax = plt.subplots(figsize=(8,5))

for i in range(30):   # 30 stochastic runs
    result = run_time_course(duration=100, step_number=200, method='stochastic', use_numbers=True)
    result.plot(y='mRNA', ax=ax, color='blue', alpha=0.3, legend=None)
    result.plot(y='sRNA', ax=ax, color='red', alpha=0.3, legend=None)
    result.plot(y='protein', ax=ax, color='green', alpha=0.3, legend=None)

ax.set_xlabel("Time (min)")
ax.set_ylabel("Molecule count")
ax.set_title("sRNA–mRNA Regulation (Stochastic Simulations)")
plt.show()
