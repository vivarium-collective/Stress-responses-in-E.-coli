
# clone the repositories
git clone https://github.com/vivarium-collective/stress-responses-tellurium.git
cd stress-responses-tellurium

#Create a Python environment and install dependencies:

conda create -n stress_env python=3.10
conda activate stress_env
pip install tellurium basico matplotlib pandas

#description
This repository contains Tellurium-based models of **sRNA regulation of RpoS** 
and sigma factor competition in *E. coli*.  
This repository contains models and scripts for simulating stress responses in Escherichia coli, focusing on oxidative stress (OxyR, SoxRS), envelope stress (RpoE), and general stress response (RpoS).
Models are written in Antimony/Tellurium for transparency and can be exported to SBML for use in COPASI, BasiCO,
It is part of the [Vivarium Collective](https://github.com/vivarium-collective).

[oxyR_SoxR.py](stress_responses_simulation_copasi/model/oxyR_SoxR.py)-OxyR regulon responding to hydrogen peroxide.

 

models/soxrs.ant → SoxRS regulon responding to superoxide stress.

models/rpos.ant → General stress response sigma factor.

models/rpoe.ant → Envelope stress response.



References

Bouillet et al. (2024) RpoS and the bacterial general stress response.

BioModels database: curated SBML models of stress responses.

Vivarium Collective: https://vivarium-collective.github.io