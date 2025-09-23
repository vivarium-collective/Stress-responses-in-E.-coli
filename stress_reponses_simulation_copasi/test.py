import sys
from basico import *
import numpy as np
import matplotlib.pyplot as plt

new_model(name='Simple Model')

add_reaction('R1', 'A -> B')

get_species().initial_concentration

set_species('B', initial_concentration=0)
set_species('A', initial_concentration=10)
get_species().initial_concentration

get_reaction_parameters()

set_reaction_parameters('(R1).k1', value=1)
get_reaction_parameters('k1')

result = run_time_course(duration=50)
result.plot()
plt.show()
print(result)