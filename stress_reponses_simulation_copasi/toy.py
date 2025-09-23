# # import sys
# # if '../..' not in sys.path:
# #     sys.path.append('../..')
# #
# # import basico.biomodels as biomodels
# #
# # glycolysis_models = biomodels.search_for_model('glycolysis')
# # for model in glycolysis_models:
# #     print ('Id: %s' % model['id'])
# #     print ('Name: %s' % model['name'])
# #     print ('Format: %s' % model['format'])
# #     print ('')
# from basico import *
# import sys
# if '../..' not in sys.path:
#     sys.path.append('../..')
#
# import basico.biomodels as biomodels
#
# # Search for RpoS-related models
# rpos_models = biomodels.search_for_model('RpoS')
# for model in rpos_models:
#     print ('Id: %s' % model['id'])
#     print ('Name: %s' % model['name'])
#     print ('Format: %s' % model['format'])
#     print ('')
#
# # Search for small RNA models
# srna_models = biomodels.search_for_model('small RNA')
# for model in srna_models:
#     print ('Id: %s' % model['id'])
#     print ('Name: %s' % model['name'])
#     print ('Format: %s' % model['format'])
#     print ('')
#
# # Search for stress response models
# stress_models = biomodels.search_for_model('stress response')
# for model in stress_models:
#     print ('Id: %s' % model['id'])
#     print ('Name: %s' % model['name'])
#     print ('Format: %s' % model['format'])
#     print ('')
from matplotlib import pyplot as plt
from basico import *
biomod = load_biomodel(10)
print(biomod)

tc = run_time_course(duration = 100)
tc.plot()
plt.show()