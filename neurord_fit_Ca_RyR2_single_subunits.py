#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os

dirname='fit_RyR2_Ca_single_subunits/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model_Z'
exp_set='po' #set of data files corresponding to model files; files may contain several molecules
mol={"RO": ["O1", "O2"]} #which molecule(s) to match in optimization
tmpdir='/tmp/RyR2_single_subunits'+dirname 
os.chdir(dirname)
kd = 1000
# Use loadconc.CSV_conc_set if data to match are csv format (typically from wet experiments)
exp = loadconc.CSV_conc_set(exp_set)

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=100
popsize=8
test_size=25 #for convergence

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(
    
    P('Ca4RyR4_open_fwd_rate',4e-7 , min=1e-12, max=1e-3,
      xpath='//Reaction[@id="RyR4Ca4"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate', 0.0025,  min=1e-5, max=1,
      xpath='//Reaction[@id="RyR4Ca4"]/reverseRate'),

    # P('Ca4RyR4_open_fwd_rate', 38.44, min=1e-3, max=1000,
    #   xpath='//Reaction[@id="RyRd"]/forwardRate'),
    # P('Ca4RyR4_open_bkw_rate', 3.02,  min=1e-3, max=1000,
    #   xpath='//Reaction[@id="RyRd"]/reverseRate'),

    # P('O1_flicker_fwd_rate', 0.0025 , min=1e-6, max=1000,
    #   xpath='//Reaction[@id="RyRf"]/forwardRate'),
    # P('O1_flicker_bkw_rate', 0.77, min=1e-3, max=1000,
    #   xpath='//Reaction[@id="RyRf"]/reverseRate'),

    # P('Ca4RyR4_O2_open_fwd_rate', 0, fixed='Ca4RyR4_open_fwd_rate',
    #   constant=1e-3,
    #   xpath='//Reaction[@id="RyRe"]/forwardRate'),
    # P('Ca4RyR4_O2_open_bkw_rate', 0, fixed='Ca4RyR4_open_bkw_rate',
    #   constant=1e-3,
    #   xpath='//Reaction[@id="RyRe"]/reverseRate'),

    # P('O2_flicker_fwd_rate',  0, fixed='O1_flicker_fwd_rate',
    #   constant=1e3,
    #    xpath='//Reaction[@id="RyRg"]/forwardRate'),
    # P('O2_flicker_bkw_rate', 0, fixed='O1_flicker_bkw_rate',
    #   constant=1e3,
    #   xpath='//Reaction[@id="RyRg"]/reverseRate'),
    # P('O2_I_flicker_fwd_rate', 11.28, min=1e-3, max=1000,

    #   xpath='//Reaction[@id="RyRh"]/forwardRate'),
    # P('O2_I_flicker_bkw_rate', 0.05, min=1e-3, max=1000,
      
    #   xpath='//Reaction[@id="RyRh"]/reverseRate'),


    
    
)

###################### END CUSTOMIZATION #######################################

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol)

############ Test fitness function
#model=dirname+'Model-CKnew-Cahz1.xml'
#sim = aju.xml.NeurordSimulation('/tmp', model=model, params=params)
#sim2=aju.xml.NeurordResult('Model_syngap_ras.h5')
#print(fitness(sim2, exp))
################

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, popsize=popsize,sigma=0.3)
mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)
if callable(fit.optimizer.result):
    result = fit.optimizer.result()
else:
    result = fit.optimizer.result
#print centroid [0] and stdev [6] of cloud of good results:
for i,p in enumerate(fit.params.unscale(result[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(result[6])[i])

save_params.save_params(fit,0,1)

