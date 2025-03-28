#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os

dirname='CaM/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model_new_reduced'
exp_set='model_old_reduced' #set of data files corresponding to model files; files may contain several molecules
mol={"CaM": ["CaM"],
     "CaMCa2C": ["CaMCa2C"],
     "CaMCa2N": ["CaMCa2N"],
     "CaMCa4": ["CaMCa4"]} #which molecule(s) to match in optimization
tmpdir='/tmp/CaM'+dirname 
os.chdir(dirname)
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
    P('CaM_Ca_2C_bind_fwd', 7.5e-8, min=1e-12, max=1e-6,
      xpath='//Reaction[@id="CaMC_bind"]/forwardRate'),
    P('CaM_Ca_2C_bind_bkw', 9.1e-3, min=1e-6, max=1,
      xpath='//Reaction[@id="CaMC_bind"]/reverseRate'),
    
    P('CaM_Ca_2C_bind_fwd2',  0, fixed='CaM_Ca_2C_bind_fwd', constant=1,
      xpath='//Reaction[@id="CaMCa2N_bind"]/forwardRate'),
    P('CaM_Ca_2C_bind_bkw2', 0, fixed='CaM_Ca_2C_bind_bkw', constant=1,
      xpath='//Reaction[@id="CaMCa2N_bind"]/reverseRate'),
    
    P('CaM_Ca_2N_bind_fwd', 0.01e-4, min=1e-7, max=1e-4,
      xpath='//Reaction[@id="CaMN_bind"]/forwardRate'),
    P('CaM_Ca_2N_bind_bkw', 1, min=1e-3, max=100,
      xpath='//Reaction[@id="CaMN_bind"]/reverseRate'),
    
    P('CaM_Ca_2N_bind_fwd2', 0, fixed='CaM_Ca_2N_bind_fwd', constant=1,
      xpath='//Reaction[@id="CaMCa2C_bind"]/forwardRate'),
    P('CaM_Ca_2N_bind_bkw2', 0, fixed="CaM_Ca_2N_bind_bkw", constant=1,
      xpath='//Reaction[@id="CaMCa2C_bind"]/reverseRate'),
    
    
    
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

