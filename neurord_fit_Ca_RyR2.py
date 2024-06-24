#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os

dirname='fit_RyR_Ca/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model'
exp_set='po' #set of data files corresponding to model files; files may contain several molecules
mol={"RO": ["Ca3RyR4O", "Ca4RyR4O" ]} #which molecule(s) to match in optimization
tmpdir='/tmp/RyR2'+dirname 
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
    
    
    P('3Ca3RyR_RyRbinding_fwd_rate', 1, min=1e-3, max=1000,
      xpath='//Reaction[@id="RyRa"]/forwardRate'),
    P('3Ca3RyR_RyRbinding_bkw_rate', 0, fixed='3Ca3RyR_RyRbinding_fwd_rate',
      constant=1000,
      xpath='//Reaction[@id="RyRa"]/reverseRate'),
    P('Ca3RyR3_RyR_fwd', 0, fixed= '3Ca3RyR_RyRbinding_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="RyRa2"]/forwardRate'),
    P('Ca3RyR3_RyR_bkw', 0, fixed='3Ca3RyR_RyRbinding_fwd_rate',
      constant=0.01,
      xpath='//Reaction[@id="RyRa2"]/reverseRate'),
    P('Ca3RyR3_CaRyR_fwd', 0, fixed= '3Ca3RyR_RyRbinding_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="RyRa3"]/forwardRate'),
    P('Ca3RyR3_CaRyR_bkw', 0, fixed='3Ca3RyR_RyRbinding_fwd_rate',
      constant=0.01,
      xpath='//Reaction[@id="RyRa3"]/reverseRate'),
   
    
    P('Ca3RyR4_flicker_fwd_rate', 1, min=1e-6, max=10,
      xpath='//Reaction[@id="RyRb"]/forwardRate'),
    P('Ca3RyR4_flicker_bkw_rate', 1, min=1e-6, max=10,
      xpath='//Reaction[@id="RyRb"]/reverseRate'),

    P('Ca4RyR4_flicker_fwd_rate', 0, fixed="Ca3RyR4_flicker_fwd_rate",
      constant=1,
      xpath='//Reaction[@id="RyRd"]/forwardRate'),
    P('Ca4RyR4_flicker_bkw_rate', 0, fixed="Ca3RyR4_flicker_bkw_rate",
      constant=1,
      xpath='//Reaction[@id="RyRd"]/reverseRate'),

    P('O_inact_fwd_rate', 5e-4, min=1e-6,max=1,
      xpath='//Reaction[@id="RyRbi"]/forwardRate'),
    P('O_inact_bkw_rate', 1.5e-3, min=1e-6,max=1,
      xpath='//Reaction[@id="RyRbi"]/reverseRate'),

    
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

