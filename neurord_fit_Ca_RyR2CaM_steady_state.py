#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os
kd = 1000
dirname='fit_RyR2CaM_Ca/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model'
exp_set='xu_meissner_ca' #set of data files corresponding to model files; files may contain several molecules
mol={"RO": ["Ca4_RyR4_4CaMCa2C_O1",
            "Ca4_RyR4_4CaMCa2C_O2",
            "Ca4_RyR4_1CaM_3CaMCa2C_O1",
            "Ca4_RyR4_1CaM_3CaMCa2C_O2",
            "Ca4_RyR4_2CaM_2CaMCa2C_O1",
            "Ca4_RyR4_2CaM_2CaMCa2C_O2",
            "Ca4_RyR4_3CaM_1CaMCa2C_O1",
            "Ca4_RyR4_3CaM_1CaMCa2C_O2",
            "Ca4_RyR4_4CaM_O1",
            "Ca4_RyR4_4CaM_O2"]} #which molecule(s) to match in optimization
tmpdir='/tmp/RyR2CaM_Ca_ss'+dirname 
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

    P('Ca4RyR4_open_fwd_rate', 5.21, min=1e-3, max=1000,
      xpath='//Reaction[@id="141"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate', 8.08,  min=1e-3, max=1000,
      xpath='//Reaction[@id="141"]/reverseRate'),
    P('Ca4RyR4_open_fwd_rate1', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="146"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate1', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1,
      xpath='//Reaction[@id="146"]/reverseRate'),
    P('Ca4RyR4_open_fwd_rate2', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="151"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate2', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1,
      xpath='//Reaction[@id="151"]/reverseRate'),
    P('Ca4RyR4_open_fwd_rate3', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="156"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate3', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1,
      xpath='//Reaction[@id="156"]/reverseRate'),
    P('Ca4RyR4_open_fwd_rate4', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1,
      xpath='//Reaction[@id="161"]/forwardRate'),
    P('Ca4RyR4_open_bkw_rate4', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1,
      xpath='//Reaction[@id="161"]/reverseRate'),


    
    P('O1_flicker_fwd_rate', 0.16 , min=1e-6, max=1000,
      xpath='//Reaction[@id="143"]/forwardRate'),
    P('O1_flicker_bkw_rate', 31.61, min=1e-3, max=1000,
      xpath='//Reaction[@id="143"]/reverseRate'),
    
    P('O1_flicker_fwd_rate1', 0, fixed='O1_flicker_fwd_rate',
      constant=1, 
      xpath='//Reaction[@id="148"]/forwardRate'),
    P('O1_flicker_bkw_rate1', 0, fixed='O1_flicker_bkw_rate',
      constant=1, 
      xpath='//Reaction[@id="148"]/reverseRate'),
    P('O1_flicker_fwd_rate2', 0, fixed='O1_flicker_fwd_rate',
      constant=1, 
      xpath='//Reaction[@id="153"]/forwardRate'),
    P('O1_flicker_bkw_rate2', 0, fixed='O1_flicker_bkw_rate',
      constant=1,  
      xpath='//Reaction[@id="153"]/reverseRate'),
    P('O1_flicker_fwd_rate3', 0, fixed='O1_flicker_fwd_rate',
      constant=1, 
      xpath='//Reaction[@id="158"]/forwardRate'),
    P('O1_flicker_bkw_rate3', 0, fixed='O1_flicker_bkw_rate',
      constant=1, 
      xpath='//Reaction[@id="158"]/reverseRate'),
    P('O1_flicker_fwd_rate4', 0, fixed='O1_flicker_fwd_rate',
      constant=1, 
      xpath='//Reaction[@id="163"]/forwardRate'),
    P('O1_flicker_bkw_rate4', 0, fixed='O1_flicker_bkw_rate',
      constant=1, 
      xpath='//Reaction[@id="163"]/reverseRate'),

    
    P('Ca4RyR4_O2_open_fwd_rate', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1e-2,
      xpath='//Reaction[@id="142"]/forwardRate'),
    P('Ca4RyR4_O2_open_bkw_rate', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1e-2,
      xpath='//Reaction[@id="142"]/reverseRate'),
    P('Ca4RyR4_O2_open_fwd_rate1', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1e-2,
      xpath='//Reaction[@id="147"]/forwardRate'),
    P('Ca4RyR4_O2_open_bkw_rate1', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1e-2,
      xpath='//Reaction[@id="147"]/reverseRate'),
    P('Ca4RyR4_O2_open_fwd_rate2', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1e-2,
      xpath='//Reaction[@id="152"]/forwardRate'),
    P('Ca4RyR4_O2_open_bkw_rate2', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1e-2,
      xpath='//Reaction[@id="152"]/reverseRate'),
    P('Ca4RyR4_O2_open_fwd_rate3', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1e-2,
      xpath='//Reaction[@id="157"]/forwardRate'),
    P('Ca4RyR4_O2_open_bkw_rate3', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1e-2,
      xpath='//Reaction[@id="157"]/reverseRate'),
    P('Ca4RyR4_O2_open_fwd_rate4', 0, fixed='Ca4RyR4_open_fwd_rate',
      constant=1e-2,
      xpath='//Reaction[@id="162"]/forwardRate'),
    P('Ca4RyR4_O2_open_bkw_rate3', 0, fixed='Ca4RyR4_open_bkw_rate',
      constant=1e-2,
      xpath='//Reaction[@id="162"]/reverseRate'),

    
    P('O2_flicker_fwd_rate',  0, fixed='O1_flicker_fwd_rate',
      constant=1e2,
       xpath='//Reaction[@id="144"]/forwardRate'),
    P('O2_flicker_bkw_rate', 0, fixed='O1_flicker_bkw_rate',
      constant=1e2,
      xpath='//Reaction[@id="144"]/reverseRate'),
    P('O2_flicker_fwd_rate1',  0, fixed='O1_flicker_fwd_rate',
      constant=1e2,
       xpath='//Reaction[@id="149"]/forwardRate'),
    P('O2_flicker_bkw_rate1', 0, fixed='O1_flicker_bkw_rate',
      constant=1e2,
      xpath='//Reaction[@id="149"]/reverseRate'),
    P('O2_flicker_fwd_rate2',  0, fixed='O1_flicker_fwd_rate',
      constant=1e2,
       xpath='//Reaction[@id="154"]/forwardRate'),
    P('O2_flicker_bkw_rate2', 0, fixed='O1_flicker_bkw_rate',
      constant=1e2,
      xpath='//Reaction[@id="154"]/reverseRate'),
    P('O2_flicker_fwd_rate3',  0, fixed='O1_flicker_fwd_rate',
      constant=1e2,
      xpath='//Reaction[@id="159"]/forwardRate'),
    P('O2_flicker_bkw_rate3', 0, fixed='O1_flicker_bkw_rate',
      constant=1e2,
      xpath='//Reaction[@id="159"]/reverseRate'),
    P('O2_flicker_fwd_rate4',  0, fixed='O1_flicker_fwd_rate',
      constant=1e2,
      xpath='//Reaction[@id="164"]/forwardRate'),
    P('O2_flicker_bkw_rate4', 0, fixed='O1_flicker_bkw_rate',
      constant=1e2,
      xpath='//Reaction[@id="164"]/reverseRate'),


   
    P('O2_I_flicker_fwd_rate', 4.37, min=1e-3, max=1000,
      xpath='//Reaction[@id="145"]/forwardRate'),
    P('O2_I_flicker_bkw_rate', 0.94, min=1e-3, max=1000,
      xpath='//Reaction[@id="145"]/reverseRate'),
    P('O2_I_flicker_fwd_rate1', 0, fixed="O2_I_flicker_fwd_rate",
      constant=1,
      xpath='//Reaction[@id="150"]/forwardRate'),
    P('O2_I_flicker_bkw_rate1', 0, fixed="O2_I_flicker_bkw_rate",
      constant=1,
      xpath='//Reaction[@id="150"]/reverseRate'),
    P('O2_I_flicker_fwd_rate2', 0, fixed="O2_I_flicker_fwd_rate",
      constant=1,
      xpath='//Reaction[@id="155"]/forwardRate'),
    P('O2_I_flicker_bkw_rate2', 0, fixed="O2_I_flicker_bkw_rate",
      constant=1,
      xpath='//Reaction[@id="155"]/reverseRate'),
    P('O2_I_flicker_fwd_rate3', 0, fixed="O2_I_flicker_fwd_rate",
      constant=1,
      xpath='//Reaction[@id="160"]/forwardRate'),
    P('O2_I_flicker_bkw_rate3', 0, fixed="O2_I_flicker_bkw_rate",
      constant=1,
      xpath='//Reaction[@id="160"]/reverseRate'),
    P('O2_I_flicker_fwd_rate4', 0, fixed="O2_I_flicker_fwd_rate",
      constant=1,
      xpath='//Reaction[@id="165"]/forwardRate'),
    P('O2_I_flicker_bkw_rate4', 0, fixed="O2_I_flicker_bkw_rate",
      constant=1,
      xpath='//Reaction[@id="165"]/reverseRate'),

 
    
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

