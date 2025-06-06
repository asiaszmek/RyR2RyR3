#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
from lxml import etree
import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params

# constants:

dirname='fit_KSRyR3_Ca/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model'
exp_set='chen_et_al' #set of data files corresponding to model files; files may contain several molecules

#which molecule(s) to match in optimization
fname_xml = os.path.join(dirname, "Rxn_RyRCaM.xml")
tree = etree.parse(fname_xml)
root = tree.getroot()
mol_list = []
for son in root:
    if son.tag == "Specie":
        if "O" in son.attrib["id"]:
            mol_list.append(son.attrib["id"])

mol={"RO": mol_list}

print(mol)

tmpdir='results'+dirname 
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


my_params = []


my_params.append(P('RyR3CaMe_f', 1.73e-13, min=1e-20,
                   max=1e-8,
                   xpath='//Reaction[@id="RyR3e"]/forwardRate'))
my_params.append(P('RyR3CaMe_r', 0.049, min=1e-3,
                   max=1e1,
                   xpath='//Reaction[@id="RyR3e"]/reverseRate'))

my_params.append(P('RyR3CaMa_f', 1.355e-17, min=1e-20,
                   max=1e-8,
                   xpath='//Reaction[@id="RyR3a"]/forwardRate'))
my_params.append(P('RyR3CaMa_r', 0.068, min=1e-3,
                   max=1e2,
                   xpath='//Reaction[@id="RyR3a"]/reverseRate'))
my_params.append(P('RyR3CaMb_f', 5e-12, min=1e-20,
                   max=1e-8,
                   xpath='//Reaction[@id="RyR3b"]/forwardRate'))
my_params.append(P('RyR3CaMb_r', 13, min=1e-3,
                   max=1e2,
                   xpath='//Reaction[@id="RyR3b"]/reverseRate'))

my_params.append(P('RyR3CaMd_f', 0, fixed="RyR3CaMa_f", constant=1e-2,
                   xpath='//Reaction[@id="RyR3d"]/forwardRate'))
my_params.append(P('RyR3CaMd_r',0, fixed="RyR3CaMa_r", constant=1e-2,
                   xpath='//Reaction[@id="RyR3d"]/reverseRate'))
my_params.append(P('RyR3CaMc_f', 0, fixed="RyR3CaMb_f", constant=1e2,
                   xpath='//Reaction[@id="RyR3c"]/forwardRate'))
my_params.append(P('RyR3CaMc_r', 0, fixed="RyR3CaMb_f", constant=1e2,
                   xpath='//Reaction[@id="RyR3c"]/reverseRate'))



my_params.append(P('RyR3CaMf_f', 5.89e-19, min=1e-20,
                   max=1e-10,
                   xpath='//Reaction[@id="RyR3f"]/forwardRate'))
my_params.append(P('RyR3CaMf_r', 0.023, 
                   min=1e-4,
                   max=1, 
                   xpath='//Reaction[@id="RyR3f"]/reverseRate'))



params = aju.optimize.ParamSet(*my_params)
                             


# ###################### END CUSTOMIZATION #######################################

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

