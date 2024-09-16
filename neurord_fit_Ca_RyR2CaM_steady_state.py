#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
from lxml import etree
import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params

# constants:
kd = 4000   # Ca affinity for RyRCaM

dirname='fit_RyR2CaM_Ca/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='model'
exp_set='xu_meissner_ca' #set of data files corresponding to model files; files may contain several molecules

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
my_params = []
counters = {}
for son in root:
    if son.tag == "Reaction":
        species = []
        for grandson in son:
            try:
                species.append(grandson.attrib["specieID"])
            except KeyError:
                continue
    
        reac_id = son.attrib["id"]
        forward_path = '//Reaction[@id="%s"]/forwardRate' % reac_id
        reverse_path = '//Reaction[@id="%s"]/reverseRate' % reac_id
        if species[-1].endswith("O1"):
            continue
        if species[-1].endswith("O2"):
            continue
        if species[-1].endswith("C1"):
            continue
        if species[-1].endswith("I"):
            continue
        if species[-1].startswith("Ca1") and species[0].startswith("RyR"):
            print(species, reac_id)
            if "Ca1" not in counters:
                counters["Ca1"] = 0
                my_params.append(P('Ca1RyR_fwd_rate', 0.001, min=1e-9,
                                    max=1,
                                    xpath=forward_path))
            else:
                my_params.append(P('Ca1RyR_fwd_rate%d'%counters["Ca1"], 0,
                                   fixed='Ca1RyR_fwd_rate',
                                   constant=1,
                                   xpath=forward_path))
            
            my_params.append(P('Ca1RyR_bkw_rate%d'%counters["Ca1"], 0,
                               fixed='Ca1RyR_fwd_rate',
                               constant=kd,
                               xpath=reverse_path))
            counters["Ca1"] += 1
        elif species[-1].startswith("Ca2") and species[0].startswith("Ca1"):
            print(species, reac_id)
                        
            if "Ca2" not in counters:
                counters["Ca2"] = 0
            my_params.append(P('Ca2RyR_fwd_rate%d'%counters["Ca2"], 0,
                                   fixed='Ca1RyR_fwd_rate',
                                   constant=0.75,
                                   xpath=forward_path))
            
            my_params.append(P('Ca2RyR_bkw_rate%d'%counters["Ca2"], 0,
                               fixed='Ca1RyR_fwd_rate',
                               constant=kd*2,
                               xpath=reverse_path))
            counters["Ca2"] += 1
        elif species[-1].startswith("Ca3") and species[0].startswith("Ca2"):
            print(species, reac_id)
                        
            if "Ca3" not in counters:
                counters["Ca3"] = 0
            my_params.append(P('Ca3RyR_fwd_rate%d'%counters["Ca3"], 0,
                                   fixed='Ca1RyR_fwd_rate',
                                   constant=0.5,
                                   xpath=forward_path))
            
            my_params.append(P('Ca3RyR_bkw_rate%d'%counters["Ca3"], 0,
                               fixed='Ca1RyR_fwd_rate',
                               constant=kd*3,
                               xpath=reverse_path))
            counters["Ca3"] += 1
        elif species[-1].startswith("Ca4") and species[0].startswith("Ca3"):
            print(species, reac_id)
                        
            if "Ca4" not in counters:
                counters["Ca4"] = 0
            my_params.append(P('Ca4RyR_fwd_rate%d'%counters["Ca4"], 0,
                                   fixed='Ca1RyR_fwd_rate',
                                   constant=0.25,
                                   xpath=forward_path))
            
            my_params.append(P('Ca4RyR_bkw_rate%d'%counters["Ca4"], 0,
                               fixed='Ca1RyR_fwd_rate',
                               constant=kd*4,
                               xpath=reverse_path))
            counters["Ca4"] += 1

        
        elif species[-1].endswith("I2"):
            print(species, reac_id)
            my_params.append(P("I2_flicker_fwd_rate", 1,
                               min=1e-3, max=1e3,
                               xpath=forward_path))
            
            my_params.append(P("I2_flicker_bkw_rate", 0.06,
                               min=1e-3, max=1e3,
                               xpath=reverse_path))
                


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

