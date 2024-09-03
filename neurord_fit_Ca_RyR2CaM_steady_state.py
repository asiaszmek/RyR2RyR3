#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
from lxml import etree
import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params

# constants:
kd = {
    "O1": 3/38.4,
    "O1C1": 0.77/0.0025,
    "O2": 3e-3/38.4e-3,
    "O2C1": 0.77e3/2.5,
    "C1I":0.05/11.28,
}


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
        if len(species) > 2:
            continue
        reac_id = son.attrib["id"]
        forward_path = '//Reaction[@id="%s"]/forwardRate' % reac_id
        reverse_path = '//Reaction[@id="%s"]/reverseRate' % reac_id
        print(reac_id, species, forward_path, reverse_path)
        if species[-1].endswith("O1"):
            if "O1" not in counters:
                counters["O1"] = 0
                my_params.append(P('Ca4RyR4_open_fwd_rate', 57.62, min=1e-3,
                                   max=1000,
                                   xpath=forward_path))
            else:
                my_params.append(P('Ca4RyR4_open_fwd_rate%d'%counters["O1"], 0,
                                   fixed='Ca4RyR4_open_fwd_rate',
                                   constant=1,
                                   xpath=forward_path))

            my_params.append(P('Ca4RyR4_open_bkw_rate%d'%counters["O1"], 0,
                               fixed='Ca4RyR4_open_fwd_rate',
                               constant=kd["O1"],
                               xpath=reverse_path))
            
            counters["O1"] += 1
        elif species[0].endswith("O1") and species[-1].endswith("C1"):
            if "O1C1" not in counters:
                counters["O1C1"] = 0
                my_params.append(P('O1_flicker_fwd_rate', 0.139 , min=1e-3,
                                   max=1000,
                                   xpath=forward_path))
            else:
                my_params.append(P('O1_flicker_fwd_rate%d'% counters["O1C1"],
                                   0, fixed='O1_flicker_fwd_rate',
                                   constant=1,
                                   xpath=forward_path))
            my_params.append(P('O1_flicker_bkw_rate%d' % counters["O1C1"], 0,
                               fixed='O1_flicker_fwd_rate',
                               constant=kd["O1C1"], 
                               xpath=reverse_path))
            
            counters["O1C1"] +=1

        elif species[-1].endswith("O2"):
            if "O2" not in counters:
                counters["O2"] = 0
                my_params.append(P('Ca4RyR4_O2_open_fwd_rate', 2.08 , min=1e-3,
                                   max=1000,
                                   xpath=forward_path))
            else:
                my_params.append(P('Ca4RyR4_O2_open_fwd_rate%d'%counters["O2"],
                               0, fixed='Ca4RyR4_O2_open_fwd_rate',
                               constant=1,
                               xpath=forward_path))
            
            my_params.append(P('Ca4RyR4_O2_open_bkw_rate%d'%counters["O2"], 0,
                               fixed='Ca4RyR4_O2_open_fwd_rate',
                               constant=kd["O2"],
                               xpath=reverse_path))
            counters["O2"] += 1

        elif species[0].endswith("O2") and species[-1].endswith("C1"):
            if "O2C1" not in counters:
                counters["O2C1"] = 0
                my_params.append(P('O2_flicker_fwd_rate', 2.57 , min=1e-3,
                                   max=1000,
                                   xpath=forward_path))
            else:
                my_params.append(P('O2_flicker_fwd_rate%d'% counters["O2C1"],
                                   0, fixed='O2_flicker_fwd_rate',
                                   constant=1,
                                   xpath=forward_path))
            my_params.append(P('O2_flicker_bkw_rate%d' % counters["O2C1"], 0,
                               fixed='O2_flicker_fwd_rate',
                               constant=kd["O2C1"], 
                               xpath=reverse_path))
            
            counters["O2C1"] +=1

            

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

