import itertools
from lxml import etree

fname = "Rxn_module_RyR3_CaM.xml"

#  Keep CaM binding constants, it's CaMBD2, which is probably the most conserved
#  Does CaM activate RyR3? It's a pretty common sentence in literature and
#  apparently it means that, when CaM is added RyR3 is open more often




kfs = {"CaM": 2.1e-8, # for Kd of 820 nM (Xu and Meissner 2004 for RyR2)
       "CaMCa2C": 3.15e-7, "CaMCa4": 3.66e-7, "2CaC": 6e-5,
       "2CaN":  0.1e-2, "RyR3Ca1": 1e-3, "RyR3Ca2": 0.75e-3,
       "RyR3Ca3":5e-4, "RyR3Ca4": 2.5e-4, "RyR3Ca4O1": 38.4, "RyR3Ca4O1C1": 0.0025,
       "RyR3Ca4O2": 38.4e-3, "RyR3Ca4O2C1":2.5, "RyR3Ca4C1I": 11.28,
       "CaMRyR3Ca1": 1e-3, "CaMRyR3Ca2": 0.75e-3,
       "CaMRyR3Ca3":5e-4, "CaMRyR3Ca4": 2.5e-4, "CaMRyR3Ca4O1": 5.21,
       "CaMRyR3Ca4O1C1": 0.16,
       "CaMRyR3Ca4O2": 5.21e-3, "CaMRyR3Ca4O2C1":0.16e3, "CaMRyR3Ca4C1I": 4.37}


krs = {"CaM": 1.73e-5, "CaMCa2C": 2.59e-5, "CaMCa4": 3.015e-6,
       "2CaC": 9.1e-3, "2CaN": 1000e-3,"RyR3Ca1": 1,
       "RyR3Ca2": 2, "RyR3Ca3": 3, "RyR3Ca4": 4,
       "RyR3Ca4O1": 3,"RyR3Ca4O1C1": 0.77, "RyR3Ca4O2": 3e-3,
       "RyR3Ca4O2C1": 0.77e3,  "RyR3Ca4C1I":0.05,
       "CaMRyR3Ca1": 1,
       "CaMRyR3Ca2": 2, "CaMRyR3Ca3": 3, "CaMRyR3Ca4": 4,
       "CaMRyR3Ca4O1": 8.08,"CaMRyR3Ca4O1C1": 32, "CaMRyR3Ca4O2": 8.08e-3,
       "CaMRyR3Ca4O2C1": 32e3,  "CaMRyR3Ca4C1I":0.95}
counter = 1

def add_reaction(root, name, what, new_name):
    global counter
    my_r = etree.SubElement(root, "Reaction",
                            name=name+"_"+what+"_"+str(counter),
                            id=name+"_"+what+"_"+str(counter))
    etree.SubElement(my_r, "Reactant", specieID=name)
    if what == "2CaC" or what == "2CaN":
        etree.SubElement(my_r, "Reactant", specieID="Ca", n="2")
    elif what in ["CaM", "CaMCa2C", "CaMCa4"]:
        etree.SubElement(my_r, "Reactant", specieID=what)
    elif "O1" in what or "O2" in what or "C1" in what or "I" in what:
        pass
    else:
        etree.SubElement(my_r, "Reactant", specieID="Ca")
    etree.SubElement(my_r, "Product", specieID=new_name)
    kf = etree.SubElement(my_r, "forwardRate")
    kf.text = str(kfs[what])
    kr = etree.SubElement(my_r, "reverseRate")
    kr.text = str(krs[what])
    q = etree.SubElement(my_r, "Q10")
    q.text = ".2"
    counter += 1

    
def generate_name(a, b, c, d=0):
    if not d:
        name = "RyR3"
    else:
        name = "Ca%d_RyR3" % d
    if a:
        name += "_%dCaM" % a
    if b:
        name += "_%dCaMCa2C" % b
    if c:
        name += "_%dCaMCa4" % c
    return name

if __name__ == "__main__":
    states = set()
    for i in range(0,5):
        for j in range(0, 5):
            for k in range(0, 5):
                if i+j+k > 4:
                    continue
                for iteration in itertools.permutations((i,j,k), r=3):
                    states.add(iteration)

    
    my_rxn_file = etree.Element("ReactionScheme")
    etree.SubElement(my_rxn_file, "Specie", name="Ca",
                         id="Ca", kdiff="200", kdiffunit="mu2/s")

    etree.SubElement(my_rxn_file, "Specie", name="CaM",
                         id="CaM", kdiff="4", kdiffunit="mu2/s")
    etree.SubElement(my_rxn_file, "Specie", name="CaMCa2C",
                         id="CaMCa2C", kdiff="4", kdiffunit="mu2/s")
    etree.SubElement(my_rxn_file, "Specie", name="CaMCa4",
                         id="CaMCa4", kdiff="4", kdiffunit="mu2/s")
    my_r = etree.SubElement(my_rxn_file, "Reaction",
                            name="CaM_Ca2C",
                            id="CaM_Ca2C")
    etree.SubElement(my_r, "Reactant", specieID="CaM")
    etree.SubElement(my_r, "Reactant", specieID="Ca", n="2")
    etree.SubElement(my_r, "Product", specieID="CaMCa2C")
    kf = etree.SubElement(my_r, "forwardRate")
    kf.text = str(6e-6)
    kr = etree.SubElement(my_r, "reverseRate")
    kr.text = str(9.1e-3)
    q = etree.SubElement(my_r, "Q10")
    q.text = ".2"

    my_r = etree.SubElement(my_rxn_file, "Reaction",
                            name="CaMCa2C_2Ca",
                            id="CaMCa2C_2Ca")
    etree.SubElement(my_r, "Reactant", specieID="CaMCa2C")
    etree.SubElement(my_r, "Reactant", specieID="Ca", n="2")
    etree.SubElement(my_r, "Product", specieID="CaMCa4")
    kf = etree.SubElement(my_r, "forwardRate")
    kf.text = str(0.1e-3)
    kr = etree.SubElement(my_r, "reverseRate")
    kr.text = str(1000e-3)
    q = etree.SubElement(my_r, "Q10")
    q.text = ".2"
    ryr_species_to_open =[]


    for l in [0, 1, 2, 3, 4]:
        for (i,  j, k) in sorted(states):
            my_specie_name = generate_name(i, j, k, l)
            etree.SubElement(my_rxn_file, "Specie", name=my_specie_name,
                             id=my_specie_name, kdiff="0", kdiffunit="mu2/s")
            if my_specie_name.startswith("Ca4_"):
                ryr_species_to_open.append(my_specie_name)

    for i, specie in enumerate(ryr_species_to_open):
        etree.SubElement(my_rxn_file, "Specie", name="%s_O1" % specie,
                         id="%s_O1" % specie, kdiff="0", kdiffunit="mu2/s")
        etree.SubElement(my_rxn_file, "Specie", name="%s_O2" % specie,
                         id="%s_O2" % specie, kdiff="0", kdiffunit="mu2/s")
        etree.SubElement(my_rxn_file, "Specie", name="%s_C1" % specie,
                         id="%s_C1" % specie, kdiff="0",
                         kdiffunit="mu2/s")
        etree.SubElement(my_rxn_file, "Specie", name="%s_I" % specie,
                         id="%s_I" % specie, kdiff="0",
                         kdiffunit="mu2/s")


    for l in [0, 1, 2, 3, 4]:
        for (i,  j, k) in sorted(states):
            my_specie_name = generate_name(i, j, k, l)
            if i+j+k < 4:
                add_reaction(my_rxn_file, my_specie_name, "CaM",
                             generate_name(i+1,j, k, l))
                add_reaction(my_rxn_file, my_specie_name, "CaMCa2C",
                             generate_name(i,j+1, k, l))
                add_reaction(my_rxn_file, my_specie_name, "CaMCa4",
                             generate_name(i,j, k+1, l))
                if my_specie_name in ryr_species_to_open:
                    for suffix in ["O1", "O2", "C1", "I"]:
                        my_name = "%s_%s" % (my_specie_name, suffix)
                        add_reaction(my_rxn_file, my_name, "CaM",
                                     "%s_%s" % (generate_name(i+1,j, k, l),
                                                suffix))
                        add_reaction(my_rxn_file, my_name, "CaMCa2C",
                                    "%s_%s" % (generate_name(i,j+1, k, l),
                                               suffix))
                        add_reaction(my_rxn_file, my_name, "CaMCa4",
                                     "%s_%s"% (generate_name(i,j, k+1, l),
                                               suffix))
            else:
                pass   
            if not i:
                pass 
            else:
                add_reaction(my_rxn_file, my_specie_name, "2CaC",
                             generate_name(i-1, j+1, k, l))
                if my_specie_name in ryr_species_to_open:
                    for suffix in ["O1", "O2", "C1", "I"]:
                        my_name = "%s_%s" % (my_specie_name, suffix)
                        add_reaction(my_rxn_file, my_name, "2CaC",
                             "%s_%s" % (generate_name(i-1, j+1, k, l),
                                        suffix))
            
            if not j:    
                pass  
            else:
                add_reaction(my_rxn_file, my_specie_name, "2CaN",
                             generate_name(i, j-1, k+1, l))
                if my_specie_name in ryr_species_to_open:
                    for suffix in ["O1", "O2", "C1", "I"]:
                        my_name = "%s_%s" % (my_specie_name, suffix)
                        add_reaction(my_rxn_file, my_name, "2CaN",
                                     "%s_%s" % (generate_name(i, j-1, k+1, l),
                                                suffix))

            if l < 4:
                new_name = generate_name(i, j, k, l + 1)
                if "CaM" in my_specie_name:
                    rxn_name = "CaMRyR3Ca%d" % (l+1)
                else:
                    rxn_name = "RyR3Ca%d" % (l+1)
                add_reaction(my_rxn_file, my_specie_name, rxn_name,
                             new_name)

                

    for i, specie in enumerate(ryr_species_to_open):
        if "CaM" not in specie:
            add_reaction(my_rxn_file, specie, "RyR3Ca4O1",
                         "%s_O1"%specie)
            add_reaction(my_rxn_file, specie, "RyR3Ca4O2",
                         "%s_O2"%specie)
            add_reaction(my_rxn_file, "%s_O1" % specie,"RyR3Ca4O1C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_O2" % specie,"RyR3Ca4O2C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_C1" % specie,"RyR3Ca4C1I",
                         "%s_I"%specie)
        else:
            add_reaction(my_rxn_file, specie, "CaMRyR3Ca4O1",
                         "%s_O1"%specie)
            add_reaction(my_rxn_file, specie, "CaMRyR3Ca4O2",
                         "%s_O2"%specie)
            add_reaction(my_rxn_file, "%s_O1" % specie,"CaMRyR3Ca4O1C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_O2" % specie,"CaMRyR3Ca4O2C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_C1" % specie,"CaMRyR3Ca4C1I",
                         "%s_I"%specie)

                
    f = open(fname, "w")
    f.write(etree.tostring(my_rxn_file, pretty_print=True).decode("utf-8"))
    f.close()
    print(ryr_species_to_open)
