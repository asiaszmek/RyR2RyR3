import itertools
from lxml import etree

fname = "Rxn_RyRCaM.xml"

kfs = {"CaM": 2.1e-8, # for Kd of 820 nM (Xu and Meissner 2004 for RyR2)
       "CaMCa2C": 3.15e-7, "CaMCa4": 3.66e-7, "2CaC": 6e-5,
       "2CaN":  0.1e-2, "RyR2Ca1": 1e-3, "RyR2Ca2": 0.75e-3,
       "RyR2Ca3":5e-4, "RyR2Ca4": 2.5e-4, "RyR2Ca4O1": 38.4, "RyR2Ca4O1C1": 0.0025,
       "RyR2Ca4O2": 38.4e-3, "RyR2Ca4O2C1":2.5, "RyR2Ca4C1I": 11.28,
       "CaMRyR2Ca1": 1e-3, "CaMRyR2Ca2": 0.75e-3,
       "CaMRyR2Ca3":5e-4, "CaMRyR2Ca4": 2.5e-4, "CaMRyR2Ca4O1": 5.21,
       "CaMRyR2Ca4O1C1": 0.16,
       "CaMRyR2Ca4O2": 5.21e-2, "CaMRyR2Ca4O2C1":0.16e2, "CaMRyR2Ca4C1I": 4.37}


krs = {"CaM": 1.73e-5, "CaMCa2C": 3.67e-6, "CaMCa4": 3.015e-6,
       "2CaC": 9.1e-3, "2CaN": 1000e-3,"RyR2Ca1": 1,
       "RyR2Ca2": 2, "RyR2Ca3": 3, "RyR2Ca4": 4,
       "RyR2Ca4O1": 3,"RyR2Ca4O1C1": 0.77, "RyR2Ca4O2": 3e-3,
       "RyR2Ca4O2C1": 0.77e3,  "RyR2Ca4C1I":0.05,
       "CaMRyR2Ca1": 1,
       "CaMRyR2Ca2": 2, "CaMRyR2Ca3": 3, "CaMRyR2Ca4": 4,
       "CaMRyR2Ca4O1": 8.08,"CaMRyR2Ca4O1C1": 32, "CaMRyR2Ca4O2": 8.08e-2,
       "CaMRyR2Ca4O2C1": 32e2,  "CaMRyR2Ca4C1I":0.95}
counter = 1

def add_reaction(root, name, what, new_name):
    global counter
    my_r = etree.SubElement(root, "Reaction",
                            name=str(counter),
                            id=str(counter))
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
        name = "RyR4"
    else:
        name = "Ca%d_RyR4" % d
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
                if i+j+k != 4:
                    continue
                for iteration in itertools.permutations((i,j,k), r=3):
                    states.add(iteration)

    
    my_rxn_file = etree.Element("ReactionScheme")
    etree.SubElement(my_rxn_file, "Specie", name="Ca",
                         id="Ca", kdiff="200", kdiffunit="mu2/s")

    ryr_species_to_open =[]


    for l in [0, 1, 2, 3, 4]:
        for (i,  j, k) in sorted(states):
            my_specie_name = generate_name(i, j, k, l)
            etree.SubElement(my_rxn_file, "Specie", name=my_specie_name,
                             id=my_specie_name, kdiff="0", kdiffunit="mu2/s")
            if my_specie_name.startswith("Ca4_"):
                ryr_species_to_open.append(my_specie_name)

    for l in [0, 1, 2, 3]:
        for (i,  j, k) in sorted(states):
            my_specie_name = generate_name(i, j, k, l)
    
            if i+j+k < 4:
                add_reaction(my_rxn_file, my_specie_name, "CaM",
                             generate_name(i+1,j, k, l))
                add_reaction(my_rxn_file, my_specie_name, "CaMCa2C",
                             generate_name(i,j+1, k, l))
                add_reaction(my_rxn_file, my_specie_name, "CaMCa4",
                             generate_name(i,j, k+1, l))
            else:
                pass   
            if not i:
                pass 
            else:
                add_reaction(my_rxn_file, my_specie_name, "2CaC",
                             generate_name(i-1, j+1, k, l))
            if not j:    
                pass  
            else:
                add_reaction(my_rxn_file, my_specie_name, "2CaN",
                             generate_name(i, j-1, k+1, l))
           
            new_name = generate_name(i, j, k, l + 1)
            if "CaM" in my_specie_name:
                rxn_name = "CaMRyR2Ca%d" % (l+1)
            else:
                rxn_name = "RyR2Ca%d" % (l+1)
            add_reaction(my_rxn_file, my_specie_name, rxn_name,
                         new_name)

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

    for i, specie in enumerate(ryr_species_to_open):
        if "CaM" not in specie:
            add_reaction(my_rxn_file, specie, "RyR2Ca4O1",
                         "%s_O1"%specie)
            add_reaction(my_rxn_file, specie, "RyR2Ca4O2",
                         "%s_O2"%specie)
            add_reaction(my_rxn_file, "%s_O1" % specie,"RyR2Ca4O1C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_O2" % specie,"RyR2Ca4O2C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_C1" % specie,"RyR2Ca4C1I",
                         "%s_I"%specie)
        else:
            add_reaction(my_rxn_file, specie, "CaMRyR2Ca4O1",
                         "%s_O1"%specie)
            add_reaction(my_rxn_file, specie, "CaMRyR2Ca4O2",
                         "%s_O2"%specie)
            add_reaction(my_rxn_file, "%s_O1" % specie,"CaMRyR2Ca4O1C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_O2" % specie,"CaMRyR2Ca4O2C1",
                         "%s_C1"%specie)
            add_reaction(my_rxn_file, "%s_C1" % specie,"CaMRyR2Ca4C1I",
                         "%s_I"%specie)

                
    f = open(fname, "w")
    f.write(etree.tostring(my_rxn_file, pretty_print=True).decode("utf-8"))
    f.close()
    print(ryr_species_to_open)
