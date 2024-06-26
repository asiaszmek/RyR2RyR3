import itertools
from lxml import etree

fname = "Rxn_module_RyR2_CaM.xml"

kfs = {"CaM": 2.4e-8, #2.1e-8 for Kd of 820 nM
       "CaMCa2C": 3.15e-7, "CaMCa4": 3.66e-7, "2CaC": 6e-5,
       "2CaN":  0.1e-2, "RyRCa1": 8.84e-3, "RyRCa2": 2.88e-3,
       "RyRCa3":0.00381, "O":0.43}


krs = {"CaM": 1.73e-5, "CaMCa2C": 3.67e-6, "CaMCa4": 1.28e-6, "2CaC": 9.1e-3,
       "2CaN": 1000e-3,"RyRCa1": 8.84e-2, "RyRCa2": 2.88e-2, "RyRCa3":0.152,
       "O": 7.43}
counter = 1

def add_reaction(root, name, what, new_name):
    global counter
    my_r = etree.SubElement(root, "Reaction",
                            name=name+"_"+what+"_"+str(counter),
                            id=name+"_"+what+"_"+str(counter))
    etree.SubElement(my_r, "Reactant", specieID=name)
    if what != "2CaC" and what != "2CaN" and not what.startswith("RyR") and what != "O":
        etree.SubElement(my_r, "Reactant", specieID=what)
    elif what == "RyRCa1" or what == "RyRCa2":
        etree.SubElement(my_r, "Reactant", specieID="Ca")
    elif what == "O":
        pass
    else:
        etree.SubElement(my_r, "Reactant", specieID="Ca", n="2")
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
        name = "RyR"
    else:
        name = "Ca%d_RyR" % d
    if a:
        name += "_%dCaM" % a
    if b:
        name += "_%dCaMCa2C" % b
    if c:
        name += "_%dCaMCa4" % c
    return name

if __name__ == "__main__":
    states = set()
    for i in range(0, 5):
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


    for l in [0, 1, 2, 3]:
        for (i,  j, k) in sorted(states):
            my_specie_name = generate_name(i, j, k, l)
            etree.SubElement(my_rxn_file, "Specie", name=my_specie_name,
                             id=my_specie_name, kdiff="0", kdiffunit="mu2/s")
            if my_specie_name.startswith("Ca3_") and not "4" in my_specie_name:
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
            if l < 3:
                new_name = generate_name(i, j, k, l + 1)
                rxn_name = "RyRCa%d" % (l+1)
                add_reaction(my_rxn_file, my_specie_name, rxn_name,
                             new_name)

    for i, specie in enumerate(ryr_species_to_open):
        etree.SubElement(my_rxn_file, "Specie", name="O%d" % i,
                         id="O%d" % i, kdiff="0", kdiffunit="mu2/s")

    for i, specie in enumerate(ryr_species_to_open):
        add_reaction(my_rxn_file, specie, "O",
                     "O%d"%i)
        
                
    f = open(fname, "w")
    f.write(etree.tostring(my_rxn_file, pretty_print=True).decode("utf-8"))
    f.close()
    print(ryr_species_to_open)
