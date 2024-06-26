import itertools
from lxml import etree

kfs = {"CaM": 2.4e-8, "CaMCa2C": 3.15e-7, "CaMCa4": 3.66e-7, "2CaC": 6e-5,
       "2CaN":  0.1e-2}


krs = {"CaM": 1.73e-5, "CaMCa2C": 3.67e-6, "CaMCa4": 1.28e-6, "2CaC": 9.1e-3,
       "2CaN": 1000e-3}
counter = 1

def add_reaction(root, name, what, new_name):
    global counter
    my_r = etree.SubElement(root, "Reaction",
                            name=name+"_"+what+"_"+str(counter),
                            id=name+"_"+what+"_"+str(counter))
    etree.SubElement(my_r, "Reactant", specieID=name)
    if what != "2CaC" and what != "2CaN":
        etree.SubElement(my_r, "Reactant", specieID=what)
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

    
def generate_name(a, b, c):
    name = "RyR"
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
                    for iteration in itertools.permutations((i,j,k),r=3):
                        states.add(iteration)

    my_rxn_file = etree.Element("ReactionScheme")
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



    
    for (i,  j, k) in sorted(states):
        my_specie_name = generate_name(i, j, k)
        etree.SubElement(my_rxn_file, "Specie", name=my_specie_name,
                         id=my_specie_name, kdiff="0", kdiffunit="mu2/s")

    
    for (i,  j, k) in sorted(states):
        my_specie_name = generate_name(i, j, k)
    
        if i+j+k < 4:
            add_reaction(my_rxn_file, my_specie_name, "CaM",
                         generate_name(i+1,j, k))
            add_reaction(my_rxn_file, my_specie_name, "CaMCa2C",
                         generate_name(i,j+1, k))
            add_reaction(my_rxn_file, my_specie_name, "CaMCa4",
                         generate_name(i,j, k+1))
        else:
            pass   
        if not i:
            pass 
        else:
            add_reaction(my_rxn_file, my_specie_name, "2CaC",
                         generate_name(i-1, j+1, k))
        if not j:    
            pass  
        else:
            add_reaction(my_rxn_file, my_specie_name, "2CaN",
                         generate_name(i, j-1, k+1))
    print(etree.tostring(my_rxn_file, pretty_print=True).decode("utf-8"))
