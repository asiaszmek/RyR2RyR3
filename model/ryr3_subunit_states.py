import itertools
import argparse
from lxml import etree


default_fname = "Rxn_module_RyR3_CaM.xml"

parser = argparse.ArgumentParser(description='Generate RyR2-CaM-Ca_cyt reactions')
parser.add_argument('--release', default=False,
                    help='Add Ca release from the ER')
parser.add_argument('--CaM_k_rev_multiplier', type=float, default=1,
                    help='Increase CaM-RyR binding k_rev to simulate oxidative conditions')
parser.add_argument('--output', default=default_fname,
                    help="Destination output")



A = 4

counter = 1

kfs = {"CaM": 2.1e-8, # for Kd of 820 nM (Xu and Meissner 2004 for RyR2)
       "CaMCa2C": 3.15e-7, "CaMCa4": 3.66e-7, "2CaC": 6e-5,
       "2CaN":  0.1e-2, "RyR3Ca1": 10e-3, "RyR3Ca2": 7.5e-3,
       "RyR3Ca3":5e-3, "RyR3Ca4": 2.5e-3, "RyR3Ca4O1": 263.66,
       "RyR3Ca4O1C1": 0.00734,
       "RyR3Ca4O2": 263.66e-2, "RyR3Ca4O2C1":0.00734e2, "RyR3Ca4C1I": 0.48,
       "II2":1.4, "release":1.6667e-3
      }


krs = {"CaM": 1.73e-5, "CaMCa2C": 2.59e-5, "CaMCa4": 3.015e-6,
       "2CaC": 9.1e-3, "2CaN": 1000e-3,"RyR3Ca1": 2.5,
       "RyR3Ca2": 5, "RyR3Ca3": 7.5, "RyR3Ca4": 10,
       "RyR3Ca4O1": 1.026,"RyR3Ca4O1C1": 0.6296, "RyR3Ca4O2": 1.026e-2,
       "RyR3Ca4O2C1": 0.6296e2,  "RyR3Ca4C1I":0.04,
       "II2":0.1525, "release": 5e-3
       }

for i in range(1, 5):
    new_specie = "CaMRyR3Ca%d"%i
    old_specie = "RyR3Ca%d"%i
    kfs[new_specie] = 0.75*kfs[old_specie]
    krs[new_specie] = 0.75*A*krs[old_specie]

open_close = ["RyR3Ca4O1", "RyR3Ca4O1C1",
              "RyR3Ca4O2", "RyR3Ca4O2C1", "RyR3Ca4C1I"]
    
for specie in open_close:
    new_specie = "CaM%s" % specie
    kfs[new_specie] = kfs[specie]
    krs[new_specie] = krs[specie]


def write_rates(root1, k_forward, k_reverse):
    kf = etree.SubElement(root1, "forwardRate")
    kf.text = str(k_forward)
    kr = etree.SubElement(root1, "reverseRate")
    kr.text = str(k_reverse)
    q = etree.SubElement(root1, "Q10")
    q.text = ".2"


def add_reaction(root, name, what, new_name):
    global counter
    multiplier = 1
    my_r = etree.SubElement(root, "Reaction",
                            name=name+"_"+what+"_"+str(counter),
                            id=name+"_"+what+"_"+str(counter))
    etree.SubElement(my_r, "Reactant", specieID=name)
    if what == "2CaC" or what == "2CaN":
        etree.SubElement(my_r, "Reactant", specieID="Ca", n="2")
    elif what in ["CaM", "CaMCa2C", "CaMCa4"]:
        etree.SubElement(my_r, "Reactant", specieID=what)
        if name.startswith("Ca"):
            Ca_no = int(name.split("_")[0][-1])
            if "CaM" not in name:
                CaM_no = 0
                if CaM_no == 0:
                    multiplier= A**(Ca_no-CaM_no)
  
    elif "O1" in what or "O2" in what or "C1" in what or "I" in what:
        pass
    elif what == "release":
        etree.SubElement(my_r, "Reactant", specieID="CaER")
        etree.SubElement(my_r, "Product", specieID="Ca")
    else:
        etree.SubElement(my_r, "Reactant", specieID="Ca")
    etree.SubElement(my_r, "Product", specieID=new_name)
    write_rates(my_r, kfs[what], krs[what]*multiplier)
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
    args = parser.parse_args()
    fname = args.output
    if args.CaM_k_rev_multiplier > 1:
        krs["CaM"] = args.CaM_k_rev_multiplier*krs["CaM"]
        krs["CaMCa2C"] = args.CaM_k_rev_multiplier*krs["CaMCa2C"]
        krs["CaMCa4"] = args.CaM_k_rev_multiplier*krs["CaMCa4"]

  
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
        if "CaMCa4" in specie:
            etree.SubElement(my_rxn_file, "Specie", name="%s_I2" % specie,
                             id="%s_I2" % specie, kdiff="0",
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
            if "CaMCa4" in specie:
                add_reaction(my_rxn_file, "%s_O1" % specie,"II2",
                             "%s_I2"%specie)
                add_reaction(my_rxn_file, "%s_O2" % specie,"II2",
                             "%s_I2"%specie)
    if args.release:
        etree.SubElement(my_rxn_file, "Specie", name="CaER",
                         id="CaER", kdiff="10", kdiffunit="mu2/s")
        for specie in ryr_species_to_open:
             add_reaction(my_rxn_file, "%s_O1" % specie, "release",
                          "%s_O1" % specie) 
             add_reaction(my_rxn_file, "%s_O2" % specie, "release",
                          "%s_O2" % specie) 
                
    f = open(fname, "w")
    f.write(etree.tostring(my_rxn_file, pretty_print=True).decode("utf-8"))
    f.close()

