import os
import glob
from lxml import etree
import numpy as np

t_end = 5000
dt = 1000
RyR_conc = 35

def save_conc(ca_conc, RyR_conc):
    fname = "IC_%f.xml" % ca_conc
    my_ic = etree.Element("InitialConditions")
    beg = etree.SubElement(my_ic, "ConcentrationSet")
    ryr = etree.SubElement(beg, "NanoMolarity", 
                           specieID="RyR", value="%d"%RyR_conc)
    ca = etree.SubElement(beg, "NanoMolarity", 
                          specieID="Ca", value="%f"% (ca_conc*1e9))
    with open(fname, "w") as f:
        f.write(etree.tostring(my_ic, pretty_print=True).decode("utf-8"))
    return fname



if __name__ == "__main__":
    original_data = np.loadtxt("../datasets_for_fitting/ryr2_mg_ca_Copello_et_al_1997.csv", skiprows=1,
                               delimiter=",")
    for data_point in original_data:
        print(data_point)
        ca_conc = round(10**(-data_point[0]), 8)
        open_channels = data_point[1]*RyR_conc
        new_file_name = "copello_ca_%f.csv" % ca_conc
        new_file = open(new_file_name, "w")
        new_file.write("time, RO\n")
        new_file.write("0, 0\n")
        for time in range(dt, t_end+dt, dt):
            new_file.write("%d, %f\n" % (time, open_channels))
        new_file.close()
        print(save_conc(ca_conc, RyR_conc))
        
