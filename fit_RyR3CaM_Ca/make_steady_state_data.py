import os
import glob
from lxml import etree
import numpy as np

t_end = 5000
dt = 1000
RyR_conc = 35
model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="%s" />
    <xi:include href="Morph.xml" />
    <xi:include href="Rxn_RyRCaM.xml" />

    <!--2D means the morphology is interpreted like a flatworm, 3D for
roundworms. The 2D case is good for testing as it is easy to visualize the
results (also, 3D may not work yet...)  -->
   
    <geometry>          2D           </geometry>
    <depth2D>           10          </depth2D>
    <distribution>      BINOMIAL     </distribution>
    <algorithm>         INDEPENDENT  </algorithm>
    <simulationSeed>    245         </simulationSeed>
    <outputQuantity>NUMBER</outputQuantity>

    <!-- run time for the calculation, milliseconds -->
    <runtime>""" + str(t_end) + """</runtime>

    <!-- set the seed to get the same spines each time testing -->
    <spineSeed>123</spineSeed>

    <discretization>
      <defaultMaxElementSide>10</defaultMaxElementSide>
      <surfaceLayers>10</surfaceLayers> 
    </discretization>
    <tolerance>0.01</tolerance>

    <outputInterval>1</outputInterval>

    <calculation>GRID_ADAPTIVE</calculation>

</SDRun>""" 


def save_conc(ca_conc, RyR_conc):
    fname = "IC_%d.xml" % ca_conc
    my_ic = etree.Element("InitialConditions")
    beg = etree.SubElement(my_ic, "ConcentrationSet")
    ryr = etree.SubElement(beg, "NanoMolarity", 
                           specieID="RyR4CaM", value="%f"%RyR_conc)
    ca = etree.SubElement(beg, "NanoMolarity", 
                          specieID="Ca", value="%f"% (ca_conc))
    with open(fname, "w") as f:
        f.write(etree.tostring(my_ic, pretty_print=True).decode("utf-8"))
    return fname



if __name__ == "__main__":
    original_data = np.loadtxt(os.path.join("..", "datasets_for_fitting",
                                            "ryr3_cam_po.csv"),
                               skiprows=1, delimiter=",")
    for data_point in original_data:
        print(data_point)
        ca_conc = int(1e9*data_point[0])
        print(ca_conc)
        open_channels = abs(data_point[1])*RyR_conc
        new_file_name = "combined_ryr3cam_%d.csv" % ca_conc
        new_file = open(new_file_name, "w")
        new_file.write("time, RO\n")
        new_file.write("0, 0\n")
        for time in range(dt, t_end+dt, dt):
            new_file.write("%d, %f\n" % (time, open_channels))
        new_file.close()
        IC_name = "IC_%d.xml" % ca_conc
        print(save_conc(ca_conc, RyR_conc))
        with open("model_%d.xml" % ca_conc, "w") as fm:
            fm.write(model_text %  IC_name)
        
