import sys
import subprocess
from datetime import date
from lxml import etree
import h5py
import numpy as np
import matplotlib.pyplot as plt

ca_conc_file = "ryr_cam_po_no_cam.csv"
ryr_op_fname = "ryr_cam_to_no_cam.csv"
ryr_cl_fname = "ryr_cam_tc_no_cam.csv"

model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="Rxn_Zahradnikova.xml" />
    <xi:include href="Morph.xml" />
    <xi:include href="%s" />
    <xi:include href="IO_RyR_simple_no_cam.xml"/>
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
    <runtime>20000</runtime>

    <!-- set the seed to get the same spines each time testing -->
    <spineSeed>123</spineSeed>

    <discretization>
      <defaultMaxElementSide>10</defaultMaxElementSide>
      <surfaceLayers>10</surfaceLayers> 
    </discretization>
    <tolerance>0.01</tolerance>

    <outputInterval>1000</outputInterval>

    <calculation>GRID_ADAPTIVE</calculation>

</SDRun>"""



IC_text = """<?xml version="1.0" encoding="utf-8"?>
<InitialConditions>
  <ConcentrationSet>
    <NanoMolarity specieID="Ca" value="%f"/>
    <NanoMolarity specieID="RyR4"      value="0.1"    />
  </ConcentrationSet>
</InitialConditions>
"""



def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')

def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))

def region_volumes(my_file):
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])
    return volumes

def get_grid_list(My_file):
    return np.array(My_file['model']['grid'])

def get_times(My_file, trial='trial0', output="__main__"):
    return np.array(My_file[trial]['output'][output]['times'])

def get_outputs(my_file):
    return my_file['model']['output'].keys()

def get_populations(my_file, trial='trial0', output='__main__'):
    return np.array(my_file[trial]['output'][output]['population'])

def get_all_species(My_file, output="__main__"):
    return [s.decode('utf-8') for s in My_file['model']['output'][output]['species']]

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs

def sum_volume(my_file, region_list):
    grid_list = get_grid_list(my_file)
    vol_sum = 0
    volumes = region_volumes(my_file)
    for region in region_list:
        if region in volumes:
            vol_sum += volumes[region]
    return vol_sum

def get_all_open(data, species):
    length = len(data[:, 0, 0])//2
    state = np.zeros(length)
    for specie in species:
        if "O" in specie:
            state +=  data[length+1:, 0, species.index(specie)]
    count = len(np.where((state[1:] - state[0:-1])==1)[0])
    sum_times = state.sum()
    print(sum_times)
    return sum_times, count


def get_numbers(my_file, output="all"):
    output_dict = get_output_regions(my_file)
    Ca_conc = []
    open_ryr3 = []
    no = 0
    sum_o = 0
    nc = 0
    sum_c = 0
 
    for trial in my_file.keys():
 
        if trial == "model":
            continue
        times = get_times(my_file, trial=trial, output=output)
        vol = sum_volume(my_file, ["dend"])
        species = get_all_species(my_file, output=output)
        data = get_populations(my_file, trial=trial, output=output)
        dt = times[1]-times[0]
        exp_len = int((times[-1]))//2
        mean_ca = data[:, 0, species.index("Ca")].mean()*10/6.023/vol
        bas_idx = species.index("RyR4")
        ryr_basal = data[0, 0, bas_idx]
        if ryr_basal > 1:
            continue
        open_sum, tot_no = get_all_open(data, species)
        p_open_ryr = dt*open_sum/ryr_basal/exp_len
        Ca_conc.append(mean_ca)
        open_ryr3.append(p_open_ryr)
        sum_o += dt*open_sum
        if tot_no > 0:
            no += tot_no
            nc += tot_no+1
        else:
            nc += 1
            no += 1
        sum_c += (exp_len - dt*open_sum)
    mean_c_t = sum_c/nc
    mean_o_t = sum_o/no
    return Ca_conc, open_ryr3, mean_o_t, mean_c_t
        

if __name__ == "__main__":
  
        
    exp_res = np.loadtxt(ca_conc_file, skiprows=1, delimiter=',')
    ca_conc_list = exp_res[:, 0]
    output = np.zeros(exp_res.shape)
    mean_times = []
    for i, ca_conc in enumerate(ca_conc_list):
        ca_conc_nM = int(np.ceil(ca_conc))
        IC_name = "Ca_%d_simple_no_cam.xml" % ca_conc_nM
        model_name = "RyR_simple_Ca_%d_no_cam.xml" % ca_conc_nM
        output_name = "RyR_simple_Ca_%d_no_cam.h5" % ca_conc_nM
        with open(IC_name, "w") as fic:
            fic.write(IC_text % ca_conc_nM)
        with open(model_name, "w") as fm:
            fm.write(model_text %  IC_name)

        process = subprocess.run(["/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
                                  "-jar",
                                  "/home/jszmek/new_neurord/neurord-3.2.3-all-deps.jar",
                                  "-Dneurord.trials=20", model_name],
                                 capture_output=True)
        print(process.returncode)
        if not process.returncode:
            my_file = h5py.File(output_name, 'r')
            conc, po, t_o, t_c = get_numbers(my_file, output="all")
            output[i, 0] = np.mean(conc)
            output[i, 1] = np.mean(po)

            mean_times.append([np.mean(conc), t_o, t_c])
            print(output[i], t_o, t_c)
            

            
    exp_open = np.loadtxt(ryr_op_fname, skiprows=1, delimiter=",")
    exp_closed = np.loadtxt(ryr_cl_fname, skiprows=1, delimiter=",")
    mean_times_a = np.array(mean_times)
    date = date.today().strftime("%y%m%d")
    res_fname1 = "po_res_%s.csv" %  date
    res_fname2 = "open_closed_times_res_%s.csv" % date
    np.savetxt(res_fname1, output, delimiter=",", header="Ca [nM], po")
    np.savetxt(res_fname2, mean_times_a, delimiter=",",
               header="Ca [nM], mean open time, mean closed time")
    fig, ax = plt.subplots(1)
    ax.set_xscale('log')
    ax.plot(exp_res[:, 0]*1e-9, exp_res[:, 1], "d", color="tab:blue", label="experimental data")
    ax.plot(output[:, 0]*1e-9, output[:, 1], "d", color="tab:red", label="model data")
    ax.legend()
    ax.set_xlabel("Concentration [M]")
    ax.set_ylabel("RyR2 open probability")
    
    fig, ax = plt.subplots(1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(exp_open[:, 0]*1e-9, exp_open[:, 1], "d", color="tab:blue",
            label="exp open")
    ax.plot(exp_closed[:, 0]*1e-9, exp_closed[:, 1], "d", color="tab:green",
            label="exp closed")
    ax.plot(mean_times_a[:, 0]*1e-9, mean_times_a[:, 1], "d",
            label="model open", color="tab:cyan")
    ax.plot(mean_times_a[:, 0]*1e-9, mean_times_a[:, 2], "d",
            label="model closed", color="tab:olive")
    
    ax.legend()
    ax.set_xlabel("Concentration [M]")
    ax.set_ylabel("Time [ms]")


plt.show()
