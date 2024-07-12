import os
import sys
import subprocess
import h5py

import numpy as np
import matplotlib.pyplot as plt


po_file = "po_%d0_uM_Ca.csv"
model_file = "model_Z_%d0_uM_Ca.xml"
output_file = "model_Z_%d0_uM_Ca.h5"


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

def get_open(my_file, output="__main__"):
    open_conc = []
    for trial in my_file.keys():
        out = []
        if not trial.startswith("trial"):
            continue
        times = get_times(my_file, trial=trial, output=output)
        vol = sum_volume(my_file, ["dend"])
        
        species = get_all_species(my_file, output=output)
        data = get_populations(my_file, trial=trial, output=output)
        specie_idx = [i for i, specie in enumerate(species) if "O" in specie]
        open_conc.append(np.array(data[:, :, specie_idx].sum(axis=(1,2))*10/6.0023/vol))
   
    return times, np.array(open_conc).mean(axis=0)
        

if __name__ == "__main__":
    fig, ax = plt.subplots(1)
    for c in range(1, 4):
        data = np.loadtxt(po_file % c, skiprows=1, delimiter=",")
        process = subprocess.run(["/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
                                  "-jar",
                                  "/home/jszmek/new_neurord/neurord-3.2.3-all-deps.jar",
                                  "-Dneurord.trials=10",  model_file % c],
                                 capture_output=True)
        if not process.returncode:
            print(output_file % c)
            my_file = h5py.File(output_file % c, 'r')
            time, model_open = get_open(my_file)

            ax.plot(data[:, 0], data[:, 1], label="Zahradnikova %d0 uM Ca" %c)
            ax.plot(time, model_open, label="Model %d0 uM Ca" %c)
    ax.legend()
    plt.show()
            
