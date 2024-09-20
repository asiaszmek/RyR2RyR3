import glob
import numpy as np


if __name__ == "__main__":
    header = "#time, Ca, CaM, CaMCa2C, CaMCa2N, CaMCa4"
    file_list = glob.glob("model_old_reduced_*_trial0_dend.txt")
    for fname in file_list:
        print(fname)
        new_fname = fname[:-16]+".csv"
        print(new_fname)
        data = np.loadtxt(fname, skiprows=1, delimiter=" ")
        np.savetxt(new_fname, data, header=header, delimiter=',')
