import os
import glob
import numpy as np

directory = os.path.join("..", "datasets_for_fitting")

files = glob.glob(os.path.join(directory, "po_?0_uM_Ca.csv"))

for fname in files:
    data = np.loadtxt(fname, skiprows=1, delimiter=',')
    data[np.where(data[:, 1] < 0)] = 0
    data[:, 1] = data[:, 1]*35
    new_fname = "po" + fname.split("po")[-1]
    np.savetxt(new_fname, data, header="time,RO", delimiter=",")
