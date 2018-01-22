import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f")
args = parser.parse_args()

frame = "frame" + args.f + ".dat"

data = pd.read_csv(frame, delimiter=" ", header=None, skiprows=1)

plt.figure(figsize=(20, 15))
plt.loglog(data[0].values, data[1].values)
plt.ylim([0.00001, 200000])
plt.xlim([0.1, 1000])
plt.title("Surface Density Evolution with Standard Photoevaporation")
plt.xlabel("Radius [AU]")
plt.ylabel("Surface Density [g / cm^2]")
plt.show()