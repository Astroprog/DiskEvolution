import numpy as np
import matplotlib.pyplot as plt

frame = 0
for i in range(0,1000):
	if i % 1 == 0:
		print("Plotting frame " + str(i));
		plt.plotfile("frame" + str(i) + ".dat", delimiter=' ', cols=(0, 1), names=('R [AU]', 'Surface Density [g/cm^2]'), marker='.', markersize='1.0')
		plt.loglog()
		plt.ylim([0.001,10000])
		plt.xlim([0.1,1000])
		plt.savefig("frame" + str(frame) + ".png")
		plt.close()
		frame = frame + 1