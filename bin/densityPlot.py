import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as manimation
import pandas as pd
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
import os


def loadFrame(fileName):
    header = pd.read_csv(fileName, delimiter=" ", header=None, nrows=1)
    time = header[1].values[0]
    data = pd.read_csv(fileName, delimiter=" ", header=None, skiprows=1)
    return data[0].values, data[1].values, time

def plotEvolution(path):
    fig = plt.figure(figsize=(12, 8))
    plt.ylim([1e-3, 1e4])
    plt.xlim([1e-3, 1e3])

    files = os.listdir(path)
    files = sorted([f for f in files if f.startswith("frame")], key=lambda f: int(os.path.splitext(f)[0][5:]))

    for file in files:
        print("Processing data: " + file)
        x, y, t = loadFrame(file)
        plt.loglog(x, y, linewidth=0.5, color='b')
        plt.xlabel("Radius [AU]")
        plt.ylabel("Column Density [g / cm^2]")

    plt.savefig("density_evolution.eps")    

def plotForVideo(path):

    FFMpegWriter = manimation.writers['ffmpeg']  
    metadata = dict(title='Photoevaporation', artist='Peter Rodenkirch',
                comment='')
    writer = FFMpegWriter(fps=60, metadata=metadata)



    fig = plt.figure(figsize=(6, 4))
    l, = plt.loglog([], [])
    plt.ylim([0.00001, 200000])
    plt.xlim([0.1, 1000])

    files = os.listdir(path)
    files = sorted([f for f in files if f.startswith("frame")], key=lambda f: int(os.path.splitext(f)[0][5:]))

    with writer.saving(fig, "test.mp4", 300):
        for file in files:
            print("Processing data: " + file)
            x, y , t = loadFrame(file)
            l.set_data(x, y)
            plt.title("Surface Density Evolution at t = " + str(t) + " years.")
            writer.grab_frame()
            
            #plt.xlabel("Radius [AU]")
            #plt.ylabel("Surface Density [g / cm^2]")
            # plt.show()

plotEvolution("./")