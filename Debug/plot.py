import numpy as np
import matplotlib.pyplot as plt


# frame = 0
# for i in range(0,1000):
#   if i % 1 == 0:
#       print("Plotting frame " + str(i));
#       plt.plotfile("frame" + str(i) + ".dat", delimiter=' ', cols=(0, 1), names=('R [AU]', 'Surface Density [g/cm^2]'), marker='.', markersize='1.0')
#       plt.loglog()
#       plt.ylim([0.001,10000])
#       plt.xlim([0.1,1000])
#       plt.savefig("frame" + str(frame) + ".png")
#       plt.close()
#       frame = frame + 1

i=0
# for fname in ("frame0.dat", "frame300.dat", "frame500.dat", "frame600.dat", "frame800.dat", "frame900.dat", "frame1000.dat", "frame1030.dat", "frame1040.dat", "frame1300.dat",):
#     legendTitle = " "
#     if i == 0:
#       legendTitle = "0 Myr"
#     elif i == 1:
#       legendTitle = "3 Myr"
#     elif i == 2:
#       legendTitle = "5 Myr"
#     elif i == 3:
#       legendTitle = "6 Myr"
#     elif i == 4:
#       legendTitle = "7.1 Myr"
#     elif i == 5:
#       legendTitle = "7.2 Myr"
#     elif i == 6:
#       legendTitle = "7.3 Myr"
#     elif i == 7:
#       legendTitle = "7.33 Myr"
#     elif i == 8:
#       legendTitle = "7.34 Myr"
#     elif i == 9:
#       legendTitle = "7.6 Myr"
        
#     plt.plotfile(fname, delimiter=' ', cols=(0, 1), skiprows=1, newfig=False, names=('R [AU]', 'Surface Density [g/cm^2]'), marker='.', markersize='1.0', label=legendTitle)
#     plt.legend()
#     i = i+1

# plt.loglog()
# plt.ylim([0.001,10000])
# plt.xlim([0.1,10000])
# plt.savefig('figure.png')


for fname in ("frame0.dat", "frame100.dat", "frame200.dat", "frame300.dat", "frame400.dat", "frame410.dat", "frame415.dat", "frame500.dat", "frame550.dat", "frame600.dat", ):
    legendTitle = " "
    if i == 0:
        legendTitle = "0 Myr"
    elif i == 1:
        legendTitle = "0.1 Myr"
    elif i == 2:
        legendTitle = "0.2 Myr"
    elif i == 3:
        legendTitle = "0.3 Myr"
    elif i == 4:
        legendTitle = "0.4 Myr"
    elif i == 5:
        legendTitle = "0.41 Myr"
    elif i == 6:
        legendTitle = "0.415 Myr"
    elif i == 7:
        legendTitle = "0.5 Myr"
    elif i == 8:
        legendTitle = "0.55 Myr"
    elif i == 9:
        legendTitle = "0.6 Myr"

    with open(fname) as f:
        data = f.read()

    data = data.split('\n')

    n = 0
    x = []
    y = []
    for row in data:
        if n == 0:
            n = n + 1
        else:
            splitted = row.split(' ')
            if len(row) > 0:
                x.append(splitted[0])
                y.append(splitted[1])

    plt.plot(x, y, label=legendTitle)
    plt.legend(loc=1, prop={'size': 8})
    i = i + 1

plt.loglog()
plt.ylim([0.01, 10000])
plt.xlim([0.1, 10000])
plt.title("Surface Density Evolution")
plt.xlabel("Radius [AU]")
plt.ylabel("Surface Density [g / cm^2]")
plt.savefig('density.eps')
