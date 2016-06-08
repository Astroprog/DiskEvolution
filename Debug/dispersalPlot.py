import numpy as np
import matplotlib.pyplot as plt

fname = "dispersal.dat"

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
            y.append(float(splitted[1]) / 1e6)

plt.plot(x, y)

# plt.ylim([0.00001, 200000])
# plt.xlim([0.1, 5000])
plt.title("Disk Lifetime in Relation to Lambda (Disk mass 0.025 m_sol)")
plt.xlabel("Lambda")
plt.ylabel("Lifetime [Myr]")
plt.savefig('dispersal2.eps')
