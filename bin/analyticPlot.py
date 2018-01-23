import numpy as np
import matplotlib.pyplot as plt

fname = "check.dat"

with open(fname) as f:
    data = f.read()

data = data.split('\n')

n = 0
x = []
y = []
for row in data:
    splitted = row.split('\t')
    if len(row) > 0:
        x.append(float(splitted[0]))
        y.append(float(splitted[1]))

plt.plot(x, y)

au = 1.495978707e13
r = 1.18921 * au
M = 1.98855e33
M_D = M * 0.015
R_s = 10.0 * au
R = r / R_s
temp = 280 / np.sqrt(R_s / au)
kb = 1.38064852e-16
mu = 2.3
mp = 1.672621898e-24
c_s = np.sqrt(kb * temp / (mu * mp))
G = 6.67408e-8
omega = np.sqrt(G * M / np.power(R_s, 3))
nu = 0.01 * np.power(c_s, 2) / omega
year = 3.1536e7

t_s = np.power(R_s, 2) / (3 * nu)

t = np.linspace(1e4, 1e7, 1000)
T = t * year / t_s + 1.0

y = M_D / (2 * np.pi * np.power(R_s, 2)) * 1 / (R * np.power(T, 1.5)) \
    * np.exp(-R / T)

plt.loglog(t, y)

# plt.ylim([0.00001, 200000])
# plt.xlim([0.1, 5000])
# plt.title("Disk Lifetime in Relation to Lambda (Disk mass 0.025 m_sol)")
plt.xlabel("Time [yr]")
plt.ylabel("Surface Density")
plt.savefig('viscous_analytic.eps')
