import numpy as np
import matplotlib.pyplot as plt

#datafiles
posRungen1 = np.loadtxt("datafiles/pos_runge_part_0n1.txt")
posRungen2 = np.loadtxt("datafiles/pos_runge_part_0n2.txt")
posRungen3 = np.loadtxt("datafiles/pos_runge_part_0n3.txt")
posRungen4 = np.loadtxt("datafiles/pos_runge_part_0n4.txt")
posEuler1 = np.loadtxt("datafiles/pos_euler_part_0n1.txt")
posEuler2 = np.loadtxt("datafiles/pos_euler_part_0n2.txt")
posEuler3 = np.loadtxt("datafiles/pos_euler_part_0n3.txt")
posEuler4 = np.loadtxt("datafiles/pos_euler_part_0n4.txt")


# constants for analytical soulution.
q = 1.0   # Charge of calcium with singlecharge.
T = 9.64852558e1 # Tesla.
m = 40.078 # mass of calcium.
B0 = 1*T  # Magnetic field.
V = 9.64852558e7
V0 = 25e-3*V #
d = 500 # Size of tube.
omega2z = 2*q*V0/(m*d**2)
omega0 = q*B0/m
v0 = 25
x0 = 20
z0 = 20
time = 50
n1 = 4000
n2 = 8000
n3 = 16000
n4 = 32000
omegap = (omega0 + np.sqrt(omega0**2 - 2*omega2z))/2
omegam = (omega0 - np.sqrt(omega0**2 - 2*omega2z))/2
Ap = (v0 + omegam*x0)/(omegam - omegap)
Am = -(v0 + omegap*x0)/(omegam - omegap)

# time arrays.
time1 = np.linspace(0, time, n1) # n = 4000
time2 = np.linspace(0, time, n2) # n = 8000
time3 = np.linspace(0, time, n3) # n = 16000
time4 = np.linspace(0, time, n4) # n = 32000

# zplane values.
z1 = z0*np.cos(np.sqrt(omega2z)*time1)
z2 = z0*np.cos(np.sqrt(omega2z)*time2)
z3 = z0*np.cos(np.sqrt(omega2z)*time3)
z4 = z0*np.cos(np.sqrt(omega2z)*time4)

#XY plane.
xy1 = Ap*np.exp(-1j*omegap*time1) + Am*np.exp(-1j*omegam*time1)
xy2 = Ap*np.exp(-1j*omegap*time2) + Am*np.exp(-1j*omegam*time2)
xy3 = Ap*np.exp(-1j*omegap*time3) + Am*np.exp(-1j*omegam*time3)
xy4 = Ap*np.exp(-1j*omegap*time4) + Am*np.exp(-1j*omegam*time4)

relerror1E = np.log(abs((posRungen1[:, 3] - z1)/z1))
relerror1R = np.log(abs((posEuler1[:, 3] - z1)/z1))

plt.plot(time1, relerror1)
plt.plot(time1, relerror1)
plt.axis('equal')
plt.show()

