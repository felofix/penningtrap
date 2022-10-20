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
h = np.array([time/n1, time/n2, time/n3, time/n4])
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

# relative errors for n = 4000.
r1 = np.zeros((4000,3))
r1[:,0] = np.real(xy1);
r1[:,1] = np.imag(xy1);
r1[:,2] = z1;
r1len = np.zeros(4000)
r1absdiff = np.zeros(4000)
r1absdiffeuler = np.zeros(4000)

for i in range(4000):
    r1len[i] = np.dot(r1[i], r1[i])
    r1absdiff[i] = np.dot(r1[i] - posRungen1[i, 1:], r1[i] - posRungen1[i, 1:])
    r1absdiffeuler[i] = np.dot(r1[i] - posEuler1[i, 1:], r1[i] - posEuler1[i, 1:])
    
# relative errors for n = 8000.
r2 = np.zeros((8000,3))
r2[:,0] = np.real(xy2);
r2[:,1] = np.imag(xy2);
r2[:,2] = z2;
r2len = np.zeros(8000)
r2absdiff = np.zeros(8000)
r2absdiffeuler = np.zeros(8000)

for i in range(8000):
    r2len[i] = np.dot(r2[i], r2[i])
    r2absdiff[i] = np.dot(r2[i] - posRungen2[i, 1:], r2[i] - posRungen2[i, 1:])
    r2absdiffeuler[i] = np.dot(r2[i] - posEuler2[i, 1:], r2[i] - posEuler2[i, 1:])
    
# relative errors for n = 16000.
r3 = np.zeros((16000,3))
r3[:,0] = np.real(xy3);
r3[:,1] = np.imag(xy3);
r3[:,2] = z3;
r3len = np.zeros(16000)
r3absdiff = np.zeros(16000)
r3absdiffeuler = np.zeros(16000)

for i in range(16000):
    r3len[i] = np.dot(r3[i], r3[i])
    r3absdiff[i] = np.dot(r3[i] - posRungen3[i, 1:], r3[i] - posRungen3[i, 1:])
    r3absdiffeuler[i] = np.dot(r3[i] - posEuler3[i, 1:], r3[i] - posEuler3[i, 1:])

# relative errors for n = 32000.
r4 = np.zeros((32000,3))
r4[:,0] = np.real(xy4);
r4[:,1] = np.imag(xy4);
r4[:,2] = z4;
r4len = np.zeros(32000)
r4absdiff = np.zeros(32000)
r4absdiffeuler = np.zeros(32000)

for i in range(32000):
    r4len[i] = np.dot(r4[i], r4[i])
    r4absdiff[i] = np.dot(r4[i] - posRungen4[i, 1:], r4[i] - posRungen4[i, 1:])
    r4absdiffeuler[i] = np.dot(r4[i] - posEuler4[i, 1:], r4[i] - posEuler4[i, 1:])

"""# PLotting relative error for runge kutta.
plt.plot(time1, np.log(r1absdiff/r1len), label = 'RK4 with n = 4000')
plt.plot(time2, np.log(r2absdiff/r2len), label = 'RK4 with n = 8000')
plt.plot(time3, np.log(r3absdiff/r3len), label = 'RK4 with n = 16000')
plt.plot(time4, np.log(r4absdiff/r4len), label = 'RK4 with n = 32000')
plt.legend()
plt.show()
"""

"""# PLotting relative error for Euler.
plt.plot(time1, r1absdiffeuler/r1len, label = 'Euler with n = 4000')
plt.plot(time2, r2absdiffeuler/r2len, label = 'Euler with n = 8000')
plt.plot(time3, r3absdiffeuler/r3len, label = 'Euler with n = 16000')
plt.plot(time4, r4absdiffeuler/r4len, label = 'Euler with n = 32000')
plt.legend()
plt.show()
"""

"""
#Max rel error part.
maxerrorrunge = np.array([np.max(r1absdiff), np.max(r2absdiff), np.max(r3absdiff), np.max(r4absdiff)])
maxerroreuler = np.array([np.max(r1absdiffeuler), np.max(r2absdiffeuler), np.max(r3absdiffeuler), np.max(r4absdiffeuler)])
relerrorrunge = 0
relerroreuler = 0
print(maxerrorrunge)
print(maxerroreuler)

for i in range(1, 4):
    relerrorrunge += 1/3*(np.log(maxerrorrunge[i]/maxerrorrunge[i-1])/(np.log(h[i]/h[i-1])))
    relerroreuler += 1/3*(np.log(maxerroreuler[i]/maxerroreuler[i-1])/(np.log(h[i]/h[i-1])))
print(relerrorrunge, relerroreuler)
"""
