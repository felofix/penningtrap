import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

halloween = ["#ffa50b", "#be0bff"]

Single1pos = np.loadtxt("datafiles/pos_runge_part_0Single1.txt")
Single1poseuler = np.loadtxt("datafiles/pos_euler_part_0Single1.txt")
Single2pos = np.loadtxt("datafiles/pos_runge_part_0Single2.txt")
Single1vel = np.loadtxt("datafiles/vel_runge_part_0Single1.txt")
Single2vel = np.loadtxt("datafiles/vel_runge_part_0Single2.txt")
Both1pos = np.loadtxt("datafiles/pos_runge_part_0Both.txt")
Both2pos = np.loadtxt("datafiles/pos_runge_part_1Both.txt")
Both1poseuler = np.loadtxt("datafiles/pos_euler_part_0Both.txt")
Both2poseuler = np.loadtxt("datafiles/pos_euler_part_1Both.txt")
Both1vel = np.loadtxt("datafiles/vel_runge_part_0Both.txt")
Both2vel = np.loadtxt("datafiles/vel_runge_part_1Both.txt")

# Position in z direction for single particle.
plt.plot(Single1pos[:, 0], Single1pos[:, 3], label='Runge', color='black')
plt.plot(Single1poseuler[:, 0], Single1poseuler[:, 3], label='Euler', color=halloween[1])
plt.xlabel('Time [t]')
plt.grid()
plt.ylabel('z(t) [$\mu m$]')
plt.legend()
plt.show()


# Position xy without interaction..
plt.plot(Single1pos[:, 1], Single1pos[:, 2], label='Particle 1', color='black')
plt.plot(Single2pos[:, 1], Single2pos[:, 2], label='Particle 2', color=halloween[0])
plt.axis('equal')
plt.grid()
plt.xlabel('x [$\mu m$]')
plt.ylabel('y [$\mu m$]')
plt.legend()
plt.show()


# Position xy with interaction..
plt.plot(Both1pos[:, 1], Both1pos[:, 2], label='Particle 1',  color='black')
plt.plot(Both2pos[:, 1], Both2pos[:, 2], label='Particle 2', color=halloween[0])
plt.grid()
plt.axis('equal')
plt.xlabel('x [$\mu m$]')
plt.ylabel('y [$\mu m$]')
plt.legend()
plt.show()


# Phase for x with and without.
plt.plot(Single1pos[:, 1], Single1vel[:, 1], label='Without interaction', color='black')
plt.plot(Both1pos[:, 1], Both1vel[:, 1], label='With interaction', color=halloween[0])
plt.axis('equal')
plt.grid()
plt.xlabel("$x [\mu m] $")
plt.ylabel("$v_x [\frac{\mum}{\mus}]$")
plt.legend()
plt.show()


# Phase for z with and without.
plt.plot(Single1pos[:, 3], Single1vel[:, 3], label='Without interaction', color='black')
plt.plot(Both1pos[:, 3], Both1vel[:, 3], label='With interaction', color=halloween[0])
plt.axis('equal')
plt.xlabel('z [$\mu m$]')
plt.grid()
plt.ylabel("$v_z [\frac{\mum}{\mus}]$")
plt.legend()
plt.show()


# 3D plot without.
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(Single1pos[:, 1], Single1pos[:, 2], Single1pos[:, 3], color='black', label='Particle 1')
ax.plot3D(Single2pos[:, 1], Single2pos[:, 2], Single2pos[:, 3], color=halloween[1], label='Particle 1')
ax.legend()
plt.show()

# 3D plot with.
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(Both1pos[:, 1], Both1pos[:, 2], Both1pos[:, 3], color='black', label='Particle 1')
ax.plot3D(Both2pos[:, 1], Both2pos[:, 2], Both2pos[:, 3], color=halloween[0], label='Particle 1')
ax.legend()
plt.show()


