from lib2to3.pgen2.token import EQUAL
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

halloween = ["#ffa50b", "#be0bff", "#1BC50A"]
# All datafiles from the Pennnig Trap
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
plt.title(f"Position of a single particle in $z$-direction with RK4 and Forward Euler")
plt.plot(Single1pos[:, 0], Single1pos[:, 3], label='RK4', color='black')
plt.plot(Single1poseuler[:, 0], Single1poseuler[:, 3], label='Forwrd Euler', color=halloween[1])
plt.xlabel('Time [t]')
plt.grid()
plt.ylabel('z(t) [$\mu m$]')
plt.legend()
plt.savefig("Sigle particle z-direction.pdf")
plt.show()

# Position xy without interaction.
plt.title(f"Position of two particles in the $xy$-plane without interaction")
plt.plot(Single1pos[:, 1], Single1pos[:, 2], label='Particle 1', color='black')
plt.plot(Single2pos[:, 1], Single2pos[:, 2], label='Particle 2', color = 'm')
plt.scatter(Single1pos[0, 1], Single1pos[0, 2], color=halloween[0])
plt.scatter(Single1pos[-1, 1], Single1pos[-1, 2], color=halloween[2])
plt.scatter(Single2pos[0, 1], Single2pos[0, 2], color=halloween[0])
plt.scatter(Single2pos[-1, 1], Single2pos[-1, 2], color=halloween[2])
plt.axis('equal')
plt.grid()
plt.xlabel('x [$\mu m$]')
plt.ylabel('y [$\mu m$]')
plt.legend()
plt.savefig("Two particles xy without interaction.pdf")
plt.show()

# Position xy with interaction..
plt.title(f"Position of two particles in the $xy$-plane with interaction")
plt.plot(Both1pos[:, 1], Both1pos[:, 2], label='Particle 1',  color='black')
plt.plot(Both2pos[:, 1], Both2pos[:, 2], label='Particle 2', color = 'm')
plt.scatter(Both1pos[0, 1], Both1pos[0, 2], color=halloween[0])
plt.scatter(Both1pos[-1, 1], Both1pos[-1, 2], color=halloween[2])
plt.scatter(Both2pos[0, 1], Both2pos[0, 2], color=halloween[0])
plt.scatter(Both2pos[-1, 1], Both2pos[-1, 2], color=halloween[2])
plt.axis('equal')
plt.grid()
plt.xlabel('x [$\mu m$]')
plt.ylabel('y [$\mu m$]')
plt.legend()
plt.savefig("Two particles xy wit interaction.pdf")
plt.show()

# Phase for x with and without.
plt.title(f"Phase space plot in $x$-direction for particle 1 with and without interaction")
plt.plot(Single1pos[:, 1], Single1vel[:, 1], label='Without interaction', color='black')
plt.plot(Both1pos[:, 1], Both1vel[:, 1], label='With interaction', color = 'm')
plt.axis('equal')
plt.grid()
plt.xlabel("$x [\mu m] $")
plt.ylabel("$v_x \left[\\frac{\mu m}{\mu s}\\right]$")
plt.legend()
plt.savefig("Phase plot x-direction with and without interaction.pdf")
plt.show()


# Phase for z with and without interaction.
plt.title(f"Phase space plot in $z$-direction for particle 1 with and without interaction")
plt.plot(Single1pos[:, 3], Single1vel[:, 3], label='Without interaction', color='black')
plt.plot(Both1pos[:, 3], Both1vel[:, 3], label='With interaction', color = 'm')
plt.xlabel('z [$\mu m$]')
plt.grid()
plt.ylabel("$v_z \left[\\frac{\mu m}{\mu s}\\right]$")
plt.legend()
plt.savefig("Phase plot z-direction with and without interaction.pdf")
plt.show()

# 3D plot without.
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(Single1pos[:, 1], Single1pos[:, 2], Single1pos[:, 3], color='black', label='Particle 1')
ax.plot3D(Single2pos[:, 1], Single2pos[:, 2], Single2pos[:, 3], color = 'm', label='Particle 2')
ax.scatter(Single1pos[0, 1], Single1pos[0, 2], Single1pos[0, 3], color=halloween[0])
ax.scatter(Single1pos[-1, 1], Single1pos[-1, 2], Single1pos[-1, 3], color=halloween[2])
ax.scatter(Single2pos[0, 1], Single2pos[0, 2], Single2pos[0, 3], color=halloween[0])
ax.scatter(Single2pos[-1, 1], Single2pos[-1, 2], Single2pos[-1, 3], color=halloween[2])
ax.legend()
ax.set_xlabel('y [$\mu m$]')
ax.set_ylabel('x [$\mu m$]')
ax.set_zlabel('z [$\mu m$]')
ax.set_title("3D-plot of two particles without interaction")
#ax.axis("equal")
plt.savefig("3D plot two particles without interaction.pdf")
plt.show()


# 3D plot with.
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(Both1pos[:, 1], Both1pos[:, 2], Both1pos[:, 3], color='black', label='Particle 1')
ax.plot3D(Both2pos[:, 1], Both2pos[:, 2], Both2pos[:, 3], color = 'm', label='Particle 2')
ax.scatter(Both1pos[0, 1], Both1pos[0, 2], Both1pos[0, 3], color=halloween[0])
ax.scatter(Both1pos[-1, 1], Both1pos[-1, 2], Both1pos[-1, 3], color=halloween[2])
ax.scatter(Both2pos[0, 1], Both2pos[0, 2], Both2pos[0, 3], color=halloween[0])
ax.scatter(Both2pos[-1, 1], Both2pos[-1, 2], Both2pos[-1, 3], color=halloween[2])
ax.legend()
ax.set_xlabel('y [$\mu m$]')
ax.set_ylabel('x [$\mu m$]')
ax.set_zlabel('z [$\mu m$]')
ax.set_title("3D-plot of two particles with interaction")
#ax.axis("equal")
plt.savefig("3D plot two particles with interaction.pdf")
plt.show()


