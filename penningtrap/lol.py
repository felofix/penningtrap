import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


Single1pos = np.loadtxt("datafiles/pos_runge_part_0Single1.txt")
Single2pos = np.loadtxt("datafiles/pos_runge_part_0Single2.txt")
Single1vel = np.loadtxt("datafiles/vel_runge_part_0Single1.txt")
Single2vel = np.loadtxt("datafiles/vel_runge_part_0Single2.txt")
Both1pos = np.loadtxt("datafiles/pos_runge_part_0Both.txt")
Both2pos = np.loadtxt("datafiles/pos_runge_part_1Both.txt")
Both1vel = np.loadtxt("datafiles/vel_runge_part_0Both.txt")
Both2vel = np.loadtxt("datafiles/vel_runge_part_1Both.txt")

"""# Position in z direction for single particle.
plt.plot(Single1pos[:, 0], Single1pos[:, 3], label='Particle 1')
plt.xlabel('Time [t]')
plt.ylabel('z(t) [$\mu m$]')
plt.legend()
plt.show()
"""

"""# Position xy without interaction..
plt.plot(Single1pos[:, 1], Single1pos[:, 2], label='Particle 1')
plt.plot(Single2pos[:, 1], Single2pos[:, 2], label='Particle 2')
plt.axis('equal')
plt.legend()
plt.show()
"""

"""# Position xy with interaction..
plt.plot(Both1pos[:, 1], Both1pos[:, 2], label='Particle 1')
plt.plot(Both2pos[:, 1], Both2pos[:, 2], label='Particle 2')
plt.axis('equal')
plt.legend()
plt.show()
"""

# Phase for x with and without.
plt.plot(Single1pos[:, 1], Single1vel[:, 1], label='Particle 1')
plt.plot(Both2pos[:, 1], Bothvel[:, 1], label='Particle 2')
plt.axis('equal')
plt.xlabel("x")
plt.ylabel("$v_x$")
plt.legend()
plt.show()


