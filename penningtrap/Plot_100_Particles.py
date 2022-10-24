import numpy as np
import matplotlib.pyplot as plt
N = 100

f01 = np.loadtxt('datafiles/f0.100000.txt')
f04 = np.loadtxt('datafiles/f0.400000.txt')
f07 = np.loadtxt('datafiles/f0.700000.txt')
fwithout = np.loadtxt('datafiles/100partwithout.txt')
fwith = np.loadtxt('datafiles/100partwith.txt')

plt.title("Particles left after 500 $\mu s$ for different frequencies with f=0.1 $\mu m$.")
plt.grid()
plt.plot(f01[:, 0], f01[:, 1]/N,'--bo', color='#BF40BF')
plt.xlabel("$\omega_V [MHz]$")
plt.ylabel("Fraction of particles")
plt.savefig("100part01.pdf")
plt.show()

plt.title("Particles left after 500 $\mu s$ for different frequencies with f=0.4 $\mu m$.")
plt.grid()
plt.plot(f04[:, 0], f04[:, 1]/N,'--bo', color='#9F2B68')
plt.xlabel("$\omega_V [MHz]$")
plt.ylabel("Fraction of particles")
plt.savefig("100part04.pdf")
plt.show()

plt.title("Particles left after 500 $\mu s$ for different frequencies with f=0.7 $\mu m$.")
plt.grid()
plt.plot(f07[:, 0], f07[:, 1]/N, '--bo', color='#5D3FD3')
plt.xlabel("$\omega_V [MHz]$")
plt.ylabel("Fraction of particles")
plt.savefig("100part07.pdf")
plt.show()

plt.title("Particles left in penning trap for frequencies around $\omega_V = 0.7 \mu m.$")
plt.grid()
plt.plot(fwithout[:, 0], fwithout[:, 1]/N, '--bo', color='#5D3FD3', label='Coulumb interactions off. ')
plt.plot(fwith[:, 0], fwith[:, 1]/N, '--bo', color='#9F2B68', label='Coulumb interactions on. ')
plt.xlabel("$\omega_V [MHz]$")
plt.ylabel("Fraction of particles")
plt.legend()
plt.savefig("100withwithout.pdf")
plt.show()


