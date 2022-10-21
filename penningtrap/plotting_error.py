import numpy as np
import matplotlib.pyplot as plt

#datafiles
# abserror
RK1abserror = np.loadtxt("datafiles/abserrorNrunge4000.000000.txt")
RK2abserror = np.loadtxt("datafiles/abserrorNrunge8000.000000.txt")
RK3abserror = np.loadtxt("datafiles/abserrorNrunge16000.000000.txt")
RK4abserror = np.loadtxt("datafiles/abserrorNrunge32000.000000.txt")
EU1abserror = np.loadtxt("datafiles/abserrorNeuler4000.000000.txt")
EU2abserror = np.loadtxt("datafiles/abserrorNeuler8000.000000.txt")
EU3abserror = np.loadtxt("datafiles/abserrorNeuler16000.000000.txt")
EU4abserror = np.loadtxt("datafiles/abserrorNeuler32000.000000.txt")

# relerror
RK1relerror = np.loadtxt("datafiles/relerrorNrunge4000.000000.txt")
RK2relerror = np.loadtxt("datafiles/relerrorNrunge8000.000000.txt")
RK3relerror = np.loadtxt("datafiles/relerrorNrunge16000.000000.txt")
RK4relerror = np.loadtxt("datafiles/relerrorNrunge32000.000000.txt")
EU1relerror = np.loadtxt("datafiles/relerrorNeuler4000.000000.txt")
EU2relerror = np.loadtxt("datafiles/relerrorNeuler8000.000000.txt")
EU3relerror = np.loadtxt("datafiles/relerrorNeuler16000.000000.txt")
EU4relerror = np.loadtxt("datafiles/relerrorNeuler32000.000000.txt")

# constants for analytical soulution.
n1 = 4000
n2 = 8000
n3 = 16000
n4 = 32000
time = 50
h = np.array([time/n1, time/n2, time/n3, time/n4])

"""# PLotting relative error for runge kutta.
plt.plot(RK1abserror[1:, 0], np.log(RK1relerror[1:, 1]), label = 'RK4 with n = 4000')
plt.plot(RK2abserror[1:, 0], np.log(RK2relerror[1:, 1]), label = 'RK4 with n = 8000')
plt.plot(RK3abserror[1:, 0], np.log(RK3relerror[1:, 1]), label = 'RK4 with n = 16000')
plt.plot(RK4abserror[1:, 0], np.log(RK4relerror[1:, 1]), label = 'RK4 with n = 32000')

# PLotting relative error for Euler.
plt.plot(EU1abserror[1:, 0], np.log(EU1relerror[1:, 1]), label = 'RK4 with n = 4000')
plt.plot(EU2abserror[1:, 0], np.log(EU2relerror[1:, 1]), label = 'RK4 with n = 8000')
plt.plot(EU3abserror[1:, 0], np.log(EU3relerror[1:, 1]), label = 'RK4 with n = 16000')
plt.plot(RK4abserror[1:, 0], np.log(EU4relerror[1:, 1]), label = 'RK4 with n = 32000')
plt.legend()
plt.show()
"""

#Max rel error part.
maxerrorrunge = np.array([np.max(RK1abserror[:, 1]), np.max(RK2abserror[:, 1]), np.max(RK3abserror[:, 1]), np.max(RK4abserror[:, 1])])
maxerroreuler = np.array([np.max(EU1abserror[:, 1]), np.max(EU2abserror[:, 1]), np.max(EU3abserror[:, 1]), np.max(EU4abserror[:, 1])])
relerrorrunge = 0
relerroreuler = 0

for i in range(1, 4):
    relerrorrunge += 1/3*(np.log(maxerrorrunge[i]/maxerrorrunge[i-1])/(np.log(h[i]/h[i-1])))
    relerroreuler += 1/3*(np.log(maxerroreuler[i]/maxerroreuler[i-1])/(np.log(h[i]/h[i-1])))
print(relerrorrunge, relerroreuler)


