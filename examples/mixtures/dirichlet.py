import numpy as np
#import matplotlib.pyplot as plt

#s= np.random.dirichlet((1, 1, 1, 1, 1), size = 20)
#s = np.random.dirichlet((2, 2, 2, 2, 2), 25)
s = np.random.dirichlet((2, 2, 2, 2, 2), 10)
np.savetxt('molar_fractions10.out', s)
#np.savetxt('molar_fractions25.out', s)
#print(s)
#print(s.shape)

#s = np.random.dirichlet((2, 2, 2, 2, 2), 25).transpose()

#plt.barh(range(25), s[0])
#plt.barh(range(25), s[1], left=s[0], color='g')
#plt.barh(range(25), s[2], left=s[0]+s[1], color='r')
#plt.barh(range(25), s[3], left=s[0]+s[1]+s[2], color='b')
#plt.barh(range(25), s[4], left=s[0]+s[1]+s[2]+s[3], color='m')
#plt.title("Lengths of Strings")
#plt.show()
