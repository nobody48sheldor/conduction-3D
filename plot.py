import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
from matplotlib import cm
import sys

style.use('dark_background')

tmax = 60


with open("res.txt", 'r') as res:
  T = res.readlines()

for i in range(len(T)):
  T[i] = float(T[i])

t = np.linspace(0,tmax,len(T))

plt.clf()
plt.ylim([min(T)-0.1,max(T)+0.1])
plt.plot(t,T)
plt.show()
