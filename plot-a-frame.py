import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
from matplotlib import cm
import sys

style.use('dark_background')

n=10

index = int(sys.argv[1])

print("plotting t =", index, "\n")



L = np.linspace(-5, 5, n)

X = []
Y = []
Z = []
V = []
W = []

for z in range(n):
    for y in range(n):
        for x in range(n):
            X.append(L[x])
            Y.append(L[y])
            Z.append(L[z])
            # V.append(function3D(L[x], L[y], L[z]))

#print(len(X))
#print(len(Y))
#print(len(Z))

with open("data/data_0.txt".format(index), 'r') as file:
    data = file.readlines()


data = data[0].split("/",n*n*n)
for i in data:
    try:
        W.append(float(i))
    except:
        None

maximum = max(W)
minimum = min(W)

with open("data/data_{}.txt".format(index), 'r') as file:
    data = file.readlines()


data = data[0].split("/",n*n*n)
for i in data:
    try:
        V.append(float(i))
    except:
        None


#print(len(V))


fig = plt.figure()
plt.clf()
ax = fig.add_subplot(111, projection = '3d')

for i in range(n**3):
    c = (V[i]+abs(minimum))/(maximum+abs(minimum))
    if abs(c) > 1:
        c = 0
    ax.scatter3D(X[i],Y[i],Z[i], color = [c*0.99, 0, 1-(0.99*c)], alpha=0.5*c+0.1, s=50*c+20, linewidth=0, antialiased=True)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.axes.set_xlim3d(left=-5, right=5)
ax.axes.set_ylim3d(bottom=-5, top=5)
ax.axes.set_zlim3d(bottom=-5, top=5)


plt.savefig("renders/render_{}.png".format(index))
#plt.show()