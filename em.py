import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from numpy.linalg import inv
from mpl_toolkits.mplot3d import Axes3D

b = 1   #square plate of side 1m
n = 10  #Discretization order
a = b / (n - 1)
k1 = 4 * 3.14159265 * 8.85e-12
x = []
y = []
m = []
for i in range(0, n - 1):
    x.append(((2 * i + 1)* a)/2)
    y.append(((2 * i + 1)* a)/2)

for i in range(0, n - 1):
    for j in range(0, n - 1):
        m.append([x[i], y[j]])
#print(m)
v = []
for i in range(0, (n-1)**2):
    f = lambda y, x : 1/np.sqrt((x - m[i][0])**2 + (y - m[i][1])**2)
    for k in range(0, (n - 1)**2):
        if(k != i):
            c = integrate.dblquad(f, m[k][0]-(a/2), m[k][0]+(a/2), lambda x:m[k][1]-(a/2),lambda x:m[k][1]+(a/2))
            v.append(c[0])
        if(k == i):
            v.append(8 * np.log(1+np.sqrt(2)))
# print(v)
z=np.empty((n-1)**2)
z.fill(k1)
v1=np.reshape(v,(-1,(n-1)**2))
# print(v1)
vinv=inv(v1)
#print(vinv)
sigma=np.matmul(vinv,z)
#print("Charge densities")
#print(sigma)
#print("Charges on each square")
#print(sigma*a*a)
res=sigma*a*a
res2=[]
x2=[]
y2=[]
for i in range(0,(n-1)**2):
    res2.append(res[i])
for i in range(0,(n-1)**2):
    x2.append(m[i][0])
    y2.append(m[i][1])
#print(x2)
#print(y2)
print(res2)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x2,y2,res2)
plt.show()
print("Capacitance of single plate")
print(np.sum(res2))
#Electric field at z=1 metres above the plate in horizontal plane
z=1
def electric(q,r,x,y):
    deno=k1*(((x-r[0])**2+(y-r[1])**2+z**2)**(1.5))
    return q*(x-r[0])/deno,q*(y-r[1])/deno
nx, ny = 30*(n-1), 30*(n-1)
x = np.linspace(-10,20,num=nx)
y = np.linspace(-10,20,num=ny)
X, Y = np.meshgrid(x, y)
Ex, Ey = np.zeros((ny, nx)), np.zeros((ny, nx))
i=0
for i in range(0,(n-1)**2):
    ex, ey = electric(res2[i],m[i],x=X, y=Y)
    Ex += ex
    Ey += ey
fig=plt.figure()
ax=fig.add_subplot(111)
color = 2 * np.log(np.hypot(Ex, Ey))
ax.streamplot(x, y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='->', arrowsize=1)
ax.set_aspect('equal')
plt.show()

