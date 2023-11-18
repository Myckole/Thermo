import numpy as np
def ant(a,b,c,t):
    lnp = a-b/(t+c)
    return np.exp(lnp)*750

a = [10.0311, 9.2806]
b = [2940.46, 2788.51]
c = [-35.93, -52.36]
t = 298
for i in [0,1]:
    print(ant(a[i],b[i],c[i],t))
