def read_file(file_name):
    with open(file_name, 'r') as file:
        data = file.readlines()
    data = [list(map(float, line.split())) for line in data]
    Pr = data[0]
    Tr = [row[0] for row in data[1:]]
    return [data,Tr,Pr]

def data_fuse(f1,f2):
    [c1,t1,p1] = read_file(f1)
    [c2,_,p2] = read_file(f2)
    # Remove Tr from c2 before combining matrices
    c2 = [row[1:] for row in c2]
    # This is necessary for horizontal matrix fusion
    data = [sublist1 + sublist2 for sublist1, sublist2 in zip(c1, c2)]
    # Trim data to rid of axes
    data = [row[1:] for row in data]
    data = data[1:]
    Pr = p1+p2
    return [data, t1, Pr]

def hunt(tr,pr):
    # Both tables use the same axes
    ti = 0
    tl = len(t0)
    pi = 0
    pl = len(p0)
    while tl-ti>1:
        inti = (ti+tl)//2
        if t0[inti]<tr:
            ti = inti
        else:
            tl = inti
    while pl-pi>1:
        inti = (pi+pl)//2
        if p0[inti]<pr:
            pi = inti
        else:
            pl = inti
    # At this point, we have x, y, x1, x2, y1, and y2
    # We still need x1y1, x1y2, x2y1, and x2y2 from both tables.
    if tl >= len(t0):
        tl = -1
        ti = -2
    if pl >= len(p0):
        pl = -1
        ti = -2
    return [[ti,tl],[pi,pl]]

def dub_int(tr,pr,arr):
    """ Double interpolation function.
    
    Parameters
    ----------
    tr, pr : int[]
        the closest values of temperature and pressure that are lower and higher in a 1D array
    arr : int[]
        The corresponding table values to every combination of tr and pr indeces in a 2D array
    """
    [t,p] = hunt(tr,pr)
    x = tr
    y = pr
    x1 = t0[t[0]]
    x2 = t0[t[1]]
    y1 = p0[p[0]]
    y2 = p0[p[1]]
    x1y1 = arr[t[0]][p[0]]
    x1y2 = arr[t[0]][p[1]]
    x2y1 = arr[t[1]][p[0]]
    x2y2 = arr[t[1]][p[1]]
    term1 = (x2-x)*(y2-y)/((x2-x1)*(y2-y1))*x1y1
    term2 = (x-x1)*(y2-y)/((x2-x1)*(y2-y1))*x2y1
    term3 = (x2-x)*(y-y1)/((x2-x1)*(y2-y1))*x1y2
    term4 = (x-x1)*(y-y1)/((x2-x1)*(y2-y1))*x2y2
    return term1+term2+term3+term4

global c0, c1, t0, p0
[c0,t0,p0]= data_fuse('data/C.7.0.txt','data/C.7.1.txt')
[c1,_,_]= data_fuse('data/C.8.0.txt','data/C.8.1.txt')

t_actual = 500+273
p_actual = 150

omega = [0.008,0.099,0.25,0.309,0.213,0.049]
tc = [190.6, 305.4, 405.6, 508.1, 553.4, 132.9]
pc = [46, 48.74, 112.77, 47.01, 40.73, 34.96]

tr = [t_actual/x for x in tc]
pr = [p_actual/x for x in pc]

with open ('data/tr.txt','w') as file:
    file.write(' '.join(map(str, t0)))

with open ('data/pr.txt','w') as file:
    file.write(' '.join(map(str, p0)))

with open ('data/c0.txt','w') as file:
    for row in c0:
        line = ' '.join(map(str,row))
        file.write(line+'\n')

with open ('data/c1.txt','w') as file:
    for row in c1:
        line = ' '.join(map(str,row))
        file.write(line+'\n')

phi0 = [dub_int(tr[x],pr[x],c0) for x in range(len(tr))]
phi1 = [dub_int(tr[x],pr[x],c1) for x in range(len(tr))]
logphi = [phi0[x]+omega[x]*phi1[x] for x in range(len(tr))]
phi = [10**x for x in logphi]
f = [p_actual*x for x in phi]

