import numpy as np
from astropy.io import ascii

CM0 = ascii.read("./program_staridentification/Catalog_mod2.txt")
N = 5
IDNx = []

RA = CM0['col1']
DE = CM0['col2']

def dist(RA1,DE1,RA2,DE2):
    u1=np.zeros(3)
    u2=np.zeros(3)
    u1[0] = np.cos(RA1) * np.cos(DE1)
    u1[1] = np.sin(RA1) * np.cos(DE1)
    u1[2] = np.sin(DE1)
    u2[0] = np.cos(RA2) * np.cos(DE2)
    u2[1] = np.sin(RA2) * np.cos(DE2)
    u2[2] = np.sin(DE2)
    dd=np.arccos(np.vdot(u1,u2))
    return dd

for i in range(len(RA)):
    A=np.zeros((len(RA),1))
    B=np.zeros((len(RA),1))
    Distance = np.hstack((A,B))
    
    for j in range(len(RA)):
        Distance[j, 0] = j + 1
        Distance[j, 1] = dist(RA[i],DE[i],RA[j],DE[j])

    Distance = sorted(Distance, key=lambda x: x[1])


    top_rows = Distance[1:N+1]
    IDN = [row[0] for row in top_rows]
    IDNx.append(IDN)

#print(Distance[:N])

np.savetxt('./program_staridentification/IDNx2.txt',IDNx,fmt='%d',delimiter=' ')