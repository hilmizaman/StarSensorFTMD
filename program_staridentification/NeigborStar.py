import numpy as np
from astropy.io import ascii

CM0 = ascii.read("./program_staridentification/Catalog_mod2.txt")
FoV = np.deg2rad((62.2**2+48.8**2)**0.5)
N = 8
IDNx = []

RA = CM0['col1']
DE = CM0['col2']

def dist(RA1,DE1,RA2,DE2):
    dd=np.arccos(np.sin(DE1)*np.sin(DE2)+np.cos(DE1)*np.cos(DE2)*np.cos(RA1-RA2))
    return dd

for i in range(len(RA)):
    A=np.zeros((len(RA),1))
    B=np.ones((len(RA),1))*10
    Distance = np.hstack((A,B))
    
    for j in range(len(RA)):
        Distance[j, 0] = j + 1
        Distance[j, 1] = dist(RA[i],DE[i],RA[j],DE[j])
        if Distance[j,1]>FoV:
            Distance[j,0]=0

    Distance = np.array(sorted(Distance, key=lambda x: x[1]))

    top_rows = np.zeros((1,N+1))
    top_rows[0,0]=i+1
    top_rows[0,1:N+1] = Distance[1:N+1, 1]
    IDN = [row for row in top_rows]
    IDNx.append(top_rows)

#print(IDNx)
# Flatten the 3D array before saving
flat_IDNx = [item for sublist in IDNx for item in sublist]

np.savetxt('./program_staridentification/IDNx_6mag_dist.txt',flat_IDNx,fmt='%.9f',delimiter=' ')