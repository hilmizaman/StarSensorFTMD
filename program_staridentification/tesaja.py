import numpy as np
import mpmath
from astropy.io import ascii

#spesifikasi CCD
l = 3280
w = 2464

f = 0.00304

myu = 1.12*(10**-6)

FOVy = np.rad2deg(2*np.arctan((myu*w/2)/f))
FOVx = np.rad2deg(2*np.arctan((myu*l/2)/f))
print(FOVx)
print(FOVy)

xtot = 2 * mpmath.tan((FOVx * mpmath.pi / 180) / 2) * f;
ytot = 2 * mpmath.tan((FOVy * mpmath.pi / 180) / 2) * f;
xpixel = l / xtot;
ypixel = w / ytot;

zz=0
#input data bintang
file_from='./result_staridentification/StarPos_sorted/WCG'+str(zz+1)+'.txt'
file_to='./result_staridentification/Result/WCG1_binary'+'.txt'

Catnew = ascii.read("./program_staridentification/Catalog_mod2.txt")
Pixstars = ascii.read(file_from)
neighbor = ascii.read('./program_staridentification/Distance_sorted.txt')
IDNx = ascii.read('./program_staridentification/IDNx_new.txt')

#Data koordinat bintang di CCD
xp = Pixstars['col1']
yp = Pixstars['col2']
#Data koordinat bintang di langit
RA = Catnew['col1']
DE = Catnew['col2']
MA = Catnew['col3']

#input berapa bintang yang digunakan
N_stars = 5

#fungsi untuk jarak antar bintang di CCD
def CCD_dist(x,y):
    xplane = (x-l/2)*xtot/l
    yplane = (-y+w/2)*ytot/w
    v = np.zeros(3)
    v[0] = xplane/np.sqrt(xplane**2 + yplane**2 +f**2)
    v[1] = yplane/np.sqrt(xplane**2 + yplane**2 +f**2)
    v[2] = f/np.sqrt(xplane**2 + yplane**2 +f**2)
    return v

#Koordinat bintang di CCD dalam satuan sudut
v = np.zeros((N_stars,3))
for i in range(N_stars):
    v[i,:] = CCD_dist(xp[i],yp[i])

#print(v)

#Metode jaring (mencari jarak ke bintang tengah)
alfa1 = np.zeros(N_stars-1)
#alfa1 nomor 0 berarti jarak 2 ke 1, nomor 2 berarti 3 ke 1, dst
for i in range(N_stars-1):
    alfa1[i] = mpmath.acos(np.vdot(v[0,:],v[i+1,:]))

print(alfa1)

#error maks untuk jarak paling dekat
error1 = 0.002

#menentukan batas atas dan bawah
upper_lim = alfa1[0]+error1
lower_lim = alfa1[0]-error1

print(upper_lim)
print(lower_lim)