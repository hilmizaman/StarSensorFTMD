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


zz = 0
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
#Data bintang tetangga
ID_sorted = neighbor['col1']
dist_12 = neighbor['col2']




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

print(v)

#Metode jaring (mencari jarak ke bintang tengah)
alfa1 = np.zeros(N_stars-1)
#alfa1 nomor 0 berarti jarak 2 ke 1, nomor 2 berarti 3 ke 1, dst
for i in range(N_stars-1):
    alfa1[i] = mpmath.acos(np.vdot(v[0,:],v[i+1,:]))

print(alfa1)

#error maks untuk jarak paling dekat
error1 = 0.001

#menentukan batas atas dan bawah
upper_lim = alfa1[0]+error1
lower_lim = alfa1[0]-error1

error2 = error1/10

def binary_search(arr, target, err):
    left = 0
    right = len(arr) - 1

    while left <= right:
        mid = (left + right) // 2
        
        # Check if target is present at mid
        if np.abs(arr[mid]-target) < err:
            return mid
        
        # If target is greater, ignore left half
        elif arr[mid] < target:
            left = mid + 1
        
        # If target is smaller, ignore right half
        else:
            right = mid - 1

    # If we reach here, the element was not present
    return -1

#mencari batas atas dan bawah 
N_upper = binary_search(dist_12,upper_lim,error2)
N_lower = binary_search(dist_12,lower_lim,error2)

print([N_lower,N_upper])


#placeholder buat lower upper di loop
low=N_lower
up=N_upper

def sort_2d_array(arr,n):
    # Combine elements of the N-th row with their corresponding elements in other rows using zip
    combined_rows = list(zip(*arr))

    # Sort the combined rows based on the N-th row
    sorted_combined_rows = sorted(combined_rows, key=lambda x: x[n])

    # Separate the sorted elements back into rows
    sorted_arr = list(zip(*sorted_combined_rows))

    return sorted_arr

IDN = np.array(ID_sorted[low:up])
dis2 = np.array(neighbor['col2'][low:up])
dis3 = np.array(neighbor['col3'][low:up])
dis4 = np.array(neighbor['col4'][low:up])
dis5 = np.array(neighbor['col5'][low:up])
Dist_345 = np.array([IDN,dis2,dis3,dis4,dis5])
#print(Dist_345)
dis=[]
Starnum=[]
for ii in range(N_stars-2):
    Star = sort_2d_array(Dist_345,n=ii+2)
    #print(Star[0][20:50])
    upper_lim1 = alfa1[ii+1]+error1*(ii+2)*2
    lower_lim1 = alfa1[ii+1]-error1*(ii+2)*2
    #mencari batas atas dan bawah 
    up = binary_search(Star[ii+2],upper_lim1,error1*(ii+2))
    low = binary_search(Star[ii+2],lower_lim1,error1*(ii+2))
    if up==-1 or low==-1:
        print('gagal bro!')
        break 
    print([low,up])
    Dist_345 = np.array([Star[0][low:up],Star[1][low:up],Star[2][low:up],Star[3][low:up],Star[4][low:up]])
    print(Dist_345)

error_tot = 99999
tot = alfa1[0]+alfa1[1]+alfa1[2]+alfa1[3]
ID_fin = []
for num in range(up-low):
    error_totz = np.abs(Dist_345[1][num]+Dist_345[2][num]+Dist_345[3][num]+Dist_345[4][num]-tot)
    #print(error_totz)
    if error_totz < error_tot:
        error_tot = error_totz
        ID_fin = Dist_345[0][num]

print(ID_fin)