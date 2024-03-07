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

OUT=[]
for zz in range(100):
    #input data bintang
    file_from='./result_staridentification/StarPos_sorted/WCG'+str(zz+1)+'.txt'
    file_to='./result_staridentification/Result/WCG_bin'+str(zz+1)+'.txt'

    Catnew = ascii.read("./program_staridentification/Catalog_mod2.txt")
    Pixstars = ascii.read(file_from)
    neighbor = ascii.read('./program_staridentification/Distance_sorted.txt')
    IDNx = ascii.read('./program_staridentification/IDNx2.txt')

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
    ID2 = IDNx['col1']
    ID3 = IDNx['col2']




    #input berapa bintang yang digunakan
    N_stars = 5

    #fungsi untuk jarak antar bintang di CCD
    def CCD_dist(x,y):
        xplane = (x)*xtot/l
        yplane = (y)*ytot/w
        v = np.zeros(3)
        v[0] = xplane/np.sqrt(xplane**2 + yplane**2 +f**2)
        v[1] = yplane/np.sqrt(xplane**2 + yplane**2 +f**2)
        v[2] = -f/np.sqrt(xplane**2 + yplane**2 +f**2)
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
    error1 = 0.00025
    error2 = error1/10
    if alfa1[0]<0.02:
        error1 = 0.0005
        error2 = error1/10
    if alfa1[0]>0.035:
        error1 = 0.0005
        error2 = error1/10
    if alfa1[0]>0.06:
        error1 = 0.001
        error2 = error1/2

    #menentukan batas atas dan bawah
    upper_lim = alfa1[0]+error1
    lower_lim = alfa1[0]-error1



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
        upper_lim1 = alfa1[ii+1]+error1*(ii+3)*2
        lower_lim1 = alfa1[ii+1]-error1*(ii+3)*2
        #mencari batas atas dan bawah 
        uppa = binary_search(Star[ii+2],upper_lim1,error1*(ii+3))
        lowwa = binary_search(Star[ii+2],lower_lim1,error1*(ii+3))
        if uppa==-1 or lowwa==-1:
            print('iter no'+str(zz+1)+'gagal bro pas ii='+str(ii))
            break 
        #print([low,up])
        up = uppa
        low = lowwa
        Dist_345 = np.array([Star[0][low:up],Star[1][low:up],Star[2][low:up],Star[3][low:up],Star[4][low:up]])
        #print(Dist_345)

    error_tot = 99999
    tot = alfa1[0]+alfa1[1]+alfa1[2]+alfa1[3]
    #print([up,low])
    ID_fin = []
    for num in range(up-low):
        error_totz = np.abs(Dist_345[1][num]+Dist_345[2][num]+Dist_345[3][num]+Dist_345[4][num]-tot)
        #print(error_totz)
        if error_totz < error_tot:
            error_tot = error_totz
            ID_fin = Dist_345[0][num]

    #print(Dist_345)
    #print(ID_fin)
    #print(IDNx['col1'][int(ID_fin-1)])
    ID_fin = int(ID_fin)
    ID_fin2 = int(IDNx['col1'][int(ID_fin-1)])
    ID_fin3 = int(IDNx['col2'][int(ID_fin-1)])
    #print([ID_fin,ID_fin2,ID_fin3])


    #mencari koordinat di langit
    xsky1 = mpmath.cos(RA[ID_fin-1])*mpmath.cos(DE[ID_fin-1]);
    ysky1 = mpmath.sin(RA[ID_fin-1])*mpmath.cos(DE[ID_fin-1]);
    zsky1 = mpmath.sin(DE[ID_fin-1]);

    xsky2 = mpmath.cos(RA[ID_fin2-1])*mpmath.cos(DE[ID_fin2-1]);
    ysky2 = mpmath.sin(RA[ID_fin2-1])*mpmath.cos(DE[ID_fin2-1]);
    zsky2 = mpmath.sin(DE[ID_fin2-1]);

    xsky3 = mpmath.cos(RA[ID_fin3-1])*mpmath.cos(DE[ID_fin3-1]);
    ysky3 = mpmath.sin(RA[ID_fin3-1])*mpmath.cos(DE[ID_fin3-1]);
    zsky3 = mpmath.sin(DE[ID_fin3-1]);

    print([np.rad2deg(RA[ID_fin-1]),np.rad2deg(DE[ID_fin-1])]) #sudah bener
    OUT.append({'ID':ID_fin,'RA':np.rad2deg(RA[ID_fin-1]),'dec':np.rad2deg(DE[ID_fin-1])})


import csv
# Specify the file name
file_name = 'Result_sementara.csv'

# Open the file in write mode
with open(file_name, 'w', newline='') as csvfile:
    # Define the header (column names)
    fieldnames = ['ID','RA', 'dec']
    
    # Create a CSV writer
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Write the rows
    writer.writerows(OUT)

print(f'Table has been written to {file_name}')


'''
vv1 = np.zeros(3)
vv2 = np.zeros(3)
vv3 = np.zeros(3)

vv1[0] = xsky1
vv1[1] = ysky1
vv1[2] = zsky1

vv2[0] = xsky2
vv2[1] = ysky2
vv2[2] = zsky2

vv3[0] = xsky3
vv3[1] = ysky3
vv3[2] = zsky3

print(vv1)
print(vv2)
print(vv3)

B1 = np.zeros((3,3))
B2 = np.zeros((3,3))
B3 = np.zeros((3,3))

B1[0][0] = xsky1*v[0,0]
B1[0][1] = xsky1*v[0,1]
B1[0][2] = xsky1*v[0,2]
B1[1][0] = ysky1*v[0,0]
B1[1][1] = ysky1*v[0,1]
B1[1][2] = ysky1*v[0,2]
B1[2][0] = zsky1*v[0,0]
B1[2][1] = zsky1*v[0,1]
B1[2][2] = zsky1*v[0,2]

B2[0][0] = xsky2*v[1,0]
B2[0][1] = xsky2*v[1,1]
B2[0][2] = xsky2*v[1,2]
B2[1][0] = ysky2*v[1,0]
B2[1][1] = ysky2*v[1,1]
B2[1][2] = ysky2*v[1,2]
B2[2][0] = zsky2*v[1,0]
B2[2][1] = zsky2*v[1,1]
B2[2][2] = zsky2*v[1,2]

B3[0][0] = xsky3*v[2,0]
B3[0][1] = xsky3*v[2,1]
B3[0][2] = xsky3*v[2,2]
B3[1][0] = ysky3*v[2,0]
B3[1][1] = ysky3*v[2,1]
B3[1][2] = ysky3*v[2,2]
B3[2][0] = zsky3*v[2,0]
B3[2][1] = zsky3*v[2,1]
B3[2][2] = zsky3*v[2,2]

B = B1+B2+B3

U, S, VH = np.linalg.svd(B, full_matrices=True)
I = np.eye(3)
R1 = np.matmul(-U,I)
R2 = np.matmul(R1,VH)

d0 = mpmath.asin(-R2[2][2])*180/np.pi
a01 = mpmath.acos(-R2[0][2]/mpmath.cos(d0*np.pi/180))*180/np.pi
a02 = mpmath.asin(-R2[1][2]/mpmath.cos(d0*np.pi/180))*180/np.pi
if a02 < 0:
    a0 = -a01
else :
    a0 = a01

p00 = mpmath.atan2(R2[2][0],R2[2][1])*180/np.pi;

RA0 = a0
DE0 = d0
Roll0 = p00

final=[RA0,DE0,Roll0]
print(final)

#OUT.append({'RA':RA0,'dec':DE0,'roll':Roll0})
#np.savetxt(file_to,final,fmt='%.3f',delimiter=' ')
'''
