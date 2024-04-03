import numpy as np
import mpmath
from astropy.io import ascii
import csv

def csv_to_txt(csv_filename, txt_filename):
    with open(csv_filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        with open(txt_filename, 'w') as txt_file:
            for row in csv_reader:
                # Assuming each element in a row is separated by a space
                txt_file.write(' '.join(row) + '\n')

# Example usage
for i in range(1000):
    path='/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/result_starcentroiding/starPositionCalculated/firstMethod/CG'
    #path='/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/result_starcentroiding/starPositionCalculated/secondMethod/WCG'
    file='.csv'
    num=str(i+1)
    path_to='./result_staridentification/StarPos_sorted/CG'
    #path_to='./result_staridentification/StarPos_sorted/WCG'
    file_to='.txt'
    file_from=path+num+file
    file_to=path_to+num+file_to

    csv_to_txt(file_from, './tulis.txt')

    # Read coordinates from the file
    with open('./tulis.txt', "r") as file:
        lines = file.readlines()

    # Parse coordinates from lines
    coordinates = np.array([list(map(float, line.strip().split())) for line in lines])

    # Calculate distances from (0,0)
    distances = np.linalg.norm(coordinates, axis=1)

    # Find the index of the point closest to (0,0)
    closest_index = np.argmin(distances)
    closest_point = coordinates[closest_index]

    # Remove the closest point from the original data
    coordinates_without_closest = np.delete(coordinates, closest_index, axis=0)

    # Calculate distances from the closest point
    distances_from_closest = np.linalg.norm(coordinates_without_closest - closest_point, axis=1)

    # Combine coordinates and distances for sorting
    data_with_distances = np.column_stack((coordinates_without_closest, distances_from_closest))

    # Sort based on distances
    sorted_data = data_with_distances[data_with_distances[:, 2].argsort()]

    # Insert the closest point at the top
    sorted_data = np.insert(sorted_data, 0, np.append(closest_point, 0), axis=0)

    np.savetxt(file_to,sorted_data[:,0:2],fmt='%.3f',delimiter=' ')

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

OUT = []
for zz in range(1000):
    #input data bintang
    file_from='./result_staridentification/StarPos_sorted/CG'+str(zz+1)+'.txt'
    file_to='./result_staridentification/Result/CG1_binary'+'.txt'

    Catnew = ascii.read("./program_staridentification/Catalog_mod2.txt")
    Pixstars = ascii.read(file_from)
    neighbor = ascii.read('./program_staridentification/Distance_sorted.txt')
    dist12 = np.array(neighbor['col2'])
    IDNx = ascii.read('./program_staridentification/IDNx2.txt')

    #Data koordinat bintang di CCD
    xp = Pixstars['col1']
    yp = Pixstars['col2']
    #Data koordinat bintang di langit
    RA = Catnew['col1']
    DE = Catnew['col2']
    MA = Catnew['col3']

    #input berapa bintang yang digunakan
    N_stars = 7

    #fungsi untuk jarak antar bintang di CCD
    def CCD_dist(x,y):
        xplane = (x)*xtot/l
        yplane = (-y)*ytot/w
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

    print(zz+1)

    #error maks untuk jarak paling dekat
    error1 = 0.002

    #menentukan batas atas dan bawah
    upper_lim = alfa1[0]+error1
    lower_lim = alfa1[0]-error1

    #print(alfa1)

    #error maks untuk jarak paling dekat
    id_sub = -9363479276.48673*alfa1[0]**5+2529712831.02641*alfa1[0]**4-239653974.935343*alfa1[0]**3+8305048.16576691*alfa1[0]**2+21270.2957076176*alfa1[0]+86.1151602836825
    #print(id_sub)
    idup = id_sub+150
    idlow= id_sub-150
    if id_sub < 125:
        idup = id_sub+250
    '''error1 = 0.0003
    error2 = error1/10
    if alfa1[0]<0.02:
        error1 = 0.0005
        error2 = error1/10
    if alfa1[0]>0.035:
        error1 = 0.0004
        error2 = error1/10
    if alfa1[0]>0.04:
        error1 = 0.00055
        error2 = error1/8
    if alfa1[0]<0.015:
        error1 = 0.0005
        error2 = error1/5
    if alfa1[0]<0.005:
        error1 = 0.001
        error2 = error1/2    
    if alfa1[0]<0.001:
        error1 = 0.001
        error2 = error1    
    if alfa1[0]>0.03:
        error1 = 0.0006
        error2 = error1/10
    if alfa1[0]>0.05:
        error1 = 0.001
        error2 = error1/5    
    if alfa1[0]>0.06:
        error1 = 0.001
        error2 = error1/2'''

    #menentukan batas atas dan bawah
    upper_lim = 2.16103801737176e-19*idup**5-2.46155782506950e-15*idup**4+1.04160334113478e-11*idup**3-1.99440998195988e-08*idup**2+2.41921700426570e-05*idup-0.00168330254604258
    lower_lim = 2.16103801737176e-19*idlow**5-2.46155782506950e-15*idlow**4+1.04160334113478e-11*idlow**3-1.99440998195988e-08*idlow**2+2.41921700426570e-05*idlow-0.00168330254604258
    if abs(upper_lim-alfa1[0])>abs(lower_lim-alfa1[0]):
        error1 = abs(upper_lim-alfa1[0])
        error2 = error1/10
    else:
        error1 = abs(lower_lim-alfa1[0])
        error2 = error1/10

    '''upper_lim = alfa1[0]+error1
    lower_lim = alfa1[0]-error1'''

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
    N_upper = binary_search(dist12,upper_lim,error2)
    N_lower = binary_search(dist12,lower_lim,error2)

    #print([N_lower,N_upper])
    if N_lower ==-1:
        N_lower = int(idlow)
        if int(idlow)<0:
            N_lower = 0
    if N_upper ==-1:
        N_upper = int(idup)
        if int(idup):
            N_upper = 5102
    #print([N_lower,N_upper])
    if N_upper > 5000:
        N_upper = 5102
    if N_lower < 100:
        N_lower = 0
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

    IDN = np.array(neighbor['col1'][low:up])
    dis2 = np.array(neighbor['col2'][low:up])
    dis3 = np.array(neighbor['col3'][low:up])
    dis4 = np.array(neighbor['col4'][low:up])
    dis5 = np.array(neighbor['col5'][low:up])
    dis6 = np.array(neighbor['col6'][low:up])
    dis7 = np.array(neighbor['col7'][low:up])

    Dist_345 = np.array([IDN,dis2,dis3,dis4,dis5,dis6,dis7])
    Disttrak = np.array([IDN,dis2,dis3,dis4,dis5,dis6,dis7])
    #print(Dist_345)
    dis=[]
    Starnum=[]
    for ii in range(N_stars-2):
        Star = sort_2d_array(Dist_345,n=ii+2)
        #print(Star[0][20:50])
        upper_lim1 = alfa1[ii+1]+error1*(ii+3)*2
        lower_lim1 = alfa1[ii+1]-error1*(ii+3)*2
        #mencari batas atas dan bawah 
        uppa = binary_search(Star[ii+2],upper_lim1,error1*(ii+3)*1)
        lowwa = binary_search(Star[ii+2],lower_lim1,error1*(ii+3)*1)
        if uppa==-1 or lowwa==-1:
            #print('iter no '+str(zz+1)+' gagal bro pas ii='+str(ii))
            break 
        #print([low,up])
        up = uppa
        low = lowwa
        Dist_345 = np.array([Star[0][low:up],Star[1][low:up],Star[2][low:up],Star[3][low:up],Star[4][low:up],Star[5][low:up],Star[6][low:up]])
        #print(Dist_345)

    error_tot = 99999
    weight = [1,1,1,1,1,1]
    tot = sum(np.multiply(alfa1,weight))
    #print([up,low])
    ID_fin = []
    for num in range(up-low):
        error_totz = 0
        mat = [Dist_345[1][num],Dist_345[2][num],Dist_345[3][num],Dist_345[4][num],Dist_345[5][num],Dist_345[6][num]]
        for kz in range(6):
            error_totz = error_totz+abs(alfa1[kz]-mat[kz])*weight[kz]
        #print([Dist_345[0][num],error_totz])
        if error_totz < error_tot:
            error_tot = error_totz
            ID_fin = Dist_345[0][num]


    #print(np.shape(alfa1))
    alfa11 = alfa1.reshape(1,6)
    #print(np.shape(alfa11))
    if error_tot>0.005:
        Disnext = Disttrak.reshape((N_upper-N_lower),N_stars)
        #print(np.shape(Disttrak))
        #print(np.shape(Disnext))
        for kz in range(N_stars-3):

            #print(kz)
            Disnextt = np.delete(Disttrak,kz+1,axis=0)
            alfax = alfa11
            Dist_3456 = np.array(Disnextt)
            #print(alfax)
            for ii in range(N_stars-3):
                Star = sort_2d_array(Dist_3456,n=ii+2)
                #print(Star[0][20:50])
                upper_lim1 = alfax[0][ii+1]+error1*(ii+3)*2
                lower_lim1 = alfax[0][ii+1]-error1*(ii+3)*2
                #mencari batas atas dan bawah 
                uppa = binary_search(Star[ii+2],upper_lim1,error1*(ii+3))
                lowwa = binary_search(Star[ii+2],lower_lim1,error1*(ii+3))
                if uppa==-1 or lowwa==-1:
                    #print('iter no '+str(zz+1)+' gagal bro pas ii='+str(ii))
                    break 
                #print([low,up])
                up = uppa
                low = lowwa
                Dist_3456 = np.array([Star[0][low:up],Star[1][low:up],Star[2][low:up],Star[3][low:up],Star[4][low:up],Star[5][low:up]])
            #print([low,up])    
            #print(Dist_3456) 
            for num in range(up-low):
                error_totz = 0
                mat = [Dist_3456[1][num],Dist_3456[2][num],Dist_3456[3][num],Dist_3456[4][num],Dist_3456[5][num]]
                for kk in range(5):
                    error_totz = error_totz+abs(alfax[0][kk]-mat[kk])*weight[kk]
                #print([Dist_3456[0][num],error_totz,kz+2])
                if error_totz < error_tot-0.0003:
                    error_tot = error_totz
                    ID_fin = Dist_3456[0][num]

            #print([ID_fin,error_tot])
            #print('BATES-----------------------')


    #print(Dist_345)
    #print(int(ID_fin-1))
    #print(IDNx['col1'][int(ID_fin-1)])
    ID_fin = int(ID_fin)
    ID_fin2 = int(IDNx['col1'][int(ID_fin-1)])
    ID_fin3 = int(IDNx['col2'][int(ID_fin-1)])
    #print([ID_fin,ID_fin2,ID_fin3])


    def skycord(ID):
        x = np.cos(RA[ID])*np.cos(DE[ID])
        y = np.sin(RA[ID])*np.cos(DE[ID])
        z = np.sin(DE[ID])
        sky = [x,y,z]
        return sky

    N_wahba = 3
    vv = np.zeros((N_wahba,3))
    vv[0][:] = skycord(ID_fin-1)
    for i in range(N_wahba-1):
        vv[i+1][:]=skycord(int(IDNx['col'+str(i+1)][int(ID_fin-1)])-1)

    #print(vv)
    #mencari koordinat di langit


    #wahba problem, maksimalkan nilai B = a*b*r'
    #print(np.outer(vv[0],v[0]))

    #BB=np.zeros((3,3,N_stars-1))
    #print(np.shape(BB))
    B=np.zeros((3,3))
    a = 1/(N_stars-1)
    #print(a)
    for hi in range(N_wahba):
        B+=a*np.outer(vv[hi],v[hi])

    #print(B)

    #print(BB)
    #print(B)

    #pakai SVD buat nyari matriks A dr B
    U, S, VT = np.linalg.svd(B)
    uplus = U @ np.diag([1,1,np.linalg.det(U)])
    vplus = VT.T @ np.diag([1,1,np.linalg.det(VT.T)])
    #matriks A di wahba problem (R karena jadi matriks rotasi jg) 
    R = vplus @ uplus.T

    #for some reason, hasil dari matriks rotasi yg dihasilkan oleh algoritma SVD ini jadi 
    # [[-21 -22 -23],[31 32 33],[-11 -12 -13]] dari matriks rotasi kalo kita putar arah ZYX sebesar [RA, -de, roll]



    d0 = -np.arcsin(R[2][2])*180/np.pi
    a0 = np.arctan2(-R[2][1],-R[2][0])*180/np.pi
    p00 = np.arctan2(-R[0][2],R[1][2])*180/np.pi

    RA0 = a0
    DE0 = d0
    Roll0 = p00
    #print([np.rad2deg(RA[ID_fin-1]),np.rad2deg(DE[ID_fin-1])]) #sudah bener
    OUT.append({'ID':ID_fin,'RA':RA0,'dec':DE0,'roll':Roll0,'dist':alfa1[0]})


import csv
# Specify the file name
file_name = 'Result_sementara2.csv'

# Open the file in write mode
with open(file_name, 'w', newline='') as csvfile:
    # Define the header (column names)
    fieldnames = ['ID','RA', 'dec','roll','dist']
    
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
