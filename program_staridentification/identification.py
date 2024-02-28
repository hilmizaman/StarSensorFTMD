import numpy as np
import mpmath
from astropy.io import ascii


OUT=[]
for zzz in range(60):
    file_from='./result_staridentification/StarPos_sorted/CG'+str(zzz+1)+'.txt'
    file_to='./result_staridentification/Result/CG'+str(zzz+1)+'.txt'

    Catnew = ascii.read("./program_staridentification/Catalog_mod.txt")
    Pixstars = ascii.read(file_from)
    IDNx = ascii.read('./program_staridentification/IDNx_new.txt')
    xp = Pixstars['col1']
    yp = Pixstars['col2']
    RA = Catnew['col1']
    DE = Catnew['col2']
    MA = Catnew['col3']

    l = 3280
    w = 2464

    FOVx = 62.2
    FOVy = 48.8

    f = 0.00345597

    xtot = 2 * mpmath.tan((FOVx * mpmath.pi / 180) / 2) * f;
    ytot = 2 * mpmath.tan((FOVy * mpmath.pi / 180) / 2) * f;
    xpixel = l / xtot;
    ypixel = w / ytot;

    x1 = xp[0]
    x2 = xp[1]
    x3 = xp[2]
    x4 = xp[3]

    y1 = yp[0]
    y2 = yp[1]
    y3 = yp[2]
    y4 = yp[3]

    v1 = np.zeros(3)
    xplane1 = (x1-l/2)*xtot/l
    yplane1 = (-y1+w/2)*ytot/w
    v1[0] = xplane1/np.sqrt(xplane1**2 + yplane1**2 +f**2)
    v1[1] = yplane1/np.sqrt(xplane1**2 + yplane1**2 +f**2)
    v1[2] = f/np.sqrt(xplane1**2 + yplane1**2 +f**2)

    v2 = np.zeros(3)
    xplane2 = (x2-l/2)*xtot/l
    yplane2 = (-y2+w/2)*ytot/w
    v2[0] = xplane2/np.sqrt(xplane2**2 + yplane2**2 +f**2)
    v2[1] = yplane2/np.sqrt(xplane2**2 + yplane2**2 +f**2)
    v2[2] = f/np.sqrt(xplane2**2 + yplane2**2 +f**2)

    v3 = np.zeros(3)
    xplane3 = (x3-l/2)*xtot/l
    yplane3 = (-y3+w/2)*ytot/w
    v3[0] = xplane3/np.sqrt(xplane3**2 + yplane3**2 +f**2)
    v3[1] = yplane3/np.sqrt(xplane3**2 + yplane3**2 +f**2)
    v3[2] = f/np.sqrt(xplane3**2 + yplane3**2 +f**2)

    v4 = np.zeros(3)
    xplane4 = (x4-l/2)*xtot/l
    yplane4 = (-y4+w/2)*ytot/w
    v4[0] = xplane4/np.sqrt(xplane4**2 + yplane4**2 +f**2)
    v4[1] = yplane4/np.sqrt(xplane4**2 + yplane4**2 +f**2)
    v4[2] = f/np.sqrt(xplane4**2 + yplane4**2 +f**2)

    alfa12 = mpmath.acos(np.vdot(v1,v2))
    alfa13 = mpmath.acos(np.vdot(v1,v3))
    alfa14 = mpmath.acos(np.vdot(v1,v4))
    alfa23 = mpmath.acos(np.vdot(v2,v3))
    alfa24 = mpmath.acos(np.vdot(v2,v4))
    alfa34 = mpmath.acos(np.vdot(v3,v4))




    u1 = np.zeros(3)
    u2 = np.zeros(3)
    u3 = np.zeros(3)
    u4 = np.zeros(3)

    ux = 1
    etresh = 0.0001
    alfa0x = np.zeros(3)
    alfa0t = np.zeros(3)
    a = len(Catnew)
    b = len(IDNx[0])
    error14 = 10.0001;
    error24 = 10.0001;
    error34 = 10.0001;
    for i in range(0,a-1):
        if i == a:
            print('recognition failed')

        u1[0] = mpmath.cos(RA[i])*mpmath.cos(DE[i])
        u1[1] = mpmath.sin(RA[i])*mpmath.cos(DE[i])
        u1[2] = mpmath.sin(DE[i])

        xx=0
        for j in range(0,b-1):
            if IDNx[i][j]-1 == 0:
                break
            jj = IDNx[i][j]-1
            u2[0] = mpmath.cos(RA[jj]) * mpmath.cos(DE[jj])
            u2[1] = mpmath.sin(RA[jj]) * mpmath.cos(DE[jj])
            u2[2] = mpmath.sin(DE[jj])

            ax12 = mpmath.acos(np.vdot(u1,u2))
            errorx12 = alfa12 - ax12
            error12 = np.absolute(errorx12)
            if error12 < etresh:
                for k in range(j+1,b-1):
                    kk = IDNx[i][k] - 1
                    if IDNx[i][k]-1 == 0:
                        break
                    u3[0] = mpmath.cos(RA[kk]) * mpmath.cos(DE[kk])
                    u3[1] = mpmath.sin(RA[kk]) * mpmath.cos(DE[kk])
                    u3[2] = mpmath.sin(DE[kk])

                    ax13 = mpmath.acos(np.vdot(u1,u3))
                    ax23 = mpmath.acos(np.vdot(u2,u3))
                    errorx13 = alfa13 - ax13
                    errorx23 = alfa23 - ax23
                    error13 = np.absolute(errorx13)
                    error23 = np.absolute(errorx23)
                    
                    if error13 < etresh and error23 < etresh:
                        for xx in range(k+1,b-1):
                            xxk = IDNx[i][xx] - 1
                            if IDNx[i][xx]-1 == 0:
                                break
                            u4[0] = mpmath.cos(RA[xxk]) * mpmath.cos(DE[xxk])
                            u4[1] = mpmath.sin(RA[xxk]) * mpmath.cos(DE[xxk])
                            u4[2] = mpmath.sin(DE[xxk])

                            ax14 = mpmath.acos(np.vdot(u1,u4))
                            ax24 = mpmath.acos(np.vdot(u2,u4))
                            ax34 = mpmath.acos(np.vdot(u3,u4))
                            errorx14 = alfa14 - ax14
                            errorx24 = alfa24 - ax24
                            errorx34 = alfa34 - ax34
                            error14 = np.absolute(errorx14)
                            error24 = np.absolute(errorx24)
                            error34 = np.absolute(errorx34)

                            if error14 < etresh and error24 < etresh and error34 < etresh:
                                break
                    if error14 < etresh and error24 < etresh and error34 < etresh:
                        break
            if error14 < etresh and error24 < etresh and error34 < etresh:
                break
        if error14 < etresh and error24 < etresh and error34 < etresh:
            break

    if i == 133:
        print('failed')



    id1 = i
    id2 = IDNx[i][j]-1
    id3 = IDNx[i][k]-1
    id4 = IDNx[i][xx]-1

    xsky1 = mpmath.cos(RA[i])*mpmath.cos(DE[i]);
    ysky1 = mpmath.sin(RA[i])*mpmath.cos(DE[i]);
    zsky1 = mpmath.sin(DE[i]);

    xsky2 = mpmath.cos(RA[id2])*mpmath.cos(DE[id2]);
    ysky2 = mpmath.sin(RA[id2])*mpmath.cos(DE[id2]);
    zsky2 = mpmath.sin(DE[id2]);

    xsky3 = mpmath.cos(RA[id3])*mpmath.cos(DE[id3]);
    ysky3 = mpmath.sin(RA[id3])*mpmath.cos(DE[id3]);
    zsky3 = mpmath.sin(DE[id3]);

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

    B1 = np.zeros((3,3))
    B2 = np.zeros((3,3))
    B3 = np.zeros((3,3))

    B1[0][0] = xsky1*v1[0]
    B1[0][1] = xsky1*v1[1]
    B1[0][2] = xsky1*v1[2]
    B1[1][0] = ysky1*v1[0]
    B1[1][1] = ysky1*v1[1]
    B1[1][2] = ysky1*v1[2]
    B1[2][0] = zsky1*v1[0]
    B1[2][1] = zsky1*v1[1]
    B1[2][2] = zsky1*v1[2]

    B2[0][0] = xsky2*v2[0]
    B2[0][1] = xsky2*v2[1]
    B2[0][2] = xsky2*v2[2]
    B2[1][0] = ysky2*v2[0]
    B2[1][1] = ysky2*v2[1]
    B2[1][2] = ysky2*v2[2]
    B2[2][0] = zsky2*v2[0]
    B2[2][1] = zsky2*v2[1]
    B2[2][2] = zsky2*v2[2]

    B3[0][0] = xsky3*v3[0]
    B3[0][1] = xsky3*v3[1]
    B3[0][2] = xsky3*v3[2]
    B3[1][0] = ysky3*v3[0]
    B3[1][1] = ysky3*v3[1]
    B3[1][2] = ysky3*v3[2]
    B3[2][0] = zsky3*v3[0]
    B3[2][1] = zsky3*v3[1]
    B3[2][2] = zsky3*v3[2]

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
    OUT.append({'RA':RA0,'dec':DE0,'roll':Roll0})
    np.savetxt(file_to,final,fmt='%.3f',delimiter=' ')

import csv

# Specify the file name
file_name = 'Result.csv'

# Open the file in write mode
with open(file_name, 'w', newline='') as csvfile:
    # Define the header (column names)
    fieldnames = ['RA', 'dec', 'roll']
    
    # Create a CSV writer
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Write the rows
    writer.writerows(OUT)

print(f'Table has been written to {file_name}')
