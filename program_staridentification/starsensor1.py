import numpy as np
import mpmath
from astropy.io import ascii
Catnew = ascii.read("./program_staridentification/Catalog_mod.txt")
Pixstars = ascii.read('./result_staridentification/StarPos_sorted/CG34.txt')
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

r012 = v2-v1
r013 = v3-v1
r023 = v3-v2
r014 = v4-v1
r024 = v4-v2

r12 = np.sqrt(np.vdot(r012,r012))
r13 = np.sqrt(np.vdot(r013,r013))
r23 = np.sqrt(np.vdot(r023,r023))
r14 = np.sqrt(np.vdot(r014,r014))
r24 = np.sqrt(np.vdot(r024,r024))

a1 = mpmath.acos((r12**2 + r13**2 - r23**2)/(2*r12*r13))*180/mpmath.pi
a2 = mpmath.acos((r12**2 + r23**2 - r13**2)/(2*r12*r23))*180/mpmath.pi
a3 = mpmath.acos((r13**2 + r23**2 - r12**2)/(2*r13*r23))*180/mpmath.pi
alfa0  = np.zeros(3)

alfa0[0] = a1
alfa0[1] = a2
alfa0[2] = a3
alfa = np.sort(alfa0,axis = 0)

a11 = mpmath.acos((r12**2 + r14**2 - r24**2)/(2*r12*r14))*180/mpmath.pi
a22 = mpmath.acos((r12**2 + r24**2 - r14**2)/(2*r12*r24))*180/mpmath.pi
a33 = mpmath.acos((r14**2 + r24**2 - r12**2)/(2*r14*r24))*180/mpmath.pi
alfa10 = np.zeros(3)
alfa10[0] = a11
alfa10[1] = a22
alfa10[2] = a33
alfa1 = np.sort(alfa10,axis = 0)

a = len(Catnew)
b = len(IDNx[0])
etresh = 0.1
A = np.zeros(a)
ux = 1
etresh = 0.0001
alfa0x = np.zeros(3)
alfa0t = np.zeros(3)
xx=0
for i in range(0, a - 1):
    if i == a:
        print('recognition failed')

    u1 = np.array([[mpmath.cos(RA[i]) * mpmath.cos(DE[i])], [mpmath.sin(RA[i]) * mpmath.cos(DE[i])], [mpmath.sin(DE[i])]])
    for j in range(0, b - 1):
        if IDNx[i][j] == 0:
            break
        jj = IDNx[i][j] - 1
        u2 = np.array([[mpmath.cos(RA[jj]) * mpmath.cos(DE[jj])], [mpmath.sin(RA[jj]) * mpmath.cos(DE[jj])],
                       [mpmath.sin(DE[jj])]])
        u = j + 1
        for k in range(u, b - 1):
            if IDNx[i][k] == 0:
                break
            kk = IDNx[i][k] - 1
            u3 = np.array([[mpmath.cos(RA[kk]) * mpmath.cos(DE[kk])], [mpmath.sin(RA[kk]) * mpmath.cos(DE[kk])],
                           [mpmath.sin(DE[kk])]])

            r12 = u2 - u1;
            r13 = u3 - u1;
            r23 = u3 - u2;

            absr12 = mpmath.sqrt(np.vdot(r12, r12))
            absr13 = mpmath.sqrt(np.vdot(r13, r13))
            absr23 = mpmath.sqrt(np.vdot(r23, r23))

            a1x = mpmath.acos((absr12 ** 2 + absr13 ** 2 - absr23 ** 2) / (2 * absr12 * absr13)) * 180 / mpmath.pi
            a2x = mpmath.acos((absr12 ** 2 + absr23 ** 2 - absr13 ** 2) / (2 * absr12 * absr23)) * 180 / mpmath.pi
            a3x = mpmath.acos((absr13 ** 2 + absr23 ** 2 - absr12 ** 2) / (2 * absr13 * absr23)) * 180 / mpmath.pi

            alfa0x[0] = a1x
            alfa0x[1] = a2x
            alfa0x[2] = a3x
            alfax = np.sort(alfa0x, axis=0)
            diff = alfax - alfa
            diffabs = np.absolute(diff)
            diffx = diffabs.sum(axis=0)
            error = diffx

            if error < etresh:
                error = etresh + 1
                xx = k
                while error > etresh:
                    xx = xx + 1
                    if xx > b - 1:
                        break
                    if IDNx[i][xx] == 0:
                        break
                    print(f"i: {i}, xx: {xx}")
                    xxk = IDNx[i][xx] - 1
                    u4 = np.array([[mpmath.cos(RA[xxk]) * mpmath.cos(DE[xxk])], [mpmath.sin(RA[xxk]) * mpmath.cos(DE[xxk])],
                                   [mpmath.sin(DE[xxk])]])

                    r14 = u4 - u1;
                    r24 = u4 - u2;

                    absr14 = mpmath.sqrt(np.vdot(r14, r14))
                    absr24 = mpmath.sqrt(np.vdot(r24, r24))

                    a1t = mpmath.acos((absr12 ** 2 + absr14 ** 2 - absr24 ** 2) / (2 * absr12 * absr14)) * 180 / mpmath.pi
                    a2t = mpmath.acos((absr12 ** 2 + absr24 ** 2 - absr14 ** 2) / (2 * absr12 * absr24)) * 180 / mpmath.pi
                    a3t = mpmath.acos((absr14 ** 2 + absr24 ** 2 - absr12 ** 2) / (2 * absr14 * absr24)) * 180 / mpmath.pi
                    alfa0t[0] = a1t
                    alfa0t[1] = a2t
                    alfa0t[2] = a3t
                    alfat = np.sort(alfa0t, axis=0)

                    error2 = alfa1 - alfat
                    errorx = np.absolute(error2)
                    errorxx = errorx.sum(axis=0)
                    error = errorxx

                    if error < etresh:
                        break
            if error < etresh:
                break
        if error < etresh:
            break
    if error < etresh:
        break


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

print(RA0)
print(DE0)
print(Roll0)
