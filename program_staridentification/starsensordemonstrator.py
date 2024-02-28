import numpy as np
import mpmath
from astropy.io import ascii
Catnew = ascii.read("Catdem.txt")
Pixstars = ascii.read('Pixstarsdem.txt')
IDNx = ascii.read('IDNxdem.txt')
xp = Pixstars['col1']
yp = Pixstars['col2']
RA = Catnew['col1']
DE = Catnew['col2']
MA = Catnew['col3']

l = 3280
w = 2464

FOVx = 62.2
FOVy = 48.8

f = 0.003044898

a = len(xp)
b = len(RA)

R0 = np.zeros(3)
R1 = np.zeros(3)
dr = np.zeros(a)
r = np.zeros((3,a))
r0 = np.zeros(3)
r1 = np.zeros(3)
drr = np.zeros(a)
u1 = np.zeros(3)
u2 = np.zeros(3)
u3 = np.zeros(3)
u4 = np.zeros(3)
alfa0 = np.zeros(3)
alfa = np.zeros(3)
alfat0 = np.zeros(3)
alfat = np.zeros(3)
ImageID0 = np.zeros(a)
ImageID = np.zeros(a)
v1 = np.zeros(3)
v2 = np.zeros(3)
v3 = np.zeros(3)
v4 = np.zeros(3)
alfax0 = np.zeros(3)
alfax = np.zeros(3)
alfatx0 = np.zeros(3)
alfatx = np.zeros(3)
error = np.zeros(b)

xtot = 2*mpmath.tan((FOVx*mpmath.pi/180)/2)*f;
ytot = 2*mpmath.tan((FOVy*mpmath.pi/180)/2)*f;
for i in range (0,a):
    xplane = (xp[i]-l/2)*xtot/l
    yplane = (-yp[i]+w/2)*ytot/w
    absR = mpmath.sqrt(xplane ** 2 + yplane ** 2 + f ** 2)
    r[0][i] = xplane/absR
    r[1][i] = yplane/absR
    r[2][i] = f/absR
    R0[0] = xplane/absR
    R0[1] = yplane/absR
    R0[2] = f/absR
    R1[0] = 0
    R1[1] = 0
    R1[2] = 1
    dr[i] = mpmath.acos(np.vdot(R0,R1))

g = np.argmin(dr)
print(dr)
for i in range(0,a):
    if i != g:
        r0[0] = r[0][g]
        r0[1] = r[1][g]
        r0[2] = r[2][g]
        r1[0] = r[0][i]
        r1[1] = r[1][i]
        r1[2] = r[2][i]
        drr[i] = mpmath.acos(np.vdot(r0,r1))
    else:
        drr[i] = 0




for i in range(0,a):
    h = np.argmin(drr)
    drr[h] = 100
    ImageID0[i] = h

ImageID = ImageID0.astype(int)
x1 = ImageID[0]
x2 = ImageID[1]
x3 = ImageID[2]
x4 = ImageID[3]

u1[0] = r[0][x1]
u1[1] = r[1][x1]
u1[2] = r[2][x1]

u2[0] = r[0][x2]
u2[1] = r[1][x2]
u2[2] = r[2][x2]

u3[0] = r[0][x3]
u3[1] = r[1][x3]
u3[2] = r[2][x3]

u4[0] = r[0][x4]
u4[1] = r[1][x4]
u4[2] = r[2][x4]

r12 = u2-u1;
r13 = u3-u1;
r23 = u3-u2;

r14 = u4-u1;
r24 = u4-u2;

absr12 = mpmath.sqrt(np.vdot(r12,r12));
absr13 = mpmath.sqrt(np.vdot(r13,r13));
absr23 = mpmath.sqrt(np.vdot(r23,r23));

alfa1 = mpmath.acos((absr12**2 + absr13**2 - absr23**2)/(2*absr12*absr13))*180/mpmath.pi;
alfa2 = mpmath.acos((absr12**2 + absr23**2 - absr13**2)/(2*absr12*absr23))*180/mpmath.pi;
alfa3 = mpmath.acos((absr13**2 + absr23**2 - absr12**2)/(2*absr13*absr23))*180/mpmath.pi;

alfa0[0] = alfa1
alfa0[1] = alfa2
alfa0[2] = alfa3

alfa = np.sort(alfa0,axis=0)

absr14 = mpmath.sqrt(np.vdot(r14,r14));
absr24 = mpmath.sqrt(np.vdot(r24,r24));

alfat1 = mpmath.acos((absr12**2 + absr14**2 - absr24**2)/(2*absr12*absr14))*180/mpmath.pi;
alfat2 = mpmath.acos((absr12**2 + absr24**2 - absr14**2)/(2*absr12*absr24))*180/mpmath.pi;
alfat3 = mpmath.acos((absr14**2 + absr24**2 - absr12**2)/(2*absr14*absr24))*180/mpmath.pi;

alfat0[0] = alfat1
alfat0[1] = alfat2
alfat0[2] = alfat3

alfat = np.sort(alfat0,axis=0)

for i in range(0,b):
    x1 = IDNx[i][0]-1
    x2 = IDNx[i][1]-1
    x3 = IDNx[i][2]-1
    x4 = IDNx[i][3]-1

    v1[0] = mpmath.cos(RA[x1]) * mpmath.cos(DE[x1])
    v1[1] = mpmath.sin(RA[x1]) * mpmath.cos(DE[x1])
    v1[2] = mpmath.sin(DE[x1])
    v2[0] = mpmath.cos(RA[x2]) * mpmath.cos(DE[x2])
    v2[1] = mpmath.sin(RA[x2]) * mpmath.cos(DE[x2])
    v2[2] = mpmath.sin(DE[x2])
    v3[0] = mpmath.cos(RA[x3]) * mpmath.cos(DE[x3])
    v3[1] = mpmath.sin(RA[x3]) * mpmath.cos(DE[x3])
    v3[2] = mpmath.sin(DE[x3])
    v4[0] = mpmath.cos(RA[x4]) * mpmath.cos(DE[x4])
    v4[1] = mpmath.sin(RA[x4]) * mpmath.cos(DE[x4])
    v4[2] = mpmath.sin(DE[x4])

    r12 = v2 - v1;
    r13 = v3 - v1;
    r23 = v3 - v2;

    r14 = v4 - v1;
    r24 = v4 - v2;

    absr12 = mpmath.sqrt(np.vdot(r12, r12));
    absr13 = mpmath.sqrt(np.vdot(r13, r13));
    absr23 = mpmath.sqrt(np.vdot(r23, r23));

    alfa1 = mpmath.acos((absr12**2 + absr13**2 - absr23**2) / (2 * absr12 * absr13)) * 180 / mpmath.pi;
    alfa2 = mpmath.acos((absr12**2 + absr23**2 - absr13**2) / (2 * absr12 * absr23)) * 180 / mpmath.pi;
    alfa3 = mpmath.acos((absr13**2 + absr23**2 - absr12**2) / (2 * absr13 * absr23)) * 180 / mpmath.pi;

    alfax0[0] = alfa1
    alfax0[1] = alfa2
    alfax0[2] = alfa3

    alfax = np.sort(alfax0,axis=0)

    absr14 = mpmath.sqrt(np.vdot(r14, r14));
    absr24 = mpmath.sqrt(np.vdot(r24, r24));

    alfat1 = mpmath.acos((absr12**2 + absr14**2 - absr24**2) / (2 * absr12 * absr14)) * 180 / mpmath.pi;
    alfat2 = mpmath.acos((absr12**2 + absr24**2 - absr14**2) / (2 * absr12 * absr24)) * 180 / mpmath.pi;
    alfat3 = mpmath.acos((absr14**2 + absr24**2 - absr12**2) / (2 * absr14 * absr24)) * 180 / mpmath.pi;

    alfatx0[0] = alfat1
    alfatx0[1] = alfat2
    alfatx0[2] = alfat3
    alfatx = np.sort(alfatx0,axis=0)

    error1 = alfa - alfax
    error1x = np.absolute(error1)
    error1xx = error1x.sum(axis=0)

    error2 = alfat - alfatx
    error2x = np.absolute(error2)
    error2xx = error2x.sum(axis=0)

    error[i] = error1xx + error2xx

v = np.argmin(error)

id1 = v
id2 = IDNx[v][1]-1
id3 = IDNx[v][2]-1

xsky1 = mpmath.cos(RA[id1])*mpmath.cos(DE[id1]);
ysky1 = mpmath.sin(RA[id1])*mpmath.cos(DE[id1]);
zsky1 = mpmath.sin(DE[id1]);

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

B1[0][0] = xsky1*u1[0]
B1[0][1] = xsky1*u1[1]
B1[0][2] = xsky1*u1[2]
B1[1][0] = ysky1*u1[0]
B1[1][1] = ysky1*u1[1]
B1[1][2] = ysky1*u1[2]
B1[2][0] = zsky1*u1[0]
B1[2][1] = zsky1*u1[1]
B1[2][2] = zsky1*u1[2]

B2[0][0] = xsky2*u2[0]
B2[0][1] = xsky2*u2[1]
B2[0][2] = xsky2*u2[2]
B2[1][0] = ysky2*u2[0]
B2[1][1] = ysky2*u2[1]
B2[1][2] = ysky2*u2[2]
B2[2][0] = zsky2*u2[0]
B2[2][1] = zsky2*u2[1]
B2[2][2] = zsky2*u2[2]

B3[0][0] = xsky3*u3[0]
B3[0][1] = xsky3*u3[1]
B3[0][2] = xsky3*u3[2]
B3[1][0] = ysky3*u3[0]
B3[1][1] = ysky3*u3[1]
B3[1][2] = ysky3*u3[2]
B3[2][0] = zsky3*u3[0]
B3[2][1] = zsky3*u3[1]
B3[2][2] = zsky3*u3[2]

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
