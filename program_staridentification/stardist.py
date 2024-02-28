import numpy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits import mplot3d

A=pd.read_table('./program_staridentification/StarCatalogue.txt',sep=' ')
B=A.min(axis=0)

#print(B)

terangmax=A.loc[A['Magnitude'].idxmin()]
#print(terangmax)

df = pd.DataFrame(A)

#magnitude treshold to see stars
tres=4

StarsAboveTres = df[df['Magnitude']<tres]

#print(StarsAboveTres)

#plot distribusi bintang di langit
plt.scatter(numpy.rad2deg(StarsAboveTres['RA'])/15,numpy.rad2deg(StarsAboveTres['DE']),s=2**(-StarsAboveTres['Magnitude'])*10, alpha=0.5)

plt.grid()
plt.xlim(-12,12)
plt.ylim(-90,90)
plt.xlabel('RA (jam)')
plt.ylabel('dec (°)')
plt.title('Katalog SAO dengan magnitudo dibawah 4')
plt.gca().invert_xaxis()
plt.savefig('PetaLangit_mag4.png')

#peta bintang asli ini brok hehe
