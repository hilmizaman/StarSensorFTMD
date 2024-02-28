import numpy
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from astropy.io import ascii

Pixstars = ascii.read('./program_staridentification/WCG34.txt')
xp = Pixstars['col1']
yp = Pixstars['col2']


#plot distribusi bintang di langit
plt.scatter(xp,yp,)
plt.grid()
plt.xlim(-12,12)
plt.ylim(-90,90)

plt.savefig('CG34.png')

#peta bintang asli ini brok hehe
