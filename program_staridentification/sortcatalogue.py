import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt


neighbor_dist = ascii.read("./program_staridentification/IDNx_6mag_dist.txt")
login=np.array(neighbor_dist)
Distance2 = np.array(sorted(login, key=lambda x: x[1]))

#np.savetxt('./program_staridentification/IDNx_6mag_dist_sorted.txt',Distance2,fmt='%.9f',delimiter=' ')
#print(Distance2)
xpoints=neighbor_dist['col1']
ypoints=neighbor_dist['col2']
ypoints=sorted(ypoints)

print(xpoints)
print(ypoints)
plt.plot(xpoints, ypoints)
plt.savefig('dist.png')