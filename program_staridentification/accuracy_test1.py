import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Specify the file path
result = '/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/Result_sementara.csv'
initial= '/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/ID_rand.csv'
# Read the CSV file into a DataFrame
res = pd.read_csv(result)
init= pd.read_csv(initial)
N=1000
# Now you can work with the DataFrame (e.g., perform analysis, manipulate data)
errornum = np.zeros((N,4))
max_error = np.deg2rad(5)
count=0
blb=np.zeros(N)
tr = np.zeros(N)
fl = np.zeros(N)
def dist(ra1,ra2,de1,de2):
    x = np.arccos(np.sin(de1)*np.sin(de2)+np.cos(de1)*np.cos(de2)*np.cos(ra1-ra2))
    x = np.rad2deg(x)
    return x
for i in range(N):
    distance = dist(np.deg2rad(init.RA[i]),np.deg2rad(res.RA[i]),np.deg2rad(init.dec[i]),np.deg2rad(res.dec[i]))
    if distance < 5:
        correct = 1
        count+=1
        tr[i] = res.dist[i]
    else:
        correct = 0
        fl[i] = res.dist[i]
    blb[i]=i+1
    errornum[i]=[int(i+1),distance,int(correct),res.dist[i]]


print(errornum) 
print(count) # Print the first few rows of the DataFrame

plt.scatter(blb,fl)
plt.show
plt.grid(color='r', linestyle='-', linewidth=2)
plt.savefig('Tesfig.png')