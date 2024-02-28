import numpy as np

# Load the text file into an array
data_array = np.loadtxt('image_matrix.txt')
matout=data_array*0
#print(matout)
tot=0

# Print the array
for i in range(np.size(data_array,0)):
    for j in range(np.size(data_array,1)):
        if data_array[i,j]>100:
            matout[i,j]=1
            tot+=1
print(i)
print(j)
print(tot)
            
