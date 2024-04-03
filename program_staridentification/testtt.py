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


