import numpy as np

# Sample matrix representing coordinates
coordinates_matrix = np.array([[1, 2], [1, 4], [1, 5], [1, -3]])

# Calculate the distance of each coordinate from the origin
distances = np.sqrt(np.sum(coordinates_matrix**2, axis=1))

# Pair each coordinate with its corresponding distance
coordinate_distance_pairs = list(zip(coordinates_matrix, distances))

# Sort the pairs based on distances and then coordinates
sorted_pairs = sorted(coordinate_distance_pairs, key=lambda x: (x[1], x[0][0], x[0][1]))

# Extract the sorted coordinates
sorted_coordinates = [pair[0] for pair in sorted_pairs]

print("Sorted coordinates based on distance from origin and then coordinates:")
print(sorted_coordinates)
