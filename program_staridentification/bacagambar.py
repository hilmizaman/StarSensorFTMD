import cv2
import numpy as np

# Specify the path to your image
image_path = "star_simulator-main/sample_images/Tes_Random_new/1_ra-85_de76_roll158.jpg"

# Read the image
image = cv2.imread(image_path)

# Check if the image is loaded successfully
if image is not None:
    # Convert the image to grayscale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Create a matrix to represent the image
    image_matrix = np.array(gray_image)

    # Save the image matrix to a text file
    np.savetxt('image_matrix.txt', image_matrix, fmt='%d')

    print("Image Matrix Shape:", image_matrix.shape)
    print("Image Matrix:")
    print(image_matrix)

else:
    print('Error: Unable to load the image.')
