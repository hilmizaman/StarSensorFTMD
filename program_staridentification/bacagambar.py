import cv2

# Specify the path to your image
image_path = "star_simulator-main/sample_images/Tes_Random/1_ra-31_de-65_roll110.jpg"

# Read the image
image = cv2.imread(image_path)

print(image.shape)
