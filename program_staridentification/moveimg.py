import shutil
import pandas as pd

# Specify the path to your CSV file
csv_file_path = 'RA_dec_roll2.csv'

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Display the DataFrame
print(df)


for i in range(100):

    # Specify the source file path and name
    source_file = 'star_simulator-main/sample_images/Tes_Random_new/'+str(i+1)+'_ra'+str(df.iloc[i,0])+'_de'+str(df.iloc[i,1])+'_roll'+str(df.iloc[i,2])+'.jpg'

    # Specify the destination folder and new file name
    destination_folder = 'result_starcentroiding/starImages2/'
    new_file_name = 'stars'+str(i+1)+'.png'

    # Construct the full destination path
    destination_path = destination_folder + new_file_name

    # Move and rename the file
    shutil.copy(source_file, destination_path)
