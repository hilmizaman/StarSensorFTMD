import csv

def csv_to_txt(csv_filename, txt_filename):
    with open(csv_filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        with open(txt_filename, 'w') as txt_file:
            for row in csv_reader:
                # Assuming each element in a row is separated by a space
                txt_file.write(' '.join(row) + '\n')

# Example usage
for i in range(100):
    path='/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/result_starcentroiding/starPositionCalculated/firstMethod/CG'
    file='.csv'
    num=str(i+1)
    path_to='./result_staridentification/StarPosCalc/CG'
    file_to='.txt'
    file_from=path+num+file
    file_to=path_to+num+file_to
    csv_to_txt(file_from, file_to)
