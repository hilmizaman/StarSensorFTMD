import csv

path='./program_starcentroiding/StarCatalogue/Below_6.0_SAO.csv'
path_to='./program_staridentification/Catalog_mod2.txt'
file_from=path
file_to=path_to

with open(file_from, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader, None)
        with open(file_to, 'w') as txt_file:
            for row in csv_reader:
                del row[0:2]

                # Assuming each element in a row is separated by a space
                txt_file.write(' '.join(row) + '\n')
