import glob
import os
from os.path import exists
import csv
from csv import DictReader
import re
import pandas as pd
from operator import itemgetter


path1 = os.getcwd()
print(path1)

#Defining csv's being loaded in
filename = 'ThermoBase_ver_1.0_2022.csv'
accn_file = 'Accessions.csv'

#create final output file
Thermobase_matches = open('Thermobase_matches.csv', 'w')
Thermobase_matches.write(f"accn,species,temp\n")

#opening thermobase csv and making it into list of dictionary
with open(filename, 'r') as f:
    f_dict_reader = DictReader(f)
    f_list_of_dict = list(f_dict_reader)
    f_list_of_dict = f_list_of_dict[:]
    #removing commas and replacing white-space with underscores in species name
    for i in f_list_of_dict:
        name = i.get('Name')
        i.update({'Name': name.replace(' ', '_').replace(',', '')})
        #print(i)
    #print(f_list_of_dict)

#opening accessions csv
with open(accn_file,'r') as data:
    test = csv.reader(data)
    cnt = 0
    #creating a list for each line of csv
    for line in test:
        cnt += 1
        #print(line)
        accn = line[0]
        #print(accn)
        species = line[2]
        #print(species)
        try:
            print(f'doing {cnt} of 730')
            match = next(h for h in f_list_of_dict if h['Name'] == species)
            print(match)
            temp = match.get('Avg. Optimum Temp (°C)')
            #deletes blank temps
            if match.get('Avg. Optimum Temp (°C)') == '':
                temp = (match.get('Min. Temp. (°C)')+match.get('Max. Temp. (°C)'))/2
            else:
                temp = temp
            #print(temp)
            Thermobase_matches.write(f'{accn},{species},{temp}\n')
        except:
            continue


f.close()
data.close()
Thermobase_matches.close()
