import glob
import os
from os.path import exists
import csv
from csv import DictReader
import re
import pandas as pd

path1 = os.getcwd()
print(path1)

#Defining csv's being loaded in
filename = '200617_TEMPURA.csv'
accn_file = 'Accessions.csv'

#create final output file
Tempura_matches = open('Tempura_matches.csv', 'w')
Tempura_matches.write(f"accn,species,temp\n")


#opening tempura csv and making it into list of dictionary
with open(filename, 'r') as f:
    f_dict_reader = DictReader(f)
    f_list_of_dict = list(f_dict_reader)
    f_list_of_dict = f_list_of_dict[:]
    #joining Genus and species with strain column to get complete name
    for i in f_list_of_dict:
        name = i.get('Genus_and_species')
        strain = i.get('Strain')
        i.update({'Genus_and_species': name + ' ' + strain})
        name_new = i.get('Genus_and_species')
        i.update({'Genus_and_species': name_new.replace(' ', '_')})
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
            match = next(h for h in f_list_of_dict if h['Genus_and_species'] == species)
            print(match)
            temp = match.get('Topt_average(ºC)')
            #print(temp)
            Tempura_matches.write(f'{accn},{species},{temp}\n')
        except:
            continue


f.close()
data.close()
Tempura_matches.close()
#print(len(list_of_dict))
