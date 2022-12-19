import glob
import os
from os.path import exists
import csv
from csv import DictReader
import re
from collections import Counter

#Procedure to analyse GC of each rrna
def GC(seq):
    seq = ''.join(seq)
    counts = Counter(seq)
    a = counts['a']
    t = counts['t']
    g = counts['g']
    c = counts['c']
    Global_dict.update(
        {'a': (Global_dict.get('a') + a), 't': (Global_dict.get('t') + t), 'g': (Global_dict.get('g') + g),
         'c': (Global_dict.get('c') + c)})

path = os.getcwd()
print(path)

f = open("rrna.csv", 'r')
f_data = f.read()
f.close()

genome_list = f_data.split('\n')
genome_list = genome_list[1:]

#create output file
GC_rrna_Analysis = open('GC_rrna_analysis.csv', 'w')
GC_rrna_Analysis.write(
    f"accn,num_rrna,rrna_GC_ratio\n")

#opening list of genomes
cnt_genome = 0
for i in genome_list:
    cnt_genome +=1
    line = i.split(',')
    accn = line[0]
    files_tfa = glob.glob(path + '/' + accn + '/' + '*rrna')
    print(accn)
    print(len((files_tfa)))
    print(f'Analysing {cnt_genome} of {len(genome_list)}')
    files_tfa_short = files_tfa[:]

    Global_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}

    #run procedure
    for gn in files_tfa_short:
        gn = open(gn, 'r')
        seq = gn.read()
        gn.close()
        seq = seq.split('\n')
        seq = seq[1:]
        GC(seq)


    #make Global GC na if no trna seq have been extracted
    if accn != '':
        if len(files_tfa) > 0:
            taxon = line[1]
            species = line[2]
            genome_size_mb = line[3]
            num_rrna = line[4]
            Global_GC = Global_dict['g'] + Global_dict['c']
            Global_GC_Ratio = Global_GC / (Global_GC + Global_dict['a'] + Global_dict['t'])
        else:
            Global_GC_Ratio = 'na'


    if accn == '':
        del i
    elif len(files_tfa) == 0:
        na = 'na'
        GC_rrna_Analysis.write(
            f'{accn},{na},{Global_GC_Ratio}\n')
    else:
        GC_rrna_Analysis.write(
            f'{accn},{len(files_tfa)},{Global_GC_Ratio}\n')

GC_rrna_Analysis.close()
