import glob
import os
from os.path import exists
import csv
from csv import DictReader
import re


path1 = os.getcwd()
print(path1)

f = open("Accessions.csv", 'r')
f_data = f.read()
f.close()

genome_list = f_data.split('\n')
genome_list = genome_list[1:]
print(genome_list)


cnt_genome = 0
#create file name per gene for 5',core and 3' (12 codons for each)
for i in genome_list:
    cnt_genome +=1
    line = i.split(',')
    accn = line[0]
    Accn_fp = open(f"{accn}" + '_fp' + '.tfa', 'w')
    Accn_fp.write('')
    Accn_tp = open(f"{accn}" + '_tp' + '.tfa', 'w')
    Accn_tp.write('')
    Accn_mid = open(f"{accn}" + '_mid' + '.tfa', 'w')
    Accn_mid.write('')
    files_tfa = glob.glob(path1 + '/' + accn + '/' + '*fta')
    print(accn)
    print(len((files_tfa)))
    print(f'Analysing {cnt_genome} of {len(genome_list)}')
    files_tfa_short = files_tfa[:5]


    cnt_gene = 0
    #get seq for each section of each gene and write new files
    for gn in files_tfa_short:
        cnt_gene += 1
        gn = open(gn, 'r')
        seq = gn.read()
        gn.close()
        seq = seq.split('\n')
        lt_get = seq[:1]
        lt_get = ''.join(lt_get)
        lt_get = lt_get.split(':')
        lt = lt_get[2]
        seq = seq[1:]
        seq = ''.join(seq)
        seq_five_prime = seq[0:36]
        seq_three_prime = seq[-36:]
        mid_point = len(seq) / 2
        mid_point = int(mid_point)
        if len(seq) % 2 == 0:
            mid_point = int(mid_point)
            seq_mid = seq[(mid_point - 18):(mid_point + 18)]
        else:
            mid_point = int(mid_point - .5)
            seq_mid = seq[(mid_point - 18):(mid_point + 18)]
        if accn != '':
            taxon = line[1]
            species = line[2]
            genome_size_mb = line[3]
            num_genes = line[4]
        if accn == '':
            del i
        else:
            Accn_fp.write(f'>{lt},{accn},{cnt_gene}\n{seq_five_prime}\n')
            Accn_tp.write(f'>{lt},{accn},{cnt_gene}\n{seq_three_prime}\n')
            Accn_mid.write(f'>{lt},{accn},{cnt_gene}\n{seq_mid}\n')

#remove redunant files
test = os.listdir(path1)
for item in test:
    if item == ('_mid.tfa'):
        os.remove(os.path.join(path1, item))


    if item == ('_fp.tfa'):
        os.remove(os.path.join(path1, item))


    if item == ('_tp.tfa'):
        os.remove(os.path.join(path1, item))





