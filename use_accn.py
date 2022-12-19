import glob
import os
from os.path import exists
import csv
from csv import DictReader
import re
from collections import Counter

#Procedure to analyse GC/GC3 of each gene
def GC(seq):
    #Global GC procedure
    seq = ''.join(seq)
    counts = Counter(seq)
    a = counts['a']
    t = counts['t']
    g = counts['g']
    c = counts['c']
    Global_dict.update(
        {'a': (Global_dict.get('a') + a), 't': (Global_dict.get('t') + t), 'g': (Global_dict.get('g') + g),
         'c': (Global_dict.get('c') + c)})
    #5' GC procedure
    seq_five_prime = seq[0:36]
    counts_5_prime = Counter(seq_five_prime)
    a_5 = counts_5_prime['a']
    t_5 = counts_5_prime['t']
    g_5 = counts_5_prime['g']
    c_5 = counts_5_prime['c']
    Five_prime_dict.update({'a': (Five_prime_dict.get('a') + a_5), 't': (Five_prime_dict.get('t') + t_5),
                            'g': (Five_prime_dict.get('g') + g_5), 'c': (Five_prime_dict.get('c') + c_5)})
    #3' GC procedure
    seq_three_prime = seq[-36:]
    counts_3_prime = Counter(seq_three_prime)
    a_3 = counts_3_prime['a']
    t_3 = counts_3_prime['t']
    g_3 = counts_3_prime['g']
    c_3 = counts_3_prime['c']
    Three_prime_dict.update({'a': (Three_prime_dict.get('a') + a_3), 't': (Three_prime_dict.get('t') + t_3),
                             'g': (Three_prime_dict.get('g') + g_3), 'c': (Three_prime_dict.get('c') + c_3)})
    #Core GC procedure
    mid_point = len(seq) / 2
    mid_point = int(mid_point)
    if len(seq) % 2 == 0:
        mid_point = int(mid_point)
        seq_mid = seq[(mid_point - 18):(mid_point + 18)]
    else:
        mid_point = int(mid_point - .5)
        seq_mid = seq[(mid_point - 18):(mid_point + 18)]
    counts_mid = Counter(seq_mid)
    a_mid = counts_mid['a']
    t_mid = counts_mid['t']
    g_mid = counts_mid['g']
    c_mid = counts_mid['c']
    Mid_dict.update(
        {'a': (Mid_dict.get('a') + a_mid), 't': (Mid_dict.get('t') + t_mid), 'g': (Mid_dict.get('g') + g_mid),
         'c': (Mid_dict.get('c') + c_mid)})
    seq_3 = seq[2:]
    #Gloabl GC3 procedure
    Third_base_seq = seq_3[::3]
    counts_third_base = Counter(Third_base_seq)
    a_third_base = counts_third_base['a']
    t_third_base = counts_third_base['t']
    g_third_base = counts_third_base['g']
    c_third_base = counts_third_base['c']
    Third_base_dict.update(
        {'a': (Third_base_dict.get('a') + a_third_base), 't': (Third_base_dict.get('t') + t_third_base),
         'g': (Third_base_dict.get('g') + g_third_base), 'c': (Third_base_dict.get('c') + c_third_base)})
    #5' GC3 procedure
    GC3_five_prime = seq_five_prime[2:]
    GC3_five_prime = GC3_five_prime[::3]
    counts_GC3_5 = Counter(GC3_five_prime)
    a_GC3_5 = counts_GC3_5['a']
    t_GC3_5 = counts_GC3_5['t']
    g_GC3_5 = counts_GC3_5['g']
    c_GC3_5 = counts_GC3_5['c']
    GC3_5_prime_dict.update({'a': (GC3_5_prime_dict.get('a') + a_GC3_5), 't': (GC3_5_prime_dict.get('t') + t_GC3_5),
                             'g': (GC3_5_prime_dict.get('g') + g_GC3_5), 'c': (GC3_5_prime_dict.get('c') + c_GC3_5)})
    #3' GC3 procedure
    GC3_three_prime = seq_three_prime[2:]
    GC3_three_prime = GC3_three_prime[::3]
    counts_GC3_3 = Counter(GC3_three_prime)
    a_GC3_3 = counts_GC3_3['a']
    t_GC3_3 = counts_GC3_3['t']
    g_GC3_3 = counts_GC3_3['g']
    c_GC3_3 = counts_GC3_3['c']
    GC3_3_prime_dict.update({'a': (GC3_3_prime_dict.get('a') + a_GC3_3), 't': (GC3_3_prime_dict.get('t') + t_GC3_3),
                             'g': (GC3_3_prime_dict.get('g') + g_GC3_3), 'c': (GC3_3_prime_dict.get('c') + c_GC3_3)})
    #Core GC3 procedure
    GC3_seq_mid = seq_mid[2:]
    GC3_seq_mid = GC3_seq_mid[::3]
    counts_GC3_mid = Counter(GC3_seq_mid)
    a_GC3_mid = counts_GC3_mid['a']
    t_GC3_mid = counts_GC3_mid['t']
    g_GC3_mid = counts_GC3_mid['g']
    c_GC3_mid = counts_GC3_mid['c']
    GC3_Mid_dict.update({'a': (GC3_Mid_dict.get('a') + a_GC3_mid), 't': (GC3_Mid_dict.get('t') + t_GC3_mid),
                         'g': (GC3_Mid_dict.get('g') + g_GC3_mid), 'c': (GC3_Mid_dict.get('c') + c_GC3_mid)})


path = os.getcwd()
print(path)

f = open("Accessions.csv", 'r')
f_data = f.read()
f.close()

genome_list = f_data.split('\n')
genome_list = genome_list[1:]

#create final output file
GC_Analysis = open('GC_Analysis.csv', 'w')
GC_Analysis.write(
    f"accn,taxon,species,genome_size_(mb),num_genes,GC_ratio,5'_gc_ratio,3'_gc_ratio,mid_gc_ratio,gc3_ratio,5'_gc3_ratio,3'_gc3_ratio,mid_gc3_ratio\n")

cnt_genome = 0
#opening list of genomes
for i in genome_list:
    cnt_genome +=1
    line = i.split(',')
    accn = line[0]
    files_tfa = glob.glob(path + '/' + accn + '/' + '*fta')
    print(accn)
    print(len((files_tfa)))
    print(f'Analysing {cnt_genome} of {len(genome_list)}')
    files_tfa_short = files_tfa[:2]

    #Create dictionary for each variable
    Global_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    Five_prime_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    Three_prime_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    Mid_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    Third_base_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    GC3_5_prime_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    GC3_3_prime_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}
    GC3_Mid_dict = {'a': 0, 't': 0, 'g': 0, 'c': 0}

    #open/read each gene seq and run procedure
    for gn in files_tfa_short:
        gn = open(gn, 'r')
        seq = gn.read()
        gn.close()
        seq = seq.split('\n')
        seq = seq[1:]
        GC(seq)

    #calculate ratios from above dictionary then read to output csv
    if accn != '':
        taxon = line[1]
        species = line[2]
        genome_size_mb = line[3]
        num_genes = line[4]
        Global_GC = Global_dict['g'] + Global_dict['c']
        Global_GC_Ratio = Global_GC / (Global_GC + Global_dict['a'] + Global_dict['t'])
        #print(f"The global GC content of {accn} is {Global_GC_Ratio}")
        Five_prime_GC = Five_prime_dict['g'] + Five_prime_dict['c']
        Five_prime_GC_Ratio = Five_prime_GC / (Five_prime_GC + Five_prime_dict['a'] + Five_prime_dict['t'])
        #print(f"The 5' GC content of {accn} is {Five_prime_GC_Ratio}")
        Three_prime_GC = Three_prime_dict['g'] + Three_prime_dict['c']
        Three_prime_GC_Ratio = Three_prime_GC / (Three_prime_GC + Three_prime_dict['a'] + Three_prime_dict['t'])
        #print(f"The 3' GC content of {accn} is {Three_prime_GC_Ratio}")
        Mid_point_GC = Mid_dict['g'] + Mid_dict['c']
        Mid_point_GC_Ratio = Mid_point_GC / (Mid_point_GC + Mid_dict['a'] + Mid_dict['t'])
        #print(f"The Mid Point GC content of {accn} is {Mid_point_GC_Ratio}")
        GC3 = Third_base_dict['g'] + Third_base_dict['c']
        GC3_Ratio = GC3 / (GC3 + Third_base_dict['a'] + Third_base_dict['t'])
        #print(f"The GC3 content of {accn} is {GC3_Ratio}")
        GC3_5 = GC3_5_prime_dict['g'] + GC3_5_prime_dict['c']
        GC3_5_Ratio = GC3_5 / (GC3_5 + GC3_5_prime_dict['a'] + GC3_5_prime_dict['t'])
        #print(f"The 5' GC3 content of {accn} is {GC3_5_Ratio}")
        GC3_3 = GC3_3_prime_dict['g'] + GC3_3_prime_dict['c']
        GC3_3_Ratio = GC3_3 / (GC3_3 + GC3_3_prime_dict['a'] + GC3_3_prime_dict['t'])
        #print(f"The 3' GC3 content of {accn} is {GC3_3_Ratio}")
        GC3_Mid = GC3_Mid_dict['g'] + GC3_Mid_dict['c']
        GC3_Mid_Ratio = GC3_Mid / (GC3_Mid + GC3_Mid_dict['a'] + GC3_Mid_dict['t'])
        #print(f"The Mid GC3 content of {accn} is {GC3_Mid_Ratio}")

    if accn == '':
        del i
    else:
        GC_Analysis.write(
            f'{accn},{taxon},{species},{genome_size_mb},{len(files_tfa)},{Global_GC_Ratio},{Five_prime_GC_Ratio},{Three_prime_GC_Ratio},{Mid_point_GC_Ratio},{GC3_Ratio},{GC3_5_Ratio},{GC3_3_Ratio},{GC3_Mid_Ratio}\n')

GC_Analysis.close()

