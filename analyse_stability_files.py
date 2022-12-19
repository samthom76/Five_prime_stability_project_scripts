import glob
import os
from os.path import exists
import re


path1 = os.getcwd()
print(path1)


f = open("Accessions.csv", 'r')
f_data = f.read()
f.close()

genome_list = f_data.split('\n')
genome_list = genome_list[1:]
print(genome_list)

#get accn for each genome
cnt_genome = 0
for i in genome_list:
    cnt_genome +=1
    line = i.split(',')
    accn = line[0]
    print(accn)
    print(f'I am analysing {cnt_genome} in {len(genome_list)}')

    #open 5'/core/3' files and force code to run in operating system to analyse seq in Vienna RNA
    for x in path1:
        if accn != '':
            fp_input = accn + '_fp.tfa'
            fp_output = accn + '/' + accn + '_fpstability.tfa'
            fp = 'RNAfold <' + fp_input + '>' + fp_output
            os.system(fp)
            tp_input = accn + '_tp.tfa'
            tp_output = accn + '/' + accn + '_tpstability.tfa'
            tp = 'RNAfold <' + tp_input + '>' + tp_output
            os.system(tp)
            mid_input = accn + '_mid.tfa'
            mid_output = accn + '/' + accn + '_midstability.tfa'
            mid = 'RNAfold <' + mid_input + '>' + mid_output
            os.system(mid)
            #remove reduntant code
            test = os.listdir(path1)
            for item in test:
                if item.endswith('.ps'):
                    os.remove(os.path.join(path1, item))


