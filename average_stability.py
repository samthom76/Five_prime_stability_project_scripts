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

#create output file
stab_file = open('stability_files.csv', 'w')
stab_file.write(f"accn,fp_stab,tp_stab,mid_stab\n")

#get accn
cnt_genome = 0
for i in genome_list:
    cnt_genome +=1
    line = i.split(',')
    accn = line[0]
    print(accn)
    #get trna stab files
    if accn != '':
        files_fp = glob.glob(path1 + '/' + accn + '/' + accn +'_fpstability.tfa')
        files_tp = glob.glob(path1 + '/' + accn + '/' + accn + '_tpstability.tfa')
        files_mid = glob.glob(path1 + '/' + accn + '/' + accn + '_midstability.tfa')

        #find gibbs free energy in 5' files
        for x in files_fp:
            files_fp = open(x, 'r')
            files_info_fp = files_fp.read()
            files_info_fp = files_info_fp.split(')\n')
            #print(files_info)
            fp_list = []
            for y in files_info_fp:
                y = y[-10:]
                #print(y)
                numbs = re.findall(r'\d+\.*\d*', y)
                numbs = ''.join(numbs)
                numbs = numbs.strip()
                if numbs != '':
                    numbs = float(numbs)
                    fp_list.append(numbs)
            #average gibbs free energy per genome
            fp_mean = sum(fp_list)/len(fp_list)
            fp_mean = 0 - fp_mean
            print(fp_mean)

        #find gibbs free energy in 3' files
        for x in files_tp:
            files_tp = open(x, 'r')
            files_info_tp = files_tp.read()
            files_info_tp = files_info_tp.split(')\n')
            #print(files_info)
            tp_list = []
            for y in files_info_tp:
                z = y[-10:]
                #print(y)
                numbs_tp = re.findall(r'\d+\.*\d*', z)
                numbs_tp = ''.join(numbs_tp)
                numbs_tp = numbs_tp.strip()
                if numbs_tp != '':
                    numbs_tp = float(numbs_tp)
                    tp_list.append(numbs_tp)
            #average gibbs free energy per genome
            tp_mean = sum(tp_list)/len(tp_list)
            tp_mean = 0 - tp_mean
            print(tp_mean)

        #find gibbs free energy in core files
        for x in files_mid:
            files_mid = open(x, 'r')
            files_info_mid = files_mid.read()
            files_info_mid = files_info_mid.split(')\n')
            #print(files_info)
            mid_list = []
            for y in files_info_mid:
                p = y[-10:]
                numbs_mid = re.findall(r'\d+\.*\d*', p)
                numbs_mid = ''.join(numbs_mid)
                numbs_mid = numbs_mid.strip()
                if numbs_mid != '':
                    numbs_mid = float(numbs_mid)
                    mid_list.append(numbs_mid)
            #average gibbs free energy per genome
            mid_mean = sum(mid_list)/len(mid_list)
            mid_mean = 0 - mid_mean
            print(mid_mean)

    if accn == '':
        del i
    else:
        stab_file.write(f"{accn},{fp_mean},{tp_mean},{mid_mean}\n")



files_fp.close()
files_tp.close()
files_mid.close()
stab_file.close()
