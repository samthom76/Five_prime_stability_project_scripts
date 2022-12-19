from pandas import *


#open thermobase csv
f = open("Thermobase_matches.csv", 'r')
f_data = f.read()
f.close()

#split csv by line
tb_genome_list = f_data.split('\n')
tb_genome_list = tb_genome_list[1:]


#create output file
Unique_temp = open('Unique_temp.csv', 'w')
Unique_temp.write(f'accn,temp\n')

#convert thermobase csv into list
tb_list = []
for i in tb_genome_list:
    line = i.split(',')
    try:
        tb_data = line[0:3]
    except:
        continue
    tb_list.append(tb_data)
#print(tb_list)

#open tempura csv
g = open("Tempura_matches.csv", 'r')
g_data = g.read()
g.close()

#split csv by line
tp_genome_list = g_data.split('\n')
tp_genome_list = tp_genome_list[1:]


#convert tempura csv into list
tp_list = []
for x in tp_genome_list:
    line = x.split(',')
    try:
        tp_accn = line[0]
        #remove duplicate that has a different temp in each data set so later code won't remove it
        if tp_accn == 'BA000002':
            del x
        tp_data = line[0:3]
    except:
        continue
    tp_list.append(tp_data)
#print(tp_list)

#creating combined list from both datasets and finding the unique values
comb_list = tp_list + tb_list
unique_list = [list(z) for z in set(tuple(z) for z in comb_list)]

#write to output file contaiing only accession code and temperature
for i in unique_list:
        try:
            accn = i[0]
            temp = i[2]
            Unique_temp.write(f'{accn},{temp}\n')
        except:
            continue

Unique_temp.close()
