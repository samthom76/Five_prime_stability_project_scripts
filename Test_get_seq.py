import glob
import os
from os.path import exists
import re

path = os.getcwd()
print(path)

files_arc = glob.glob(path + '/Archaea_filtered_genomes/' + '*embl')

print(files_arc)

files_bact = glob.glob(path + '/Bacterial_filtered_genomes/' + '*embl')

print(files_bact)

files_all = files_arc + files_bact

print(f'There are {len(files_arc)} Archeal files')
print(f'There are {len(files_bact)} Bacterial files')


files_all_short = files_all[:]
print(files_all_short)

error_file = open("errors.xls", "w")
if not exists(path + '/'+'errors'):
    os.mkdir(path + '/'+'errors')

accn_file = open('Accessions.csv', 'w')
accn_file.write(f'accn,taxon,species,genome_size_(mb),num_genes\n')

#Looping through every file
cnt_genome = 0
for fl in files_all_short:
    cnt_genome +=1
    print(f'Analysing {cnt_genome} of {len(files_all_short)}')
    #Assigning Taxon
    if re.search('Arch', fl):
        taxon = 'Archaea'
    else:
        taxon = 'Bacteria'

    #Extracting accesion number from file name (fl)
    start = fl.rfind('/')
    nd = fl.rfind('.')
    accn = fl[start + 1: nd]
    accn = accn.strip()
    print(accn)



    #Creating directory if not already present
    if not exists(path + '/'+ accn):
        os.mkdir(path + '/'+ accn)

    #Opening file
    f_open = open(fl, 'r')
    f_r = f_open.read()
    f_open.close()

    #Extracting Genome Species

    spec = f_r.find('XX\nOS')
    spec_1 = spec + 6
    spec_end = f_r.find('\n',spec_1)
    species = f_r[spec_1:spec_end]
    species = species.strip()
    species = re.sub('\s','_',species)

    #Getting and Cleaning sequence
    seq_start = f_r.find('XX\nSQ')
    #print(seq_start)
    seq_act_start = f_r.find('\n', seq_start + 10)
    seq_raw = f_r[seq_act_start:]
    seq = re.sub('[^A-Za-z]+','',seq_raw)
    seq = seq.strip()
    seq = seq.lower()
    #print(f'the length of the raw seq is {len(seq_raw)} the length of the cleaned seq is {len(seq)}')
    #print(seq)

    #Pull out genes
    gene_annot = f_r[:seq_start]
    del f_r
    gene_list = gene_annot.split('FT   gene            ')
    print(len(gene_list))
    gene_list = gene_list[1:]
    try:
        gene_list[0]
    except IndexError:
        gene_list = gene_annot.split('FT   CDS             ')
        #print(len(gene_list))
        gene_list = gene_list[1:]
        #print(gene_list[0])

    #Gene Sequence extraction:
    cnt = 0
    for gn in gene_list:
        cnt +=1
        #gene_start = gn.find('\n')
        gene_start = gn.find('/')
        top_ln = gn[:gene_start]
        #print(top_ln)
        gn_struc = re.sub('[^0-9\.,]','',top_ln)
        exon_list = gn_struc.split(',')
        #if len(exon_list) > 1:
            #print('Found a Join?')
        #print(exon_list)
        cds =''
        for x in exon_list:
                numbs = re.search('([0-9]+?)\.\.([0-9]+?)$',x)
                st = numbs.group(1)
                st = int(st)
                ndg = numbs.group(2)
                ndg = int(ndg)

        exon = seq[st-1:ndg]
        cds = cds + exon

        #Reversing Complement CDS
        if 'complement' in top_ln:
            complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
            cds = ''.join(complement.get(base, base) for base in reversed(cds))
            #print('I have been reversed')


        #Get Locus Tag
        try:
            lt = re.search('/locus_tag=\"(.*?)\"',gn)
            lt1 = lt.group(1)
        except:
            lt1 = str(cnt)

        #Get Translation Table
        try:
            tt = re.search('/transl_table=(.*?)\n',gn)
            tt1 = tt.group(1)
        except:
            continue

        #find errors in code

        #multiple of 3 long
        fail_list = []
        mult = len(cds)
        if mult % 3 != 0:
           fail_list.append('3x_')
           #print('Not a multiple of 3')

        #Do they end with a stop
        stop_codon = cds[-3:]
        stops = ["taa", "tga", "tag"]
        if not any(a in stop_codon for a in stops):
            fail_list.append('Inc_stop_')
            #print('Incorrect Stop Codon')

        #Do they start with ATG
        start_codon = cds[1:3]
        #print(start_codon)
        if start_codon != 'tg':
            fail_list.append('Inc_start_')
            #print('Incorrect start codon')

        #Do they have an internal stop codon
        n = 3
        cds_no_stop = cds[:-3]
        codon =[cds_no_stop[x:x+n] for x in range(0,len(cds_no_stop),n)]
        stops = ["taa", "tga", "tag"]
        if any(a in codon for a in stops):
            fail_list.append('Int_stop_')
            #print('Internal Stop Codon')

        #Do they have only ATGC
        allowed_cds = 'atgc'
        if all(x in allowed_cds for x in cds) is False:
            fail_list.append('Non_ATGC')
            #print('Seq contains non ATGC base')

        fail_string = "".join(fail_list)

        if len(lt1) ==0:
            lt1 = cnt

        if len(fail_list) >0:
            error_file.write(f'{accn}\t{lt1}\t{fail_string}\n')
            f_w = open(f"errors/{taxon}_{accn}_{lt1}.fta", "w")
            f_w.write(f'>{taxon}:{species}:{lt1}:{fail_string}:{accn}\n{cds}')
            f_w.close()
        else:
            f_w = open(f"{accn}/{lt1}.fta", "w")
            f_w.write(f'>{taxon}:{species}:{lt1}:{tt1}\n{cds}')
            f_w.close()
    accn_file.write(f'{accn},{taxon},{species},{(len(seq))/1000000},{len(gene_list)}\n')





error_file.close()
accn_file.close()
