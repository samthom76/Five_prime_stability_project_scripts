# Five_prime_stability_project_scripts
README

The aim of the project is to extract the gene and structural rna sequences for a collection of filtered
bacterial and archeal genomes and then analyse their GC/GC3 content and structural stability. For more
detail code is annotated explaining individual steps in each proccess.
Extra Step = Install Pandas module. If using pycharm this can be installed through 'preferences' if
using terminal install from these instructions: https://pandas.pydata.org/docs/getting_started/install.html
Extra Step: Install Vienna RNA package instructions can be found here: https://www.tbi.univie.ac.at/RNA/

test_get_seq.py

Aim = Extract the sequence for each gene within each genome, creating a FASTA file with their DNA 
sequence as well as other key information such as taxon, species, locus tag, length of seq(mb), position in 
gene list. Also create folders for each genome using accession numbers, that all future analysis can be
stored in. The sequences will be tested using the below
Tests = Test that each sequence has only ATGC bases, is a multiple of 3, has a correct start (ATG) 
codon, ends with a stop codon (TAA, TGA, TAG) and does not have an internal stop codon
Input files = Archeal filtered genomes, Bacterial filtered genomes
Output = Folder for each genome named by Accession number
Output = FASTA file for each gene within each genome - these are stored in the folders for each
genome named by accession. The files will be named in the format locus_tag.fta
Output = Error file that collects and removes sequences that do not pass the tests
Output = Accessions.csv - list of accession numbers along with taxon, species, genome size, number of genes
all by accession

get_trna_seq.py

Aim = Extract the sequence for each trna within each genome, creating a FASTA file with their 
sequence as well as other key information such as taxon, species, length of seq(mb), position in 
trna list. Store analysis in respective accession folders. The sequences will be tested using the below
Tests = Test that each sequence has only ATGC bases - CHECK IF THIS IS IN OTHER CODE
Input files = Archeal filtered genomes, Bacterial filtered genomes
Output = FASTA file for each trna within each genome - these are stored in the folders by accession for
each genome. The files will be named in the format locus_tag.fta
Output = trna.csv - list of accession numbers along with taxon, species, genome size, number of trna
all by accession

get_rna_seq.py

Aim = Extract the sequence for each rrna within each genome, creating a FASTA file with their 
sequence as well as other key information such as taxon, species, length of seq(mb), position in 
rrna list. Store analysis in respective accession folders. The sequences will be tested using the below
Tests = Test that each sequence has only ATGC bases - CHECK IF THIS IS IN OTHER CODE
Input files = Archeal filtered genomes, Bacterial filtered genomes
Output = FASTA file for each rrna within each genome - these are stored in the folders by accession for
each genome. The files will be named in the format locus_tag.fta
Ouput = rrna.csv - list of accession numbers along with taxon, species, genome size, number of rrna
all by accession

use_accn.py

Aim = Use the gene sequences extracted to work out the GC and GC3 content for each genome. The gene's 
GC and GC3 will also be analysed at the first 36(5') mid 36 and last 36(3') bases. A csv will also be
created showing the above results by accession number.
Input = FASTA file for each gene within each genome - these are stored in the folders for each
Input = Accession.csv
genome named by accession. The files will be named in the format locus_tag.fta
Output = GC_Analysis.CSV - csv by accesions with the variables detailed above.


use_trna.py
Aim = Use the trna sequences extracted to work out the GC content for each trna. A csv will also be
created showing the above results by accession number.
Test = Some genomes do not have their trna's annotated so no sequences will be extratced, for these
genomes 'na' will appear in the csv for their trna GC content.
Input = FASTA file for each trna within each genome - these are stored in the folders by accession for
Input = trna.csv 
each genome. The files will be named in the format locus_tag.fta
Output = GC_trna_Analysis.CSV - csv by accesions with the variables detailed above.


use_rrna.py
Aim = Use the rrna sequences extracted to work out the GC content for each rrna. A csv will also be
created showing the above results by accession number.
Test = Some genomes do not have their rrna's annotated so no sequences will be extratced, for these
genomes 'na' will appear in the csv for their rrna GC content.
Input = FASTA file for each rrna within each genome - these are stored in the folders by accession for
each genome. The files will be named in the format locus_tag.fta
Input = rrna.csv
Output = GC_rrna_Analysis.CSV - csv by accesions with the variables detailed above.

make_stability_files.py
Aim = for each genome, create three files with the sequence of each gene first 36(5') mid 36 and last 
36(3') bases. The output file will be formated as a FASTA file formated as such: >Locus tag, Accession,
gene count. Then the sequence on the line below. Then the above will be repeated for each gene in the 
genome. These files can then be used to analyse the stability of the sequences
Input = FASTA file for each gene within each genome - these are stored in the folders for each
Input = Accession.csv
Output = 3 files for each accession number/ genome in the format: accession_fp.fta, accession_tp.fta, 
accession_mid.fta

analyse_stability_files.py
Aim = Use the Vienna package to analyse the stability of the sequences collected above
Input = 3 files for each accession number/ genome in the format: accession_fp.fta, accession_tp.fta, 
accession_mid.fta
Input = Accesions.csv
Output = 3 files for each accession number/ genome in the format: accession_fpstability.fta, 
accession_tpstability.fta, accession_midstability.fta


average_stability_files.py
Aim = Calculate the average stability of the sequences for each genome from the above files and read
out to a csv in the format accesion, 5' stability, 3' stability, mid stability
Input = 3 files for each accession number/ genome in the format: accession_fpstability.fta, 
accession_tpstability.fta, accession_midstability.fta
Input = Accesions.csv
Output = stability_files.csv - contains the mean stability of the 3 variables in the format outlined
above

Tempura.py
Aim = Find as many optimal growth temperature of the bacteria/archea as possible, by matching with the csv downloaded from the tempura database
Extra Step: Download tempura csv from here: http://togodb.org/db/tempura
Input = 200617_TEMPURA.csv
Input = Accessions.csv
Output = Tempura_matches.csv

Thermobase.py
Aim = Find as many optimal growth temperature of the bacteria/archea as possible, by matching with the csv downloaded from the thermobase database
Extra Step: Download Thermobase csv from here: http://togodb.org/db/tempura
Input = ThermoBase_ver_1.0_2022.csv
Input = Accessions.csv
Output = Thermobase_matches.csv

Unique_temp.py
Aim = Match Thermobase_matches.csv and Tempura_matches.csv together so only unique temperatures are found
Input = Thermobase_matches.csv
Input = Tempura_matches.csv
Output = Unique_temp.csv

crawlerkmerforunkowntemps.py
Note: This code was written and sourced from Kate Daniels
Aim = Use CnnPOGTP to predict optimal growth temperatures for all genomes in the analysis
Extra Step: Install packages listed in code
Input = Archeal filtered genomes, Bacterial filtered genomes
Output = unknown_tempresults.txt

join_csv.py
Aim = Merge the csv's created in previous steps to create a 'final' csv for analysis in R
Input = GC_Analysis.csv
Input = GC_trna_Analysis.csv
Input =  GC_rrna_Analysis.csv
Input = GC_Analysis_final.csv
Input = Unique_temp.csv
Input = unknown_tempresults.txt
Output= analysis.csv



