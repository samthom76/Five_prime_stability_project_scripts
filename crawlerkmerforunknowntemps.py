import time
from selenium import webdriver
import glob
import os
import re
import pyperclip
from selenium.webdriver.common.by import By
from os.path import exists
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.action_chains import ActionChains
totalaccnlist = []
path1 = os.getcwd()


#creating file for if predictions fail
failedpredictionsfile = path1 + '/temperatures/failedpredictions.txt'
if exists(failedpredictionsfile):
    os.remove(failedpredictionsfile)
failed_predictions = open(failedpredictionsfile, 'a')

#open master csv file
info_in = open('Final_analysis.csv', 'r')

info_contents = info_in.read()
info_lines = info_contents.split('\n')
info_lines_2 = info_lines[1:-1]
info_in.close()
for x in info_lines_2:
    y = x.split(',')
    accn = y[0]
    totalaccnlist.append(accn)


#open file containing known temperatures
megatempfile = path1 + '/temperatures/megatemp.txt'
megatemp_in = open(megatempfile, 'r')
megatemp_contents = megatemp_in.read()
megatemp_lines = megatemp_contents.split('\n')
megatemp_lines_2 = megatemp_lines[1:-1]
megatemp_in.close()
knowntempaccnlist = []


for x in megatemp_lines_2:
    i = x.split('\t')
    knowntempaccn = i[0]
    knowntempopT = i[1]
    knowntempaccnlist.append(knowntempaccn)

unknowntemplist = set(totalaccnlist).difference(knowntempaccnlist)

# set up file for predicted temperatures
unknowntemppredictionsfile = path1 + '/temperatures/unknown_tempresults.txt'

#remove duplicates and set genomes that need to be predicted
if exists(unknowntemppredictionsfile):
    os.remove(unknowntemppredictionsfile)
unknown_predictions = open(unknowntemppredictionsfile, 'a')
unknown_predictions.write('accn\tpredictedtemp\n')

path1 = os.getcwd()

#opening bacterial and archael genome files
all_gen = []
bact_gen = glob.glob(path1 + '/bacterial_filtered_genomes/' + '*.embl')
arc_gen = glob.glob(path1 + '/Archaea_filtered_genomes/' + '*.embl')
all_gen = bact_gen + arc_gen
all_gen_short = all_gen[:]
M = 0
print('analysing ' + str(len(all_gen_short)) + ' genomes')

#running through genomes extracting full seqeunce (introns and exons)
for f in all_gen_short:
    print(f)
    M = M + 1
    step = M/len(all_gen_short)
    print(M, 'out of', len(all_gen_short))
    f_o = open(f, "r")
    f_r = f_o.read()
    f_o.close()
    position = f_r.find('XX\nSQ')
    position2 = f_r.find('\n', position + 7)
    rawseq = f_r[position2:]
    genome = re.sub('[^a-zA-Z]+', '', rawseq)
    genome = genome.lower()
    genome = str(genome)
    firstpos = f.rfind('/')
    lastpos = f.rfind('.')
    accn = f[firstpos + 1:lastpos]
    totalaccnlist.append(accn)
    outstuff = '>' + accn + '\n' + genome

    #put sequence into CnnPOGTP website and extarct result
    if accn in unknowntemplist:
        try:
            pyperclip.copy(outstuff)

            web = webdriver.Chrome()
            url = 'http://www.orgene.net/CnnPOGTP/'
            web.get(url)
            time.sleep(2)
            actions = ActionChains(web)

            geno = web.find_element(By.XPATH, '/html/body/div[2]/div/div/div[1]/div/div/form[1]/textarea')
            geno.click()
            actions.key_down(Keys.COMMAND).perform()
            actions.send_keys("v").perform()
            actions.key_up(Keys.COMMAND).perform()
            web.implicitly_wait(10)
            RadioButtonPeriod = web.find_element(By.XPATH, '/html/body/div[2]/div/div/div[1]/div/div/form[1]/input').click()

            web.implicitly_wait(500)
            get_confirmation_text_dix_text = web.find_element(By.CSS_SELECTOR, '.Introduction_font3')
            temp = get_confirmation_text_dix_text.text
        except BaseException as error:
            print('the code failed here', M, 'out of', len(all_gen_short), accn)
            failed_predictions.write(accn + '\n')
        #write succesful result to output predicted temperature file
        else:
            if temp.endswith('℃'):
                print('this worked', M, 'out of', len(all_gen_short), accn)
                unknown_predictions.write(accn + '\t' + temp + '\n')
            else:
                print('THIS DIDNT WORK')


unknown_predictions.close()
failed_predictions.close()
