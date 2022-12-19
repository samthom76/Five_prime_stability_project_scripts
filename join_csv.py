import pandas as pd

#process here is repeated for each merge - annotations are for the frist time that step is shown
#read in 2 data sets being merged
data_1 = pd.read_csv('GC_Analysis.csv')
data_2 = pd.read_csv('GC_trna_analysis.csv')

#merge data set on accession, join is inner to only get matching values from both data sets
GC_Analysis_trna = pd.merge(data_1, data_2, on='accn', how= 'inner')


#writing output to csv
GC_Analysis_trna.to_csv('GC_Analysis_trna.csv')

data_3 = pd.read_csv('GC_Analysis_trna.csv')
data_4 = pd.read_csv('GC_rrna_analysis.csv')

GC_Analysis_total = pd.merge(data_3, data_4, on='accn', how= 'inner')

#removing extra collumn created in this process
GC_Analysis_total = GC_Analysis_total.iloc[: , 1:]

GC_Analysis_total.to_csv('GC_Analysis_total.csv')

data_5 = pd.read_csv('GC_Analysis_total.csv')
data_6 = pd.read_csv('stability_files.csv')

GC_Analysis_final = pd.merge(data_5, data_6, on='accn', how= 'inner')

GC_Analysis_final = GC_Analysis_final.iloc[: , 1:]

GC_Analysis_final.to_csv('GC_Analysis_final.csv')

data_7 = pd.read_csv('GC_Analysis_final.csv')
data_8 = pd.read_csv('Unique_temp.csv')

Final_analysis = pd.merge(data_7, data_8, on='accn', how= 'outer')

Final_analysis = Final_analysis.iloc[: , 1:]

Final_analysis['temp'].fillna('na', inplace = True)

Final_analysis.to_csv('Final_analysis.csv')

data_9 = pd.read_csv('Final_analysis.csv')
data_10 = pd.read_csv('unknown_tempresults.txt')

Analysis = pd.merge(data_9, data_10, on='accn', how= 'inner')

Analysis = Analysis.iloc[: , 1:]

Analysis.to_csv('analysis.csv')





