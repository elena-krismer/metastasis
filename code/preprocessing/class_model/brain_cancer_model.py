import pandas as pd
import scipy.io
import csv
from model_object import model, load_obj

aU133 = load_obj("aU133_dict")

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13276
"""
This study attempted to define the molecular characteristics of the GBM surrounding tissue. To this end, GBM tumor 
samples were obtained from 5 patients who underwent total tumor resection. Surrounding tumor mass tissue was retrieved 
in all cases from not infiltrated white matter sited at 2 cm from the macroscopic tumor border. Furthermore, control 
white matter biopsies were harvested from patients operated on for deep intracerebral cavernomas. Each sample was 
hybridized onto Affymetrix human U133 arrays. For each patient, tumor core sample and surrounding tissue were harvested 
and are identified with the same suffix number. In 2 cases, (patients 3 and 4), two tumor peripheral tissue samples were 
harvested and are identified with the same number followed by "R" (replicate).
"""
"""""
with open('../../../data/gene_expression/GSE13276_RMA.txt', 'r') as csv_file:
    with open("../../../data/gene_expression/GSE13276_RMA_geneID.txt", "w") as output:
        for row in csv.reader(csv_file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] in aU133:
                row[0] = aU133.get(row[0])
                output.write(str('\t'.join(row) + '\n'))
            else:
                output.write(str('\t'.join(row) + '\n'))
    output.close()
    csv_file.close()


df = pd.read_csv('../../../data/gene_expression/GSE13276_RMA_geneID.txt', delimiter="\t")

# glioblastoma 'GSM335185.CEL':'GSM335189.CEL'
# gbm_surrounding_tissue 'GSM335190.CEL':'GSM335196.CEL'
# white_matter 'GSM335197.CEL':'GSM335199.CEL'
pd.Index(df.columns, dtype=object)
df_glioblastoma= df.loc[:, 'GSM335185.CEL':'GSM335189.CEL'].copy()
df_glioblastoma['mean'] = df_glioblastoma.mean(axis=1)
glioblastoma = model('glioblastoma', df_glioblastoma.index.to_numpy(),
                           df_glioblastoma['mean'].to_numpy())

df_gbm_surrounding_tissue= df.loc[:, 'GSM335190.CEL':'GSM335196.CEL'].copy()
df_gbm_surrounding_tissue['mean'] = df_gbm_surrounding_tissue.mean(axis=1)
gbm_surrounding_tissue = model('gbm_surrounding_tissue', df_gbm_surrounding_tissue.index.to_numpy(),
                           df_gbm_surrounding_tissue['mean'].to_numpy())

df_white_matter= df.loc[:, 'GSM335197.CEL':'GSM335199.CEL'].copy()
df_white_matter['mean'] = df_white_matter.mean(axis=1)
white_matter = model('white_matter', df_white_matter.index.to_numpy(),
                           df_white_matter['mean'].to_numpy())

# save as mat for matlab
# issues when sampling blastoma model
scipy.io.savemat('../../../obj/models/glioblastoma_class_corrupt_GSE13276.mat', mdict={'glioblastoma': glioblastoma})

scipy.io.savemat('../../../obj/models/gbm_surrounding_tissue_class.mat', mdict={'gbm_surrounding_tissue': gbm_surrounding_tissue})

# issues when sampling model
scipy.io.savemat('../../../obj/models/white_matter_class.mat', mdict={'white_matter': white_matter})
"""



with open('../../../data/gene_expression/GSE50161_RMA.txt', 'r') as csv_file:
    with open("../../../data/gene_expression/GSE50161_RMA_geneID.txt", "w") as output:
        for row in csv.reader(csv_file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] in aU133:
                row[0] = aU133.get(row[0])
                output.write(str('\t'.join(row) + '\n'))
            else:
                output.write(str('\t'.join(row) + '\n'))
    output.close()
    csv_file.close()


df = pd.read_csv('../../../data/gene_expression/GSE50161_RMA_geneID.txt', delimiter="\t")

# ependymoma
pd.Index(df.columns, dtype=object)
df_ependymoma= df.loc[:, 'GSM1214834.CEL':'GSM1214879.CEL'].copy()
df_ependymoma['mean'] = df_ependymoma.mean(axis=1)
ependymoma = model('ependymoma', df_ependymoma.index.to_numpy(),
                           df_ependymoma['mean'].to_numpy())

pd.Index(df.columns, dtype=object)
df_glioblastoma= df.loc[:, 'GSM1214880.CEL':'GSM1214913.CEL'].copy()
df_glioblastoma['mean'] = df_glioblastoma.mean(axis=1)
glioblastoma = model('glioblastoma', df_glioblastoma.index.to_numpy(),
                           df_glioblastoma['mean'].to_numpy())

df_med = df.loc[:, 'GSM1214915.CEL':'GSM1214935.CEL'].copy()
df_med['mean'] = df_med.mean(axis=1)
medullablastoma = model('medullablastoma', df_med.index.to_numpy(),
                           df_med['mean'].to_numpy())

df_pa= df.loc[:, 'GSM1214949.CEL':'GSM1214963.CEL'].copy()
df_pa['mean'] = df_pa.mean(axis=1)
pilocyticastrocytoma = model('pilocyticastrocytoma', df_pa.to_numpy(),
                           df_pa['mean'].to_numpy())

# save as mat for matlab
scipy.io.savemat('../../../obj/models/ependymoma_class.mat', mdict={'ependymoma': ependymoma})
scipy.io.savemat('../../../obj/models/glioblastoma_class.mat', mdict={'glioblastoma': glioblastoma})
scipy.io.savemat('../../../obj/models/pilocyticastrocytoma_class.mat', mdict={'pilocyticastrocytoma': pilocyticastrocytoma})
scipy.io.savemat('../../../obj/models/medullablastoma_class.mat', mdict={'medullablastoma':medullablastoma})