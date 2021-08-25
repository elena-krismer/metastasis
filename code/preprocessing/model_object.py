import pickle
import pandas as pd
import numpy, scipy.io
import csv
import gzip

# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)


def load_obj(name):
    with open('../../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# create class model with expression data and gene names
class model:
    def __init__(self, name, gene, value):
        self.name = name
        # array with gene names
        self.gene = gene
        # expression data
        self.value = value


# dict for Affymetrix Human Genome U133 Plus 2.0 Array - GeneID
plusU133_dict = {}

# create dictionary
with open('../../data/misc/geneID_PlusU133.txt', 'r') as csv_file:
   for row in csv.reader(csv_file, delimiter='\t'):
       if row[1] != '':
            plusU133_dict[row[1]] = row[0]

save_obj(plusU133_dict, "plusU133_dict")

aU133 = load_obj("aU133_dict")


# tranlate affymetrix ID to geneID
flag = False
with gzip.open('../../data/gene_expression/raw/GSE11078_breast.txt.gz', 'rt') as csv_file:
    with open("../../data/gene_expression/GSE11078.txt", "w") as output:
        for row in csv.reader(csv_file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] == "!series_matrix_table_end":
                flag = False
            if flag:
                if row[0] in plusU133_dict:
                    row[0] = plusU133_dict.get(row[0])
                    output.write(str('\t'.join(row) + '\n'))
                else:
                    output.write(str('\t'.join(row) + '\n'))
            if row[0] == "!series_matrix_table_begin":
                flag = True
    output.close()
    csv_file.close()


with open('../../data/gene_expression/GSE2603_RMA.txt', 'r') as csv_file:
    with open("../../data/gene_expression/GS2603_RMA_geneID.txt", "w") as output:
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

# load translated expression data
df = pd.read_csv('../../data/gene_expression/GSE11078.txt', delimiter="\t")
df = df[df.ID_REF.str.contains("ENSG")]
# extract the gene ID
geneID = df.iloc[:, 0]
# lung GSM279964, GSM279974, GSM279975, GSM279977, GSM279978
lung_tissue = ['GSM279964', 'GSM279974', 'GSM279975', 'GSM279977', 'GSM279978']
df_lung = df.loc[:, lung_tissue].copy()
# combine the values by calculating the mean
mean = df_lung.mean(axis=1)
# save expression data as class model
lung_metastasis = model('lung_metastasis', geneID.to_numpy(), mean.to_numpy())

# Breast Cancer Met-1/8/9/10
brain = ['GSM279958', 'GSM279969', 'GSM279971', 'GSM279972']
df_brain = df.loc[:, brain].copy()
mean = df_brain.mean(axis=1)
brain_metastasis = model('brain_metastasis', geneID.to_numpy(), mean.to_numpy())


# MDA-MB cell line
# df MDA-MB-231 cell line
df_2 = pd.read_csv('../../data/gene_expression/GS2603_RMA_geneID.txt', delimiter="\t")
print(df_2.index)
df_2 = df_2[df_2.index.str.contains("ENSG")]
pd.Index(df_2.columns, dtype=object)
df_breast= df_2.loc[:, 'GSM50017.CEL':'GSM50032.CEL'].copy()
df_breast['mean'] = df_breast.mean(axis=1)
df_breast['std'] = df_breast.std(axis=1)
MDA_MB_231 = model('MDA_MB_231', df_2.index.to_numpy(),
                           df_breast['mean'].to_numpy())

# save as mat for matlab
scipy.io.savemat('../../obj/models/lung_metastasis.mat', mdict={'lung_metastasis': lung_metastasis})

scipy.io.savemat('../../obj/models/brain_metastasis.mat', mdict={'brain_metastasis': brain_metastasis})

scipy.io.savemat('../../obj/models/MDA_MB_231.mat', mdict={'MDA_MB_231': MDA_MB_231})
