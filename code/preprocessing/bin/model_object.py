import pickle
import pandas as pd
import numpy, scipy.io

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


# load translated expression data
df = pd.read_csv('../../../data/gene_expression/bin/lung_metastasis_quantile_normalized.csv', delimiter=",")
# extract the gene ID
geneID = df.iloc[:, 0]
# combine the values by calculating the mean
mean = df.iloc[:, 1:5].mean(axis=1)
# save expression data as class model
lung_metastasis = model('lung_metastasis', geneID.to_numpy(), mean.to_numpy())


df = pd.read_csv('../../../data/gene_expression/raw/primary_breast_cancer_raw.csv', delimiter=";")
# extract the gene ID
geneID = df.iloc[:, 0]
value = df.iloc[:, 1].to_numpy()
value = [s.replace(',', '.')for s in value]
primary_breast_cancer = model('primary_breast_cancer', geneID.to_numpy(),
                               value)


df = pd.read_csv('../../../data/gene_expression/raw/regular_lungtissue_raw.csv', delimiter=";")
# extract the gene ID
geneID = df.iloc[:, 0]
value = df.iloc[:, 1].to_numpy()
value = [s.replace(',', '.')for s in value]
print(value)
lungtissue = model('lungtissue', geneID.to_numpy(), value)


# save as mat for matlab
#scipy.io.savemat('../../obj/models/lung_metastasis.mat', mdict={'lung_metastasis': lung_metastasis})

#scipy.io.savemat('../../obj/models/primary_breast_cancer.mat', mdict={'primary_breast_cancer': primary_breast_cancer})

scipy.io.savemat('../../../obj/models/bin/lungtissue.mat',
                 mdict={'lungtissue': lungtissue})
