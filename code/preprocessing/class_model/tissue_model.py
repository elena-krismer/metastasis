from code.preprocessing.model_creation.model_object import model
import pandas as pd
import scipy.io

df = pd.read_csv('../../../data/gene_expression/gtex_median_tissue_tpm.txt', delimiter="\s")
geneID = df.iloc[:,0]

brain = model('brain', geneID.to_numpy(),
                           df.loc[:,'"brain"'].to_numpy())
scipy.io.savemat('../../../obj/models/brain_class.mat', mdict={'brain': brain})

breast = model('breast', geneID.to_numpy(),
                           df.loc[:,'"breast"'].to_numpy())
scipy.io.savemat('../../../obj/models/breast_class.mat', mdict={'breast': breast})

lung = model('lung', geneID.to_numpy(),
                           df.loc[:,'"lung"'].to_numpy())
scipy.io.savemat('../../../obj/models/lung_class.mat', mdict={'lung': lung})

