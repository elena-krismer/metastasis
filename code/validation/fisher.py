import pandas as pd

import scipy.stats as stats


"""
perform fisher exact test
compare "flux mode" (secretion/uptake) of different 
cancer cell lines and model
"""

# mcnemar code chunk maybe for later
""""
# pairwaise comparison of classifier models
def mc_nemar(classifier_1, classifier_2, model_1, model_2):
    table = confusion_matrix(classifier_1, classifier_2)
    c1, c2 = (model_2 + "- correct"), (model_2 + "- wrong")
    i1, i2 = (model_1 + "- correct"), (model_1 + "- wrong")
    table_pd = pd.DataFrame(table, columns=(c1, c2), index=(i1, i2))
    # calculate mcnemar test
    result = mcnemar(table, exact=True)
    alpha = 0.5
    print("Comparison of ", model_1,
          "and", model_2)
    print(table_pd)
    print("Statistic= %.3f, p-value= %.3f" % (result.statistic, result.pvalue))
    if result.pvalue > alpha:
        print(" -> 0 Hypothesis can not be rejected\n")
    else:
        print(" -> Reject 0 Hypothesis\n")
"""

# save cell lines with p-value lower than 0.05
significant_list = []


def fishertest(model, cellline):
    length = len(model)
    # uptake, secretion
    cellline_row = [sum(cellline), (int(length) - int(sum(cellline)))]
    model_row = [sum(model), (length - sum(model))]
    data = [cellline_row,
            model_row]
    c1, c2 = "uptake", "secretion"
    i1, i2 = cellline.name, "model"
    table_pd = pd.DataFrame(data, columns=(c1, c2), index=(i1, i2))
    odd, pvalue = stats.fisher_exact(data)
    print(table_pd)
    print("ODDs: ", odd, "p-value: ", pvalue, "\n")
    if pvalue < 0.05:
        # save significant cell lines
        significant_list.append(cellline.name)


df = pd.read_csv("../../data/wetlab_data/flux_binary_secretion_uptake.csv", sep=',', encoding="utf-8",
                 skipinitialspace=True)
# get column with model values


model = df[df.columns[126]]
i = 1
while i < 121:
    cellline = df[df.columns[i]]
    fishertest(model, cellline)
    i += 1




