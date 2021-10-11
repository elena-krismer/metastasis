import pandas as pd
import pickle
import numpy as np
import warnings

def load_obj(name):
    with open('../../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

subsystem_dict = load_obj("subsystem_dict")

breast_df = pd.read_csv("../../results/sampling/sampling_breast_cancer.csv")
breast_reactions = breast_df.columns.to_list()

brain_df = pd.read_csv("../../results/sampling/sampling_brain_metastasis.csv")
brain_reactions = brain_df.columns.to_list()

lung_df = pd.read_csv("../../results/sampling/sampling_lung_metastasis.csv")
lung_reactions = lung_df.columns.to_list()

print(brain_reactions)
merged_list = breast_reactions + lung_reactions + brain_reactions
common_reactions = list(set(merged_list))
print(common_reactions)


subsystem_list = [*subsystem_dict]
used_reactions = []

df = pd.read_csv("../../results/sampling_statisticalcomparison/lung_brain_breast_ratio.tsv", delimiter= "\t")

with open("../../results/sampling_statisticalcomparison/lung_brain_breast_ratio.tsv", "w") as output:
    output.write("reaction\tminFlux_lung\tminFLux_brain\t;minFlux_breast\tmaxFlux_lung\tmaxFLux_brain\tmaxFlux_breast\tmed_lung\tmed_brain\tmed_breast\tratio_lung_breast\tratio_brain_breast\n")
    for subsystem in subsystem_list:
        reaction_in_subsystem = subsystem_dict.get(subsystem)
        output.write(str(subsystem + '\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx'+'\n'))
        for reaction in reaction_in_subsystem:
            if reaction in common_reactions:
                used_reactions.append(reaction)
                if reaction in lung_reactions:
                    lung_sample = lung_df[reaction].to_numpy()
                    min_lung, med_lung, max_lung = np.percentile(lung_sample, 25),np.percentile(lung_sample, 50),\
                                               np.percentile(lung_sample, 75)
                else:
                    min_lung, med_lung, max_lung = "NA", "NA", "NA"

                if reaction in brain_reactions:
                    brain_sample = brain_df[reaction].to_numpy()
                    min_brain, med_brain, max_brain =  np.percentile(brain_sample, 25),np.percentile(brain_sample, 50), \
                                                   np.percentile(brain_sample, 75)
                else:
                    min_brain, med_brain, max_brain = "NA", "NA", "NA"

                if reaction in breast_reactions:
                    breast_sample = breast_df[reaction].to_numpy()
                    min_breast, med_breast, max_breast = np.percentile(breast_sample, 25),  np.percentile(breast_sample, 50),\
                                                     np.percentile(breast_sample, 75)
                else:
                    min_breast, med_breast, max_breast = "NA", "NA", "NA"

                if reaction in breast_reactions and reaction in lung_reactions:
                    with np.errstate(divide='ignore'):
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            ratio_lung_breast = med_lung.astype("float")/med_breast.astype("float")
                else:
                    ratio_lung_breast = "NA"

                if reaction in breast_reactions and reaction in brain_reactions:
                    with np.errstate(divide='ignore'):
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            ratio_brain_breast = med_brain.astype("float")/med_breast.astype("float")
                else:
                    ratio_brain_breast = "NA"

                row_text = [reaction, min_lung, min_brain,min_breast , max_lung,
                            max_brain, max_breast, med_lung,med_brain, med_breast, ratio_lung_breast,
                            ratio_brain_breast]
                row_text = [str(a) for a in row_text]
                output.write(str('\t'.join(row_text) + '\n'))
            else:
                continue
    misc_reactions = list(set(common_reactions) - set(used_reactions))
    for reaction in misc_reactions:
        if reaction in lung_reactions:
            lung_sample = lung_df[reaction].to_numpy()
            min_lung, med_lung, max_lung = np.percentile(lung_sample, 25), np.percentile(lung_sample, 50), \
                                           np.percentile(lung_sample, 75)
        else:
            min_lung, med_lung, max_lung = "NA", "NA", "NA"

        if reaction in brain_reactions:
            brain_sample = brain_df[reaction].to_numpy()
            min_brain, med_brain, max_brain = np.percentile(brain_sample, 25), np.percentile(brain_sample, 50), \
                                              np.percentile(brain_sample, 75)
        else:
            min_brain, med_brain, max_brain = "NA", "NA", "NA"

        if reaction in breast_reactions:
            breast_sample = breast_df[reaction].to_numpy()
            min_breast, med_breast, max_breast = np.percentile(breast_sample, 25), np.percentile(breast_sample, 50), \
                                                 np.percentile(breast_sample, 75)
        else:
            min_breast, med_breast, max_breast = "NA", "NA", "NA"

        if reaction in breast_reactions and reaction in lung_reactions:
            with np.errstate(divide='ignore'):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    ratio_lung_breast = med_lung.astype("float") / med_breast.astype("float")
        else:
            ratio_lung_breast = "NA"

        if reaction in breast_reactions and reaction in brain_reactions:
            with np.errstate(divide='ignore'):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    ratio_brain_breast = med_brain.astype("float") / med_breast.astype("float")
        else:
            ratio_brain_breast = "NA"

        row_text = [reaction, min_lung, min_brain, min_breast, max_lung,
                    max_brain, max_breast, med_lung, med_brain, med_breast, ratio_lung_breast,
                    ratio_brain_breast]
        row_text = [str(a) for a in row_text]
        output.write(str('\t'.join(row_text) + '\n'))




output.close()

