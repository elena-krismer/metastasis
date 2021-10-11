import pickle
import csv

# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f,  protocol = 2)

def load_obj(name ):
    with open('../obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

subsystem_dict = {}
# create dict with reactions belong to each subsystem
with open('../../data/misc/sub_reactions.txt', 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        if len(row) != 0:
            subsystem = row[0]
            reactions = row[1:]
            subsystem_dict[subsystem] = reactions

csv_file.close()
save_obj(subsystem_dict, "subsystem_dict")

# add CMP-Neu5AC biosynthesis https://www.genome.jp/module/M00922+C00128
# R00414 UDP-N-acetyl-D-glucosamine 2-epimerase
# R01804
# R01117


# Hexosamine Pathway
# Fructose 6 -phosphate -> UDP-N-acetylglucosamine MAM03111c

# MAR04300, MAR04628, MAR04631, MAR04158

#Sialic acid biosynthesis
# UDP-N-acetylglucosamine -> N-acetylneuraminic acid (Sialic acid) MAM02543c
# MAR04526, MAR04528, MAR4529, MAR04530

#Acitvation
# -> CMP- N-acetylneuraminate (CMP-sialic acid)
# MAR08377
# ALternative: MAR04531,
# - >MAR04532 -> MAR08377
# MAR04532 -> MAR04533
# MARO4532 -> MAR04534

#no activation MAM02524c N-acetyl-Dmannosamine
# no acitvation MAM02819c pyruvate

"""

'MAM01965c' glucose#
'MAM01975c' glutamine # MAR04300
'MAM01261c' acetyl coa# MAR04628
'MAM02696c' PEP
'MAM03130c' UTP
'MAM01371c' ATP
'MAM02040c' H20
'MAM01623c' CTP


#fructose 6 phosphate MAM01845c

'MAM01592c' CMP-N-acetylneuraminate
'MAM01974c' glutamate
'MAM01597c' CoA
'MAM01285c' ADP
'MAM02759c' PPi
'MAM02039c' H
'MAM02751c' Pi
'MAM03106c' UDP# MAR004526
"""