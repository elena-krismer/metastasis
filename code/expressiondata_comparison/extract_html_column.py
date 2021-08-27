import pandas as pd
import urllib3
import urllib.request as urllib2
import re
import pickle
import csv


# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)


def load_obj(name):
    with open('../obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


page1 = urllib2.urlopen(
    'file:///Users/s202425/Downloads/brain_results_200605_rma.compact.print-all-prots-info.html').read()

substring_newline_html = 'td></tr>\\n'
html_colored_id_1 = '<u><font color=CC0000>'
html_colored_id_2 = '</font></u>'
html_low_expression = 'expression=<font color=33FF00>infra_expressed'
html_high_expression = 'expression=<font color=FF33FF>over_expressed'
# read html page into string/list, separat by new line
table_list = re.split(r'</td><td>|<tr><td>', str(page1))

gene_uniprot_dict = {}
id_list = []
low_expressed_proteins = []
high_expressed_proteins = []
i = 0
length_list = len(table_list) - 3
while i < length_list:
    if substring_newline_html in table_list[i]:
        if html_colored_id_1 in table_list[i + 1]:
            # extract ids which are hyperlinked highlighted
            id = re.search('<u><font color=CC0000>(.+?)</font></u>', table_list[i + 1]).group(1)
        else:
            id = table_list[i + 1]
            # id_list.append()
        id_list.append(id)
        if html_low_expression in table_list[i]:
            low_expressed_proteins.append(id)
        elif html_high_expression in table_list[i]:
            high_expressed_proteins.append(id)
    i += 1

# low expressions proteins brain
# 'P62266', 'Q02878', 'O00303', 'P61247', 'P09382', 'Q6NW38', 'P06576', 'Q6X734', 'Q15185', 'Q5VXW5',
# 'Q9UER7', 'Q9NSL4', 'P46779', 'P11142', 'P49207', 'O75106', 'Q96PE7', 'Q9UQ16', 'Q14203', 'Q8N1B4',
# 'P84090', 'O60841', 'Q00613', 'P04264', 'Q53YD7', 'Q9Y581', 'O15371', 'Q13332', 'Q9Y6S7', 'Q9BYK2',
# 'Q9UPT6', 'P63220', 'Q6IPE9', 'Q5SWX2', 'P05218', 'O15213', 'Q580W6', 'P35900', 'P17600', 'Q5JVP3',
# 'P13639', 'Q9NRH1', 'Q5D044', 'P62857', 'Q5H8X1', 'Q16775', 'O14776', 'Q4VXY7'

# high expressions proteins brain
# 'Q8IYD1', 'P13693', 'Q9H4K1', 'Q5H8X0', 'Q9UHX1', 'Q86V19', 'Q4VBY5', 'P78364', 'Q9NVI6', 'Q5TBM7',
# 'Q2XSC7', 'Q6P2G1', 'Q9UIY3', 'P19801', 'P15927', 'Q8IWL2', 'P05388', 'P62906', 'Q9BPX7', 'Q5VXU9',
# 'P63173', 'P05387', 'P35250', 'O14563', 'O14561', 'Q15056', 'Q13666', 'P02751', 'Q6DKH2', 'Q92878',
# 'Q13736', 'Q9BRP3', 'Q7Z5Y5', 'Q06609', 'Q5TEI9', 'Q13969', 'Q5TF75', 'P57678', 'Q3LRW1', 'Q6FG99',
# 'Q3LRW4', 'Q3LRW5', 'O00499', 'Q92911', 'Q5JVP0', 'Q9BSQ9', 'P41091', 'Q14489', 'P62829', 'Q8IYB7',
# 'Q8WTY3', 'P15863', 'Q9UKM9', 'Q6GMQ3', 'Q79822', 'P08708', 'Q5CAQ7', 'P22676', 'P45973', 'Q15372',
# 'Q9NR50', 'GI:46430941', 'Q9NVY0', 'Q9NXU5', 'Q2PP14', 'Q6IAD0', 'P54132', 'Q96C72', 'Q14134', 'Q9Y698',
# 'Q3MHD9', 'Q75MD9', 'P08670', 'Q99613', 'Q6ICR8', 'O15525', 'Q9NZV6', 'Q71UM5', 'Q6IPY3', 'Q9BTW9',
# 'Q5JWL1', 'P23588', 'Q5JP94', 'P60866', 'Q5JQ37', 'P08865', 'Q04637', 'Q5SZW8', 'Q8N6R0', 'P38398',
# 'O95696', 'Q8NEL2', 'O14639', 'P55010', 'Q7L684', 'P62945', 'Q03181', 'Q02543', 'Q5SXN6', 'Q9UNE7',
# 'P04183', 'Q16623', 'P49137', 'Q9NXJ5', 'Q04864', 'Q5PY61', 'Q9H1X3', 'O00255', 'Q96L34', 'P62847',
# 'Q53TT5', 'Q8NDI4', 'Q5QPL9', 'Q5QPM0', 'Q5QPM1', 'O00231', 'Q08380', 'Q8N6E1', 'Q86XS4', 'P14678',
# 'P42677', 'Q7Z6X0', 'P61601', 'Q5T6N5', 'P21673', 'P49368', 'P15880', 'Q9BXP5', 'Q4VXZ3


with open('../../data/misc/genestableID_uniprot.txt', 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        if row[1] != '':
            if row[1] in id_list:
                gene_uniprot_dict[row[1]] = row[0]

save_obj(gene_uniprot_dict, "gene_uniprot_dict")

low_geneID = [gene_uniprot_dict[k] for k in low_expressed_proteins if k in gene_uniprot_dict.keys()]
high_geneID = [gene_uniprot_dict[k] for k in high_expressed_proteins if k in gene_uniprot_dict.keys()]
brain_gene_expression = {"low_expression": low_geneID, "high_expression": high_geneID}
save_obj(brain_gene_expression, "brain_gene_expression_dict")
