import pickle

def load_obj():
    with open('/Users/s202425/Documents/GitHub/metastasis/obj/dicts/brain_protein_expression_dict.pkl', 'rb') as f:
        return pickle.load(f)

dict = load_obj()

cmas_gene = "ENSG00000111726"

low_expression_list = dict.get("low_expression")
high_expression_list = dict.get( "high_expression")

if cmas_gene in high_expression_list:
    print("CMAS highly expressed")
print(high_expression_list)

if cmas_gene in low_expression_list:
    print("CMAS lowly expressed")

print(dict.keys())