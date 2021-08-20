import csv 
import pandas as pd

def subsystem_df(path, subsystem_list, cell_type):
    with open(path, 'r') as file:
        min_cell_type = 'min' + cell_type
        max_cell_type = 'max' + cell_type
        df = pd.DataFrame(columns=['reaction', min_cell_type, 'min_breast',
                                   max_cell_type, 'max_breast', 'prop1', 'prop2'], index=[1, 2])
        for row in csv.reader(file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] in subsystem_list:
                reaction = row[0]
                min1, min2 = row[15], row[16]
                max1, max2 = row[17], row[18]
                prop1_value, prop2_value = row[19], row[20]
                print(reaction, min1, min2, max1, max2, prop1_value, prop2_value)
                df2 = pd.Series([reaction, min1, min2, max1, max2, prop1_value,
                prop2_value], index=df.columns)
                df = df.append(df2, ignore_index=True)
        file.close()
        df.dropna(subset=["reaction"], inplace=True)
        return df

def dataframe_merge(df1, df2):
  outer_merged = pd.merge(df1, df2, how="outer", on=["reaction"])
  return outer_merged

def subset_merge(path1, path2, subsystem_list,celltype1, celltype2):
  df1 = subsystem_df(path1, subsystem_list, celltype1)
  df2 = subsystem_df(path2, subsystem_list, celltype2)
  df = dataframe_merge(df1, df2)
  return df
 

