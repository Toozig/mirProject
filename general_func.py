import pickle
import pandas as pd
import numpy as np


def gen_open(path):
    file = open(path, "rb")
    df = pickle.load(file)
    file.close()
    return df


def get_chrom_dict():
    return gen_open("DB/CNV_DICT(CANCER-dict_chromosom_to_CNV).pickle")

def get_grouped_CNV_info():
    df =  gen_open("DB/grouped_by_same_cnv.pickle")
    groups  = gen_open("DB/groups_of_CNV.pickle")
    return df, groups

def get_CNV_map():
    return gen_open("DB/new_CNV_map.pickle")

def get_CCLE_table():
    return gen_open("DB/CCLE_ABSOLUTE_combined.pickle")

def get_miRNA_location():
    return gen_open("DB/miRNA_loc_37.pickle")

def get_chromesome_length():
    return pd.read_csv("DB/tables/chromosome_length.csv")

def create_miRNA_location_pickle():
    df = pd.read_csv("DB/tables/38to37miRNA.csv")
    with open('miRNA_loc_37.pickle','wb') as file:
        pickle.dump(df, file)

def get_chromosomes_hist():
    return gen_open("DB/histograms.pickle")

def save_to_pickle(name, df):
    with open('%s.pickle' % name,'wb') as file:
        pickle.dump(df, file)

def get_groups():
    return pd.read_pickle("grouped_by_same_cnv.pickle")

def get_groups_dict():
    file = open("groups_dictionary.pickle", 'rb')
    return pickle.load(file)

def fix_CNV_map():
    map = get_CNV_map()

    mask =[not x.startswith("MIMAT") for x in map.to_numpy()[:,0]]
    map2= map.loc[mask]
    save_to_pickle("new_CNV_map", map2)

map = get_CNV_map()
map.to_csv("final_table.csv")