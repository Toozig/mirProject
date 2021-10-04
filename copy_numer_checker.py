import copy
import pandas as pd
import gffpandas.gffpandas as gffpd
import pickle
import re
CHROMOSOME = 'Chr'
ID_REGEX = "ID=(.+?);"
MI_RNA = 'attributes'

ID = 'ID'

DIPLOID = 2

NUM_OF_CHROMOSOME = 24

CANCER_ID = "DB/CNCER_ID.csv"

SAMPLE = "sample"
ID_CCLE = "CCLE_ID"
START = 'Start'
END = 'End'
CHROMOSOME = 'Chromosome'
TOTAL_CNV = 'Modal_Total_CN'


class RangeDict(dict):
    def __getitem__(self, item):
        for key in self:
            if item.start > key.stop :
                continue
            if item.stop < key.start :
                continue
            if self.range_subset(item, key):  # todo make it more effciant
                return super().get(key)
        return DIPLOID


    @staticmethod
    def range_subset(range1, range2):
        """Whether range1 is a subset of range2."""
        return range1.start in range2 and range1[-1] in range2


def get_chromosome_CNV_map(datafram, type_set):
    """
    create a dictionary for each cancer type to map each area
    and the CNV of it
    :param datafram:
    :type DataFrame
    :return: dictionary
    :type dict
    """
    # Gets all cancer ID

    # create range dictionary for each cancer type
    overall_dict = {}
    for cancer_type in type_set:
        chromosome_list = []
        sample_df = datafram[datafram[SAMPLE] == cancer_type]
        # creates dictionary for each chromosome
        for chromosome in range(1, NUM_OF_CHROMOSOME):
            chrom_df = sample_df[sample_df[CHROMOSOME] == chromosome]
            chr_dict = RangeDict({range(chrom_df[START][idx],
                                        chrom_df[END][idx]):
                                      chrom_df[TOTAL_CNV][idx] for idx
                                  in chrom_df.index})
            chromosome_list.append(chr_dict)
        overall_dict[cancer_type] = chromosome_list

    return overall_dict


def create_CNV_table(CNV_dict, location_df, type_set):
    """

    :param CNV_dict:
    :param location_table:
    :param type_set:
    :type type_set list
    :return:
    """
    columns = copy.copy(type_set)
    columns.insert(0, MI_RNA)
    columns.insert(1, CHROMOSOME)
    table = pd.DataFrame(columns=columns)
    table.set_index(MI_RNA)
    mi_id = location_df[MI_RNA].tolist()
    for id in range(len(mi_id)):
        loc_range = range(
            location_df.loc[location_df[MI_RNA] == mi_id[id]]['start'].values[
                0],
            location_df.loc[location_df[MI_RNA] == mi_id[id]]['end'].values[
                0] + 1)
        chromosome_num = int(
            location_df.loc[location_df[MI_RNA] == mi_id[id]]['seq_id'].values[
                0])
        table.loc[id] = re.findall(ID_REGEX,mi_id[id]) + [0] *( len(CNV_dict ) +1)
        table.loc[id][CHROMOSOME] = chromosome_num
        # checks the CNV for each type of cancer.
        for cancer_type in CNV_dict.keys():
            # get's the current chromosom dict
            cancer_dict = CNV_dict[cancer_type][chromosome_num - 1]
            CNV = cancer_dict[loc_range]
            table.loc[id][cancer_type] = CNV
    with open("main_table", "wb") as handle:
        pickle.dump(table, handle)
    return table


def create_table():
    type_set = pd.read_csv('DB/CNCER_ID.csv')[ID_CCLE].tolist()
    miRNA_df = gffpd.read_gff3("DB/hsa_converted.gff").df
    file = pickle.load( open("DB/CCLE_ABSOLUTE_combined.pickle", "rb" ))
    a =  pickle.load( open("DB/CNV_DICT(CANCER-dict_chromosom_to_CNV).pickle", "rb" ))
    b = create_CNV_table(a, miRNA_df, type_set)
    b.to_excel("miRNA_CNV_table_with_chromosome.xlsx")
    attr = pd.read_csv("DB/ID_TEST")
    file.insert(0,'ID', attr["attributes"])
    file.to_excel("result.xlsx")



if __name__ == '__main__':


    create_table()