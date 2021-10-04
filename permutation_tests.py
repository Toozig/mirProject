from  mir_to_gene_Network import  gene_mir_network_analysis
import numpy as np
import pandas as pd
from scipy import stats
TEST_SIZE = 0.2


def one_permutation(data, labels):
    """
    simulate one permutation
    :param data: the data for the permutation
    :param labels:  the lable of the data
    :return: the error rate for the permutation
    """
    np.random.shuffle(labels)
    _, error = gene_mir_network_analysis.get_RF_model(data, labels)
    return error


def single_tailed_test(real_error, error_lst, R):
    return  np.sum(np.asarray(error_lst) >= real_error) / R



def comp_to_norm_biostatistic(real_error, error_lst):
    mean = np.mean(error_lst)
    sd = np.std(error_lst)
    z = (real_error - mean) / sd
    return stats.norm.sf(abs(z))


def normail_dist_compare_site(real_error, error_lst):
    p = stats.ttest_1samp(error_lst, real_error)
    return p.pvalue





def permutation_test(permutation_num):
    """
    H0 - The e
    :param permutation_num:
    :return:
    """
    data , labels = gene_mir_network_analysis.classification_from_network()
    _, real_error = gene_mir_network_analysis.get_RF_model(data, labels)
    error_lst = []
    for i in range(permutation_num):
        error_lst.append(one_permutation(data, labels))
    df = pd.DataFrame()
    df["single tailed test"] = single_tailed_test(real_error, error_lst, permutation_num)
    df["compare to norm biostatistic"] =  comp_to_norm_biostatistic(real_error, error_lst)
    df["normail_dist_compare_site"] = normail_dist_compare_site(real_error, error_lst)
    print(df)
    print(error_lst)
    df.to_csv("permutation(%s)_test_P_value.csv" % permutation_num)



permutation_test(1000)