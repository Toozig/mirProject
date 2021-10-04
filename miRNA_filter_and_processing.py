import general_func as gf
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from yellowbrick.cluster import KElbowVisualizer
import matplotlib.pyplot as plt


def statistic_filter(df: pd.DataFrame, threshold: float, statistic_func)-> pd.DataFrame:
    """
    This function filters a given CNV table according to a given threshold
    :param df: a miRNA / cancer CNV table
    :param threshold: the threshold of the variance to filter
    :return: a filtered DF
    """

    stat_calculation = statistic_func(1).to_numpy()
    mask = stat_calculation > threshold
    idx = np.argwhere(mask).flatten()
    return df.iloc[idx]


def variance_filtering(df: pd.DataFrame, sd_threshold:float):
    return statistic_filter(df, sd_threshold, df.drop(['ID'], axis=1).var)


def sd_filtering(df: pd.DataFrame, sd_threshold:float):
    return statistic_filter(df, sd_threshold, df.drop(['ID'], axis=1).std)


def group_miRNA(df):
    """
    group miRNA with the same CNV into 1 group
    :param df: the cancer / cnv DF
    :return: grouped df , dictionary of the groups
    """
    df = df.reset_index(drop=True)
    new_df = df.copy()
    counter = 0
    miRNA = df.iloc[:, 1:].to_numpy()
    idx = set(list(range(miRNA.shape[0])))
    names = df.iloc[:, 0].to_numpy().astype(str)
    groups_dict = {}
    grouped_df = pd.DataFrame(columns=df.columns)
    while (idx):
        print(len(idx))
        miR = miRNA[list(idx)[0]]
        similar = np.sum(miRNA - miR[np.newaxis, ...], 1)
        sim_idx = np.argwhere(similar == 0).flatten()
        idx = idx - set(sim_idx)
        if not len(sim_idx) - 1:
            continue
        groups_dict[counter] = names[sim_idx]
        new_df = new_df.drop(sim_idx)
        new_line = ['group %s' % counter] + miR.tolist()
        grouped_df.loc[counter] = new_line
        counter += 1
    grouped_df = pd.concat([new_df, grouped_df]).reset_index(drop=True)
    return grouped_df, groups_dict


# %%
def cluster_by_sd(df: pd.DataFrame, threshold:float):
    df =  sd_filtering(df, threshold)
    np_df = df.drop('ID', 1).to_numpy()
    kmean = KMeans()
    n_clusters = KElbowVisualizer(kmean,k=25)
    n_clusters.fit(np_df)
    n_clusters.show()
    kmean = KMeans(n_clusters.elbow_value_)
    pca = PCA(3)
    pca.fit(np_df)
    kmean.fit(np_df.T)
    plt.scatter(pca.components_[0], pca.components_[1], c=kmean.labels_)
    plt.xlabel("pca 1 (%s" %(pca.explained_variance_ratio_[0]*100)+'% explaind)')
    plt.ylabel("pca 2 (%s" %(pca.explained_variance_ratio_[1] *100)+'% explaind)')
    plt.title("t = %s" %threshold)
    plt.legend()
    plt.show()
    print(np_df.shape)
    fig = sns.clustermap(np_df)
    plt.show()





# %%

if __name__ == '__main__':
    df = gf.get_groups()
    cluster_by_sd(df, 0)



