import pickle
import sys
import os
from yellowbrick.cluster import KElbowVisualizer
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd


MAX_RANGE = 4
MIN_RANGE = 3
ALL = 100
PRECENTAGE = 2
ROW = 1
DATA_PATH = 1


def check_means(data, range_min, range_max, metric_sys):

    model = KMeans()
    visualizer = KElbowVisualizer(model, k=(range_min, range_max), metric=metric_sys)
    visualizer.fit(data)
    visualizer.show()
    return visualizer.elbow_value_




def load_file(file_path):
    if file_path.endswith(".xlsx") or file_path.endswith(".xls") \
            or file_path.endswith(".csv"):
        data = pd.read_excel(file_path)
        data = data.to_numpy()
    else:
        data = pickle.load(open(file_path, "rb"))
        data = data.to_numpy()
    return data


def main():
    """
    arg[1] - data path
    arg[2] precentage of data to evaluate // 0% = 100%
    arg[3] min range
    arg[4] max range
    :return:
    """
    if not len(sys.argv - 1):
        print( "file_path precentage_of_sample")
    data = load_file(sys.argv[DATA_PATH])
    samples_num = data.shape[ROW] * float(sys.argv[PRECENTAGE])
    samples_num  = samples_num if samples_num else data.shape[ROW]
    samples_idx = np.random.choice(range(data.shape[ROW]), int(round(samples_num)), replace=False)
    samples = data[:][samples_idx]
    check_means(samples, int(sys.argv[MIN_RANGE]), int(sys.argv[MAX_RANGE]), 'distortion')
    check_means(samples, int(sys.argv[MIN_RANGE]), int(sys.argv[MAX_RANGE]), 'calinski_harabasz')






# if __name__ == '__main__':
#     main()


