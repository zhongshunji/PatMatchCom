"""
This script partitions the dataset into Dictionary Data (D) and Input Data (T) in a 50:50 ratio.
"""
import pandas as pd
import numpy as np

# This function converts a df to 1 column and saves it as a csv file
def df_to_1column(df, csv_filename):

    flattened_data = df.values.flatten(order='F')

    df_flattened = pd.DataFrame({'Flattened': flattened_data})

    df_flattened.to_csv(csv_filename, index=False, header=False)


def dataset_split(df, filename1, filename2):
    
    half_num_rows = int(df.shape[0] / 2)  

    df1 = df[:half_num_rows]  
    df2 = df[half_num_rows:]  

    df_to_1column(df1, filename1)
    df_to_1column(df2, filename2)


if __name__ == "__main__":
    path10 = "./raw_datasets/RSSI-2.csv"
    df10 = pd.read_csv(path10, usecols=[4], header=None)  
    dataset_split(df10, 'RSSI-2_D.csv', 'RSSI-2_T.csv')

    import wfdb
    path13 = './raw_datasets/100'
    record = wfdb.rdrecord(path13, channel_names=['MLII'])  
    df13 = pd.DataFrame(record.p_signal, columns=['Signal'])
    dataset_split(df13, 'MIT-ECG_D.csv', 'MIT-ECG_T.csv')

    path16 = './raw_datasets/Samsung-Galaxy-S3.csv'
    df16 = pd.read_csv(path16, usecols=[3, 4, 5])
    dataset_split(df16, 'S3-XYZ_D.csv', 'S3-XYZ_T.csv')

    path1 = "./raw_datasets/0.csv"
    df1 = pd.read_csv(path1)
    dataset_split(df1, 'KPI-0_D.csv', 'KPI-0_T.csv')

    path23 = "./raw_datasets/Wafer_TRAIN.tsv"
    df23 = pd.read_csv(path23, header=None, usecols=range(1, 153), sep='\t')
    df23 = df23.T
    df_to_1column(df23, 'Wafer_D.csv')
    path24 = "./raw_datasets/Wafer_TEST.tsv"
    df24 = pd.read_csv(path24, header=None, usecols=range(1, 153), sep='\t')
    df24 = df24.T
    df_to_1column(df24, 'Wafer_T.csv')

    path21 = "./raw_datasets/Wine_TRAIN.tsv"
    df21 = pd.read_csv(path21, header=None, usecols=range(1, 235), sep='\t')
    df_to_1column(df21, 'Wine_D.csv')
    path22 = "./raw_datasets/Wine_TEST.tsv"
    df22 = pd.read_csv(path22, header=None, usecols=range(1, 235), sep='\t')
    df_to_1column(df22, 'Wine_T.csv')

    print("Dataset splitting completed")
    






