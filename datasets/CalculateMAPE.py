"""
This script is used to calculate the Mean Absolute Percentage Error (MAPE) between the original data and the decompressed data.
It reads two CSV files with only one column, one containing the original data and the other containing the decompressed data,
then calculates and prints the MAPE.

"""

import numpy as np
import pandas as pd

# Calculate the Mean Absolute Percentage Error (MAPE)
def calculate_MAPE(y_true, y_decompressed):
    """
    Calculate the Mean Absolute Percentage Error (MAPE).

    Parameters:
        y_true (numpy.ndarray): The true values of the original data.
        y_decompressed (numpy.ndarray): The decompressed values.

    Returns:
        float: The MAPE value.
    """
    mask = y_true != 0
    return np.mean(np.abs((y_true[mask] - y_decompressed[mask]) / y_true[mask]))


if __name__=="__main__":

    """Modify the following paths to point to the paths of the two CSV files to be compared"""

    path1 = 'RSSI-2_T.csv'
    path2 = 'RSSI-2_T_decompressed.csv'

    # path1 = 'MIT-ECG_T.csv'
    # path2 = 'MIT-ECG_T_decompressed.csv'

    # path1 = 'S3-XYZ_T.csv'
    # path2 = 'S3-XYZ_T_decompressed.csv'

    # path1 = 'KPI-0_T.csv'
    # path2 = 'KPI-0_T_decompressed.csv'

    # path1 = 'Wafer_T.csv'
    # path2 = 'Wafer_T_decompressed.csv'

    # path1 = 'Wine_T.csv'
    # path2 = 'Wine_T_decompressed.csv'

    df1 = pd.read_csv(path1, header=None, usecols=[0])  
    df2 = pd.read_csv(path2, header=None, usecols=[0])

    original_values = df1.iloc[:, 0].to_numpy()
    decompressed_values = df2.iloc[:, 0].to_numpy()

    MAPE = calculate_MAPE(original_values, decompressed_values)
    MAPE_percentage = MAPE * 100
    print(f'MAPE: {MAPE_percentage:.2f}%')

