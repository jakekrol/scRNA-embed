#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys

# Get the .npz file from command line arguments
print('test')
npz_file = sys.argv[1]

# Load the .npz file
data = np.load(npz_file, allow_pickle=True)

# Extract the data, index, and columns
matrix_data = data['data']
index = data['index']
columns = data['columns']

# Create the DataFrame
restored_df = pd.DataFrame(matrix_data, index=index, columns=columns)

# Print the head and shape of the DataFrame
print("Head of the DataFrame:")
print(restored_df.head())

print("\nShape of the DataFrame:")
print(restored_df.shape)