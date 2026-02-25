import pandas as pd
import numpy as np

df = pd.DataFrame({"organisms": ['', 'nan', np.nan, None, 'Panax']})
print("Original:")
print(df)

df['organisms'] = df['organisms'].replace(['', 'nan', None], 'Not Result').fillna('Not Result')
print("After replace:")
print(df)
