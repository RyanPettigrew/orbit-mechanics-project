import json
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np

with open('russia_data.json', 'r') as file:
    data = json.load(file)

# need to do for china too
with open('china_data.json', 'r') as file:
    data2 = json.load(file)

extracted_data = []
for item in data:
    extracted_data.append({
        'object_name': item['OBJECT_NAME'],
        'object_id': item['OBJECT_ID'],
        'epoch': item['EPOCH'],
        'mean_motion': item['MEAN_MOTION'],
        'eccentricity': item['ECCENTRICITY'],
        'inclination': item['INCLINATION'],
        'ra_of_asc_node': item['RA_OF_ASC_NODE'],
        'arg_of_pericenter': item['ARG_OF_PERICENTER'],
        'mean_anomaly': item['MEAN_ANOMALY'],
        'norad_cat_id': item['NORAD_CAT_ID'],
    })

df = pd.DataFrame(extracted_data)
print(df.head())

relevant_data = df[['eccentricity', 'inclination', 'mean_motion']]
normalized_data = (relevant_data - relevant_data.mean()) / relevant_data.std()
kmeans = KMeans(n_clusters=4, random_state=0).fit(normalized_data) # find clusters
df['cluster'] = kmeans.labels_
selected_debris = df.groupby('cluster').apply(lambda x: x.sample(1)).reset_index(drop=True)
print(selected_debris)
