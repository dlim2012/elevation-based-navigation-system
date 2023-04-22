import numpy as np
import pandas as pd
from pyrosm import OSM
from pyrosm import get_data
from tqdm import tqdm
import requests

import math
import os

open_elevation_url = "http://192.168.1.20:1000/api/v1/lookup"

pd.set_option('display.max_columns', None)

# filename = '/mnt/ssd/osm/us-secondary'
# filename = '/mnt/ssd/osm/us-tertiary'

# osm = OSM(f"{filename}.osm.pbf")
osm = OSM(get_data('test_pbf')); filename = 'test'
# osm = OSM(get_data("Helsinki")); filename = "Helsinki"
save_dir = "maps"

nodes, edges = osm.get_network(network_type="driving", nodes=True)

print(len(nodes), len(edges))

nodes = nodes[['id', 'lon', 'lat']]
edges = edges[['u', 'v', 'oneway']]

nodes = nodes.to_dict()
nodes['elevation'] = {}

batch_size = 1000
indices = list(nodes['id'].keys())
# for i, index in enumerate(nodes['id'].keys()):
for i in tqdm(range(0, len(indices), batch_size)):
    response = requests.post(
        open_elevation_url,
        json={
            'locations':[
                {
                    "latitude": nodes['lat'][index],
                    "longitude": nodes['lon'][index]
                } for index in indices[i: i+batch_size]
            ]
        }
    )
    assert response.status_code == 200
    for j, index in enumerate(indices[i: i+batch_size]):
        nodes['elevation'][index] = response.json()['results'][j]['elevation']


coords = {
    nodes['id'][index]: (nodes['lat'][index], nodes['lon'][index])
        for index in nodes['id'].keys()
}
elevation = {
    nodes['id'][index]: nodes['elevation'][index]
        for index in nodes['id'].keys()
}


nodes = pd.DataFrame.from_dict(nodes)
print("nodes: ", len(nodes))
nodes.to_csv(f"{os.path.join(save_dir, filename)}.nodes.csv", float_format="%.7f")

def distance(u, v, coords):
    (lat1, lon1) = coords[u]
    (lat2, lon2) = coords[v]
    lat1, lon1, lat2, lon2 = lat1 * math.pi / 180, lon1 * math.pi / 180, lat2 * math.pi / 180, lon2 * math.pi / 180
    return math.acos(math.sin(lat1)*math.sin(lat2)+math.cos(lat1)*math.cos(lat2)*math.cos(lon2-lon1))*6371*1000
def distance2(u, v, coords):
    (lat1, lon1) = coords[u]
    (lat2, lon2) = coords[v]
    return math.sqrt((lat1-lat2)**2+(lon1-lon2)**2)

edges['oneway'] = np.where(edges['oneway'] == "yes", True, False)
edges = edges.to_dict()
edges['length'] = {}
edges['elevation'] = {}
for i in range(len(edges['u'])):
    u, v = edges['u'][i], edges['v'][i]
    edges['length'][i] = distance(u, v, coords)
    edges['elevation'][i] = abs(elevation[u] - elevation[v])
    edges['oneway'][i] = 1 if edges['oneway'][i] else 0

edges = pd.DataFrame.from_dict(edges)
print("edges:", len(edges))
edges.to_csv(f"{os.path.join(save_dir, filename)}.edges.csv", float_format="%.5f")
