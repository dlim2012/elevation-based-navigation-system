import numpy as np
import pandas as pd
from pyrosm import OSM
from pyrosm import get_data

pd.set_option('display.max_columns', None)

# filename = '/mnt/ssd/osm/us-secondary'
# filename = '/mnt/ssd/osm/us-tertiary'

# osm = OSM(f"{filename}.osm.pbf")
osm = OSM(get_data('test_pbf')); filename = 'test'
nodes, edges = osm.get_network(network_type="driving", nodes=True)

print(len(nodes), len(edges))

nodes = nodes[['id', 'lon', 'lat']]

edges = edges[['u', 'v', 'length', 'oneway']]
edges['oneway'] = np.where(edges['oneway'] == "yes", True, False)

nodes.to_csv(f"{filename}.nodes.csv", float_format="%.7f")
edges.to_csv(f"{filename}.edges.csv", float_format="%.3f")