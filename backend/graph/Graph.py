import ctypes as c
import numpy as np
from random import randint

class Graph:

    def __init__(self, map_name="maps/test"):
        self.map_name = map_name
        self.lib = c.cdll.LoadLibrary('libelena.so')
        self.c_graph_handle = c.POINTER(c.c_char)

        # create a graph
        n_edges = (c.c_int * 1)()
        self.lib.createGraph.argtypes = [c.POINTER(c.c_char), c.POINTER(c.c_int)]
        self.lib.createGraph.restype=self.c_graph_handle
        self.c_graph = self.lib.createGraph(map_name.encode('utf-8'), n_edges)
        self.edges_buffer_size = n_edges[0]

        # getAllEdges
        self.lib.getAllEdges.argtypes = [
            self.c_graph_handle,
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='C_CONTIGUOUS'),
        ]
        self.lib.getAllEdges.res_type = c.c_int

        # getRandomNodeId
        self.lib.getRandomNodeId.argtypes = [self.c_graph_handle, c.c_int]
        self.lib.getRandomNodeId.restype = c.c_long

        # getNearestNode
        self.lib.getNearestNode.argtypes = [
            self.c_graph_handle,
            c.POINTER((c.c_double * 2))
        ]
        self.lib.getNearestNode.restype = c.c_long

        # aStarAlgorithm
        self.lib.aStarAlgorithm.argtypes = [
            self.c_graph_handle,
            c.c_long,
            c.c_long,
            c.POINTER((c.c_long * self.edges_buffer_size)), # can remove if not use (node ids)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
            c.POINTER(c.c_double),
            c.POINTER(c.c_double)
            # ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 2) * 100)
        ]
        self.lib.aStarAlgorithm.restype = c.c_int

        # dijkstraAlgorithm
        self.lib.dijkstraAlgorithm.argtypes = [
            self.c_graph_handle,
            c.c_long,
            c.c_long,
            c.POINTER((c.c_long * self.edges_buffer_size)), # can remove if not use (node ids)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
            c.POINTER(c.c_double),
            c.POINTER(c.c_double)
            # ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 2) * 100)
        ]
        self.lib.dijkstraAlgorithm.restype = c.c_int

        # findElevationBasedPath
        self.lib.findElevationBasedPath.argtypes = [
            self.c_graph_handle,
            c.c_long,
            c.c_long,
            c.c_double,
            c.c_bool,
            c.POINTER((c.c_long * self.edges_buffer_size)), # can remove if not use (node ids)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
            c.POINTER(c.c_double),
            c.POINTER(c.c_double)
            # ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 2) * 100)
        ]
        self.lib.findElevationBasedPath.restype = c.c_int

    def getAllEdges(self):
        edges_buffer = np.empty((self.edges_buffer_size, 2, 2))
        size = self.lib.getAllEdges(self.c_graph, edges_buffer)
        return edges_buffer[:size]

    def getRandomNodeId(self):
        seed = randint(0, 2147483647)
        return self.lib.getRandomNodeId(self.c_graph, seed)

    def getNearestNode(self, x, y):
        coords_buffer = (c.c_double * 2)(x, y)
        node_id = self.lib.getNearestNode(self.c_graph, coords_buffer)
        return node_id, (coords_buffer[0], coords_buffer[1])

    def aStarAlgorithm(self, start_id, end_id):
        path_buffer = (c.c_long * self.edges_buffer_size)()
        coords_buffer = np.empty((self.edges_buffer_size, 2))
        length, elevation = (c.c_double * 1)(), (c.c_double  * 1)()
        size = self.lib.aStarAlgorithm(
            self.c_graph,
            start_id,
            end_id,
            path_buffer,
            coords_buffer,
            length,
            elevation
        )
        return coords_buffer[:size], (length[0], elevation[0])

    def dijkstraAlgorithm(self, start_id, end_id):
        path_buffer = (c.c_long * self.edges_buffer_size)()
        coords_buffer = np.empty((self.edges_buffer_size, 2))
        length, elevation = (c.c_double * 1)(), (c.c_double  * 1)()
        size = self.lib.dijkstraAlgorithm(
            self.c_graph,
            start_id,
            end_id,
            path_buffer,
            coords_buffer,
            length,
            elevation
        )
        return coords_buffer[:size], (length[0], elevation[0])


    def findElevationBasedPath(self, start_id, end_id, max_weight_ratio=1.5, minimize=True):
        path_buffer = (c.c_long * self.edges_buffer_size)()
        coords_buffer = np.zeros((self.edges_buffer_size, 2))
        length, elevation = (c.c_double * 1)(), (c.c_double  * 1)()
        size = self.lib.findElevationBasedPath(
            self.c_graph,
            start_id,
            end_id,
            max_weight_ratio,
            minimize,
            path_buffer,
            coords_buffer,
            length,
            elevation
        )
        return coords_buffer[:size], (length[0], elevation[0])



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    map_name = 'maps/test'
    # map_name = 'maps/Helsinki'


    # This project
    g = Graph(map_name)
    # Find nearest node
    x, y = 26.9500, 60.5300
    node_id, (xx, yy) = g.getNearestNode(x, y)
    source, target = g.getRandomNodeId(), g.getRandomNodeId()

    # Get random node ids
    source, target = g.getRandomNodeId(), g.getRandomNodeId()


    a_star, (length, elevation) = g.aStarAlgorithm(source, target)
    print(f"length: {length:.2f}, elevation: {elevation:.0f}")

    dijkstra, (length, elevation) = g.dijkstraAlgorithm(source, target)
    print(f"length: {length:.2f}, elevation: {elevation:.0f}")

    elena, (length, elevation) = g.findElevationBasedPath(source, target)
    print(f"length: {length:.2f}, elevation: {elevation:.0f}")


    # Get all edges to plot map
    edges = g.getAllEdges()

    # Plot results
    for edge in edges:
        (x1, y1), (x2, y2) = edge
        plt.plot([x1, x2], [y1, y2], color='lightgray')

    for path, color in [(a_star, 'darkorange'), (dijkstra, 'dodgerblue'), (elena, 'forestgreen')]:
        for i in range(len(path) - 1):
            plt.plot([path[i][0], path[i + 1][0]], [path[i][1], path[i + 1][1]], color=color)

    plt.show()
