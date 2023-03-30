import ctypes as c
import numpy as np
from random import randint

class Graph:

    def __init__(self, version_name="maps/test.maxGroup"):
        self.version_name = version_name
        self.lib = c.cdll.LoadLibrary('libelena.so')
        c_graph_handle = c.POINTER(c.c_char)
        self.path_buffer_size = 100000

        # createGraph
        self.lib.createGraph.restype=c_graph_handle

        # shortestPath
        self.lib.shortestPath.argtypes = [
            c_graph_handle,
            c.c_long,
            c.c_long,
            c.POINTER((c.c_long * self.path_buffer_size)), # can remove if not use (node ids)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
            c.POINTER(c.c_double)
            # ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 2) * 100)
        ]
        self.lib.shortestPath.restype = c.c_int

        # getRandomNodeId
        self.lib.getRandomNodeId.argtypes = [c_graph_handle, c.c_int]
        self.lib.getRandomNodeId.restype = c.c_long

        # getNearestNode
        self.lib.getNearestNode.argtypes = [c_graph_handle, c.POINTER((c.c_double * 2))]
        self.lib.getNearestNode.restype = c.c_long

        # findRandomPath
        self.lib.findRandomPath.argtypes = [
            c_graph_handle,
            c.c_long,
            c.c_long,
            c.POINTER((c.c_long * self.path_buffer_size)), # can remove if not use (node ids)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
            c.POINTER(c.c_double * 1),
            c.c_int
            # ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 2) * 100)
        ]
        self.lib.findRandomPath.restype = c.c_int


        # create a graph
        self.c_graph = self.lib.createGraph(version_name.encode('utf-8'))


        #
        # self.path_buffer = (c.c_long * self.buffer_size)()
        # self.coords_buffer = np.zeros((self.buffer_size, 2))

    def RandomNodeId(self):
        seed = randint(0, 2147483647)
        return self.lib.getRandomNodeId(self.c_graph, seed)

    def getNearestNode(self, x, y):
        coords_buffer = (c.c_double * 2)(x, y)
        node_id = self.lib.getNearestNode(self.c_graph, coords_buffer)
        return node_id, (coords_buffer[0], coords_buffer[1])

    def shortestPath(self, start_id, end_id):
        path_buffer = (c.c_long * self.path_buffer_size)()
        coords_buffer = np.zeros((self.path_buffer_size, 2))
        length = (c.c_double * 1)()
        size = self.lib.shortestPath(
            self.c_graph,
            start_id,
            end_id,
            path_buffer,
            coords_buffer,
            length
        )
        return length[0], coords_buffer[:size]

    def randomPath(self, start_id, end_id, length=float('inf')):
        seed = randint(0, 2147483647)
        path_buffer = (c.c_long * self.path_buffer_size)()
        coords_buffer = np.zeros((self.path_buffer_size, 2))
        length = (c.c_double * 1)(length)
        size = self.lib.findRandomPath(
            self.c_graph,
            start_id,
            end_id,
            path_buffer,
            coords_buffer,
            length,
            seed
        )
        return length[0], coords_buffer[:size]

g = Graph()
g.shortestPath(1364765716, 3350088311)