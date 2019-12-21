from Queue import *
import networkx as nx

# class EdgeNode:
#     def __init__(self, node, weight = 0):
#         self.node = node
#         self.weight = weight


# class Graph:

#     def __init__(self, nVertices):
#         self.numberVertices = nVertices;
#         self.adjacency_list = {}
#         self.vertices = []

#         for x in range(1, nVertices + 1):
#             self.adjacency_list[x] = []
#             self.vertices.append(x)

#         self.dist = {}
#         for x in range(1, nVertices + 1):
#             self.dist[x] = float('inf')

#         self.pred = {}
#         for x in range(1, nVertices):
#             self.pred[x] = None


# def add_edge(graph, x, y, weight):
#     graph.adjacency_list[x].append(EdgeNode(y, weight))
#     graph.adjacency_list[y].append(EdgeNode(x, weight))


parent = dict()
rank = dict()

def make_set(v):
    parent[v] = v
    rank[v] = 0

def find(v):
    if parent[v] == v:
        return parent[v]
    else:
        return find(parent[v])

def union(x, y):
    xroot = find(x)
    yroot = find(y)

    if xroot != yroot:
        if rank[xroot] > rank[yroot]:
            parent[yroot] = xroot
        else:
            parent[xroot] = yroot
        if rank[xroot] == rank[yroot]:
            rank[yroot] = rank[yroot] + 1



def kruskal(graph):

    mst = []
#    for vertice in graph.vertices:
#        make_set(vertice)
    for vertices in nx.nodes(graph):
        make_set(vertices)
    # edges = []
    # for key in graph.adjacency_list:
    #     for edge in graph.adjacency_list[key]:
    #         edges.append((edge.weight, key, edge.node))
    # alledge = list(graph.edges())
    # alledge = nx.get_edge_attributes(graph, 'weight')
    # for key in alledge:
    #     print key
    #     edges.append((graph.get_edge_data(key[0], key[1]), key[0], key[1]))
    #     # edges.append((alledge[key], key[0], key[1]))
    edges = [ (edata['weight'], u, v) for u,v,edata in graph.edges(data=True) 
               if 'weight' in edata]
    # print len(edges)
    edges.sort()
    # print edges

    for edge in edges:
        weight, v1, v2 = edge
        if find(v1) != find(v2):
            mst.append(edge)
            union(v1, v2)

    return mst

def neighbors(list, node):

    neighbors = []

    for tuple in list:
        if node == tuple[1] or node == tuple[2]:
            neighbors.append(tuple[1])
            neighbors.append(tuple[2])
    return neighbors



#Breadth first search algorithm
def bfs(graph, start, goal):
    open_set = Queue()
    open_set.put(start)
    came_from = {}
    came_from[start] = None

    while not open_set.empty():
        current = open_set.get()
        if current == goal:
            break
        for next in neighbors(graph, current):
            if next not in came_from:
                open_set.put(next)
                came_from[next] = current


    current = goal
    path = []
    while current != start:
        path.append(current)
        current = came_from[current]


    path.append(start)
    path.reverse()

    return path


def get_highest(mst, path):

    edges = list()
    for i in range(len(path)):
        if i < len(path) - 1:
            edges.append((path[i], path[i+1]))

    cost = 0
    temp = 0

    for edge in edges:
        for tuple in mst:
            if edge[0] == tuple[1] and edge[1] == tuple[2]:
                temp = tuple[0]
                if temp > cost:
                    cost = temp

    return cost
