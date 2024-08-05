#Youssef Ehab Maher
#20210465
class Graph:
    def __init__(self):
        self.graph = {}

    def add_edge(self, node1, node2, cost):
        if node1 not in self.graph:
            self.graph[node1] = []
        self.graph[node1].append((node2, cost))

    def print_solution(self, dist,source):
        for vertex, distance in dist.items():
            if(vertex!=source):
                print('distance from',source,'-->',vertex,'is',distance)

    def shortest_dist(self, source):
        dist = {vertex: float("inf") for vertex in self.graph}
        dist[source] = 0
        for _ in range(len(self.graph) - 1):
            for node1 in self.graph:
                for node2, cost in self.graph[node1]:
                    if dist[node1] + cost < dist[node2]:
                        dist[node2] = dist[node1] + cost
        self.print_solution(dist,source)

#add edge as first node--> second node, cost
#nodes must be strings

graph1 = Graph()
graph1.add_edge('1', '2', 7)
graph1.add_edge('1', '3', 8)
graph1.add_edge('1', '4', 5)
graph1.add_edge('2', '5', 12)
graph1.add_edge('3', '5', 8)
graph1.add_edge('3', '6', 9)
graph1.add_edge('4', '5', 9)
graph1.add_edge('4', '6', 13)
graph1.add_edge('5', '7', 9)
graph1.add_edge('6', '7', 6)
graph1.add_edge('7','7',0)

graph1.shortest_dist('1')
