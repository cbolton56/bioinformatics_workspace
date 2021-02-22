""" A Graph class implemented in python. 
A graph represented in Python using dictionary objects.
Each instance of the graph has a member variable __graph_dict that
is a dictionary object. The keys of this object are nodes of the graph. 
Each key value pair has a value of list. The list can be filled with other nodes.
Examples:
        Basic layout:  

        __graph_dict = {"key": [["edge", "edge"], float_height_of_node, int_sequence_count]}

        Empty graph: 

        __graph_dict = {}

        Graph object with multiple nodes that may have multiple edges: 

        __graph_dict = {"suzy": ["bob", "carol"], "bob":["suzy"], "carol": ["suzy"]}      """
class Graph(): 
    # Constructor
    def __init__(self, graph_dict = None):
        if graph_dict == None: 
            graph_dict = {}
        self.__graph_dict = graph_dict

    #accessor function returns the nodes in a graph
    def nodes(self):
        #nodes are represented by keys in the dictionary object
        return list(self.__graph_dict.keys())

    def edges(self):
        return self.get_edges()

    def add_node(self, node): 
        #If the node isn't in the graph, this method adds the node.
        #Adds a key to the graph_dict object given by node paramter. 
        #Empty by default
        if node not in self.__graph_dict: 
            self.__graph_dict[node] = [[], 0, 0]
    
    #Takes a list of nodes and adds to the graph
    def add_nodes(self, nodes): 
        for node in nodes:
            self.add_node(node)
    
    def add_edge(self, edge): 
        #Assumes edge is of type tuple, i.e.: (node, node, edge-weight)
        if edge[0] not in self.__graph_dict: 
            self.__graph_dict[edge[0]] = [[], 0, 0] 
        if edge[1] not in self.__graph_dict: 
            self.__graph_dict[edge[1]] = [[], 0, 0] 
        self.__graph_dict[edge[0]][0].append([edge[1], edge[2]])
        #self.__graph_dict[edge[1]][0].append([edge[0], edge[2]])

    def add_edges(self, edges): 
        #method for adding multiple edges at one time
        #takes a list of tuples
        for edge in edges:
            self.add_edge(edge)
    
    def disconnected_graph(self): 
        #Method returns true or false value. 
        #True means that the graph object is a connected graph and all the vertices are connected to at least one other vertex. 
        #False means there is at least one node that is still disconnected from the rest of the graph
        for node in self.__graph_dict.keys(): 
            if len(self.__graph_dict[node][0]) == 0: 
                return True
        return False

    def connected_graph(self): 
        connected_vertices = 0
        for node in self.__graph_dict.keys(): 
            if len(self.__graph_dict[node][0]) > 0: 
                connected_vertices += 1
        if connected_vertices == len(self.__graph_dict.keys()): 
            return True 
        return False 
        
    def change_height(self, node, distance): 
        self.__graph_dict[node][1] += distance
    
    def return_height(self, node): 
        return self.__graph_dict[node][1]
    
    def return_edge_weight(self, nodes):
        species1 = nodes[0]
        species2 = nodes[1]
        edge_list = self.__graph_dict[species1][0]
        idx = [x[0] for x in edge_list].index(species2)
        return self.__graph_dict[species1][0][idx][1]

    def update_edge_weight(self, nodes, weight): 
        #where "nodes" is a tuple object of the weight we want updated ("human-chimp", "human")
        #Find the string values of the species, unpacking the nodes variable
        species1 = nodes[0]
        species2 = nodes[1]
        #print("species1: ", species1, "species2: ", species2)
        #Find the list of edges for the first species
        edge_list_1 = self.__graph_dict[species1][0] 
        #print("edge list 1: ", edge_list_1) 
        #Find the index value of the edge we want to add weight to (second species)
        idx = [x[0] for x in edge_list_1].index(species2)
        #Add the weight at that index
        self.__graph_dict[species1][0][idx][1] += weight 
        #Repeat for the other node we are adding weight to.
        edge_list_2 = self.__graph_dict[species2][0]
        #print("edge list 2: ", edge_list_2)
        idx = [x[0] for x in edge_list_2].index(species1)
        self.__graph_dict[species2][0][idx][1] += weight

    def return_smaller_subtree(self, node1, node2): 
        node1_sequence = self.get_sequence_count(node1)
        node2_sequence = self.get_sequence_count(node2)
        if node1_sequence > node2_sequence: 
            return node2
        else: 
            return node1

    def increment_sequence_count(self, node, sequences): 
        self.__graph_dict[node][2] += sequences

    def get_sequence_count(self, node): 
        return self.__graph_dict[node][2]

    def get_edges(self): 
        #Return all edge pairs as a list of tuples
        edges = []
        for node in self.__graph_dict: 
            for neighbor in self.__graph_dict[node][0]: 
                if (neighbor, node) not in edges: 
                    edges.append((node, neighbor))
        return edges
    
    def return_head(self): 
        species = ""
        sequence_count = 0
        for key in self.__graph_dict.keys(): 
            if self.__graph_dict[key][2] > sequence_count: 
                species = key
                sequence_count = self.__graph_dict[key][2]
        return species

    #Returns the direct neighbors of the node being passed in (one level out)
    def get_neighbors(self, node):
        neighbors = self.__graph_dict[node][0]
        return neighbors
    
    #Calculates degree of the node passed in 
    def get_degree(self, node): 
        degree = len(self.__graph_dict[node][0])
        return degree

    def return_node_count(self): 
        return len(self.__graph_dict.keys())

    def to_newick_UPGMA(self): 
        next_node = self.return_head()
        counter = 0
        newick = ""
        while(1): 
            children = self.get_neighbors(next_node)
            if len(children) > 1: 
                leaf_node = self.return_smaller_subtree(children[0][0], children[1][0])
                internal = [x[0] for x in children if x[0] != leaf_node]
                internal_node = internal[0]
                newick += "(" +str(next_node) + ":" + str(self.return_height(next_node)) + ")" + str(internal_node) + ":" + str(self.return_height(internal_node)) +  ", " + str(leaf_node) + ":" + str(self.return_height(internal_node))  +  ")" 
                next_node = internal_node
                counter += 1
            else: 
                break
        str1 = "(" * counter
        newick = str1 + newick
        return newick 

    def to_newick_NJ(self):
        next_node = self.return_head()
        counter = 0
        newick = ""
        while(1): 
            children = self.get_neighbors(next_node)
            if len(children) > 1: 
                leaf_node = self.return_smaller_subtree(children[0][0], children[1][0])
                internal = [x[0] for x in children if x[0] != leaf_node]
                internal_node = internal[0]
                newick = str(internal_node) + ", " + str(leaf_node) + ")" + newick
                next_node = internal_node
                counter += 1
            else: 
                break 
        str1 = "(" * counter 
        newick = str1 + newick
        return newick  
    

    def print_graph(self): 
        print(self.__graph_dict)

    #Converts a graph object to a JSON file. returns object to be passed into the JSON file for this application
    def to_json(self): 
        #Data is the main dictionary object to be converted to a JSON file
        data = {"nodes": [], "edges" : []}

        #Creates lists of all nodes and edges (edges are stored in tuples)
        nodes = list(self.nodes())
        edges = list(self.edges())

        #Node dictionary object to be added to the data object: 
        # {"id": the id of the node, "label", the name of the node, "x": the x coordinate
        # of the node on a plane, "y": the y coordinate of the node on a plane, "size": how
        # large the node is when displayed}
        #These individual node dictionaries will be appended to the node list in the data dict.
        for n in nodes:
            d = {"id": n, "label": n, "x": 0, "y": 0, "size": 3}
            data["nodes"].append(d)
        
        #Edge dictionary object to be added to the data object: 
        # {"id": the id of the node, "source": node the edge originates at, "target": end node of edge}
        #Note that the graph is not directed so source and target nodes are arbitrary.
        for idx, e in enumerate(edges):
            temp = {"id": idx, "source": e[0], "target": e[1]}
            data["edges"].append(temp)
        return data 

    #Overloaded __str__() method returns string representation of object
    def __str__(self):
        label = "nodes: "
        for k in self.__graph_dict:
            label += str(k) + " , "
        label += "\nedges: "
        for edge in self.get_edges():
            label += str(edge) + " "
        return label

    def return_graph(self): 
        return self.__graph_dict
