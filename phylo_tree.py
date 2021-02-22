import graph as gx
import copy
import pandas as pd
import sys 

class Phylo_Tree():
    def __init__(self, G, distance_matrix, species_dict, file_name, algo_type):
        distance_matrix_copy = copy.deepcopy(distance_matrix)
        G_copy = copy.deepcopy(G)
        species_dict_copy = copy.deepcopy(species_dict)
        self.distance_matrix = distance_matrix_copy
        self.species_dict = species_dict_copy
        self.G = G_copy 
        self.algo_type = algo_type
        self.file_name = file_name
        self.newick_form = ""
        if self.algo_type == "UPGMA":
           self.newick_form = self.UPGMA(self.G, self.distance_matrix, self.species_dict)
           #self.newick_form = self.G.to_newick_UPGMA()
        elif self.algo_type == "NEIGHBOR":
            self.neighbor_joining(self.G, self.distance_matrix, self.species_dict)
            self.newick_form = self.G.to_newick_NJ()
        else: 
            print("Algorithm type unkown. Enter a valid algorithm to compute tree.")
            sys.exit(1) 
        print(self.newick_form)
        self.write_to_file(self.file_name, self.newick_form)

    def write_to_file(self, file_name, newick_form): 
        outfile = open(file_name, "w+")
        print(self.newick_form, file = outfile)
        outfile.close()

    def neighbor_joining(self, G, distance_matrix, dictionary):
        inverse_dict = {v: k for k, v in dictionary.items()}
        while(len(distance_matrix) > 1): 
            minimum, c_i, c_j, u_c_1, u_c_2 = self.return_neighbor_distance(dictionary, distance_matrix)
            #print("c_i: ", inverse_dict[c_i], "c_j: ", inverse_dict[c_j])
            new_node = self.add_new_node(G, c_i, c_j, inverse_dict, dictionary, minimum)
            #print(new_node)
            self.update_edge_weights(G, minimum, new_node, c_i, c_j, u_c_1, u_c_2, inverse_dict, distance_matrix)
            self.update_distance_matrix(G, new_node, distance_matrix, inverse_dict, c_i, c_j, dictionary)
            self.delete_rows(dictionary, inverse_dict, c_i, c_j, distance_matrix, new_node)

    def update_edge_weights(self, G, minimum, new_node, c_i, c_j, u_c_1, u_c_2, inverse_dict, distance_matrix):
        species1 = inverse_dict[c_i]
        species2 = inverse_dict[c_j]
        temp = sorted([c_i, c_j])
        distance = distance_matrix[temp[0]][temp[1]]
        weight1 = (distance / 2) + ((u_c_1 - u_c_2) / 2)
        G.update_edge_weight((new_node, species1), weight1)
        weight2 = (distance / 2) + ((u_c_2 - u_c_1) / 2)
        G.update_edge_weight((new_node, species2), weight2)

    def compute_new_distances(self, G, new_node, distance_matrix, inverse_dict, c_i, c_j, dictionary): 
        for row, idx in enumerate(range(len(distance_matrix))): 
            if idx != c_i and idx != c_j: 
                new_dist = 0
                temp = sorted([idx, c_i])
                new_dist += distance_matrix[temp[0]][temp[1]]
                temp = sorted([idx, c_j])
                new_dist += distance_matrix[temp[0]][temp[1]]
                new_dist = new_dist / 2
                temp = sorted([idx, c_i])
                distance_matrix[temp[0]][temp[1]] = new_dist

    def return_neighbor_distance(self, dictionary, distance_matrix): 
        index_1 = None
        index_2 = None
        U_C_1 = None
        U_C_2 = None
        separation_dist = None
        minimum = 1000
        for i, row in enumerate(range(len(distance_matrix))): 
            for j, col in enumerate(range(len(distance_matrix))):
                if i > j: 
                    if distance_matrix[i][j] < minimum:
                        minimum = distance_matrix[i][j]
                        u_c_1 = self.return_sum_distance(i, distance_matrix) 
                        u_c_2 = self.return_sum_distance(j, distance_matrix) 
                        separation_dist = distance_matrix[row][col] - u_c_1 - u_c_2
                        index_1 = i
                        index_2 = j
                        U_C_1 = u_c_1
                        U_C_2 = u_c_2
        return(separation_dist, index_1, index_2, U_C_1, U_C_2)

    def return_sum_distance(self, i, distance_matrix): 
        sum = 0
        for idx, element in enumerate(range(len(distance_matrix))): 
            sum += distance_matrix[i][idx]
        return sum / (len(distance_matrix))


    def UPGMA(self, G, distance_matrix, dictionary):
        def newickify(node_to_children, root_node): 
            visited_nodes = set()

            def newick_render_node(name, distance): 
                assert name not in visited_nodes, "Error: circular"
                if name not in node_to_children: 
                    #leaves
                    return F'{name}'
                else: 
                    visited_nodes.add(name)
                    children = node_to_children[name]
                    children_strings = [newick_render_node(child, children[child]) for child in children.keys()]
                    children_strings = ",".join(children_strings)
                    return F'({children_strings}){name}'

        #This tree_dict stuff is, admittedly, a pretty terrible work around for not having a binary tree implementation. You live and you learn.
        tree_dict = {}  
        leaves = G.nodes()
        """
        for leaf in leaves: 
            tree_dict[leaf] = {None, None} """
        n_string = ""
        print("leaves: ", leaves)
        newick_str = "("
        node_set = []
        counter = 0
        inverse_dict = {v: k for k, v in dictionary.items()}
        while(len(distance_matrix) > 1):
            minimum, c_i, c_j = self.return_upgma_distance(dictionary, distance_matrix)
            newick_info = [inverse_dict[c_i], inverse_dict[c_j]]
            #self.print_results(distance_matrix, dictionary)
            new_node = self.add_new_node(G, c_i, c_j, inverse_dict, dictionary, minimum)
            #print("species1: ", inverse_dict[c_i], "ci: ", c_i, " species2: ", inverse_dict[c_j],", cj: ", c_j, ", min: ", minimum, "new_node: ", new_node)
            self.update_node_height_UPGMA(G, minimum, new_node)
            self.update_distance_matrix(G, new_node, distance_matrix, inverse_dict, c_i, c_j, dictionary)
            self.delete_rows(dictionary, inverse_dict, c_i, c_j, distance_matrix, new_node)
            #self.to_newick(newick_str, newick_info, new_node, dictionary)
            newick_str += "(" + str(newick_info[0]) + ", " + str(newick_info[1]) + ")" 
            if new_node not in node_set:
                tree_dict[new_node] = {str(newick_info[0]): 0, str(newick_info[1]): 0} 
                newick_str += new_node + ")"
                node_set.append(new_node)
            counter += 1
        x = counter * "("
        x += newick_str
        #G.print_graph()
        #print("tree dictionary: ")
        #print(tree_dict)

        def newickify(node_to_children, root_node):
            visited_nodes = set()
            def newick_render_node(name, distance): 
                assert name not in visited_nodes, "Error: circular"
                if name not in node_to_children: 
                    #leaves
                    return F'{name}'
                else: 
                    visited_nodes.add(name)
                    children = node_to_children[name]
                    children_strings = [newick_render_node(child, children[child]) for child in children.keys()]
                    children_strings = ",".join(children_strings)
                    return F'({children_strings}){name}'
            newick_s = newick_render_node(root_node, 0) + ';'
            assert visited_nodes == set(node_to_children.keys()), "Error"
            return newick_s

        head = G.return_head()
        print("*" * 30)
        print(newickify(tree_dict, root_node = head))
        print("*" * 30)

        return x

    def to_newick(self, newick_str, newick_info, new_node, dictionary): 
        newick_str = "(" + newick_str
        newick_str += str(newick_info[0]) + ", " + str(newick_info[1]) + ", " + new_node + ")" 
        print(newick_str)

    
    def return_upgma_distance(self, dictionary, distance_matrix):
        index_1 = None
        index_2 = None
        minimum = 1000
        for i, row in enumerate(range(len(distance_matrix))): 
            for j, column in enumerate(range(len(distance_matrix))): 
                if j > i:
                    if distance_matrix[i][j] < minimum:
                        minimum = distance_matrix[i][j]
                        index_1 = i
                        index_2 = j
        return (minimum, index_1, index_2)

    def update_node_height_UPGMA(self, G, minimum, new_node): 
        minimum = minimum / 2
        G.change_height(new_node, minimum)

    def add_new_node(self, G, c_i, c_j, inverse_dict, dictionary, minimum):
        temp = sorted([c_i, c_j])
        #Look up the species name of the index values we are dealing with and combine into a string
        new_node = inverse_dict[temp[0]] + "-" + inverse_dict[temp[1]]
        #Add a new node with the combined string name (example: chimp_human) to the graph object
        G.add_node(new_node)
        #Connect the new node to the two child nodes
        G.add_edges([(new_node, inverse_dict[c_i], 0), (new_node, inverse_dict[c_j], 0)])
        #Find how many actual sequences (non internal nodes) are in the two child nodes that we just combined
        #Add together to get the count of sequences in both clusters.
        sequences = G.get_sequence_count(inverse_dict[c_i]) + G.get_sequence_count(inverse_dict[c_j])
        #Add the sequence number to the new node
        G.increment_sequence_count(new_node, sequences) 
        #Return the string name of the new node
        return new_node

    def update_distance_matrix(self, G, new_node, distance_matrix, inverse_dict, c_i, c_j, dictionary):
        for row, idx in enumerate(range(len(distance_matrix))):
            if idx != c_i and idx != c_j: 
                c_i_sequence_count = G.get_sequence_count(inverse_dict[c_i])
                c_j_sequence_count = G.get_sequence_count(inverse_dict[c_j])
                temp = sorted([idx, c_i])
                distance1 = distance_matrix[temp[0]][temp[1]]
                temp = sorted([idx, c_j])
                distance2 = distance_matrix[temp[0]][temp[1]]
                numerator = distance1 * c_i_sequence_count + distance2 * c_j_sequence_count 
                denominator = c_i_sequence_count + c_j_sequence_count
                dist = numerator / denominator
                #Overwrite c_i with new value
                temp = sorted([idx, c_i])
                distance_matrix[temp[0]][temp[1]] = dist

    def delete_rows(self, dictionary, inverse_dict, c_i, c_j, distance_matrix, new_node):
        #Delete column c_j from distance matrix (c_i was overwritten with the new values)
        [j.pop(c_j) for j in distance_matrix]
        #Delete row c_j from distance matrix
        distance_matrix.pop(c_j)
        #Update the dict value for the new node
        dictionary[new_node] = c_i
        del dictionary[inverse_dict[c_i]]
        del dictionary[inverse_dict[c_j]]
        del inverse_dict[c_i]
        del inverse_dict[c_j]
        #Delete the value of c_j in the dictionary and update the values of all other dictionary items
        sorted_dict = {k: v for k, v in sorted(dictionary.items(), key = lambda item: item[1])}
        dict_list = list(sorted_dict.keys())
        for idx, element in enumerate(dict_list):
            dictionary[element] = idx
            inverse_dict[idx] = element

    def print_results(self, results, dictionary):
        print(dictionary)
        pd.set_option("display.max_rows", None)
        pd.set_option("display.max_columns", None)
        df = pd.DataFrame(results, columns= [i for i in range(len(results))])
        print(df)

