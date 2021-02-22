import os
import sys
import pandas as pd
import graph as gx
import phylo_tree as px

FILE_PATH = "/home/turing/mhou/csci652/data/sars-cov2/UCSC/wuhCor1.strainName119way.maf"
ASCII_LENGTH = 128

def main():
    #Opens the main file that contains the alignment blocks
    infile = open(FILE_PATH)
    lines = infile.readlines()
    #Temp stores all of the lines that start with 's' to isolate the sequences
    temp = []
    species = []
    dictionary = {}
    #Only store lines in temp that contain sequence information
    counter = 0 
    for line in lines:
        if line[0] == "s": 
            temp.append(line)
            counter += 1
    """The create_dict function makes a dictionary to store the species names and the index value of each spices in the distance matrix.
    example: {"wigeon": 0} means the distances for wigeon are stored at index 0 in the distance matrix. """
    create_dict(temp, dictionary)
    infile.close()
    #Open the file that stores the position information for each gene in the sequence.
    genefile = open("/home/turing/mhou/csci652/data/sars-cov2/UCSC/sarsCov2structure.txt")
    lines = genefile.readlines()
    """Gene_dict is a dictionary for storing the name of the gene and the position information for each gene. 
    For example: {"S": (1232, 4254)}  The tuple stores start and end information. """
    gene_dict = {}
    for line in lines: 
        #print(line)
        arr = line.split()
        #print(arr)
        if len(arr):
            gene_dict[arr[0]] = (int(arr[1]), int(arr[2]))
    genefile.close()
    counter = 0
    "Compute the gene block, which is an array of sequences with columns remove that correspond with gaps in the reference sequence"
    sequences = return_gene_block(temp)
    "Each gene is stored in a dictionary that holds beginning and ending values as a tuple, ie: {S: (52, 322)}"
    for key in gene_dict.keys():
        #Declare a new distance matrix for each gene
        distance_matrix = [[0 for x in range(len(temp))] for x in range(len(temp))]
        #Assign values to c_begin and c_end
        c_beg = gene_dict[key][0]
        c_end = gene_dict[key][1]
        #print("c_beg: ", c_beg)
        #print("c_end", c_end)
        #We compare each gene to one another in the temporary array, skipping when the genomes are the same
        for idx1, i in enumerate(sequences): 
            for idx2, j in enumerate(sequences): 
                if i != j: 
                    #count_matrix stores how many counts of each character in the ASCII table match in both sequences
                    count_matrix = [[0 for x in range(ASCII_LENGTH)] for x in range(ASCII_LENGTH)]
                    #compute_matches finds the count matrix for each sequences between the beginning and ending positions of the genes
                    compute_matches(i, j, count_matrix, dictionary, c_beg, c_end)
                    #find the substitution rate
                    substitution_rate = compute_substitution_rate(count_matrix)
                    #add the substitution rate to the distance matrix at the dorrect column for the species
                    edit_distance_matrix(distance_matrix, substitution_rate, idx1, idx2, dictionary)
        #print the distance matrices to a separate file for each gene
        print_distance_matrix(distance_matrix, key, dictionary)
        #cluster_zero_distances_test(distance_matrix, dictionary) 
        #construct the phylogenetic tree
        G = construct_graph(dictionary)
        #the file name is send to the phylo tree class; each gene will include the key name for reference
        filename = key + "_gene_tree.txt"
        #construct the phylogenetic tree
        px.Phylo_Tree(G, distance_matrix, dictionary, filename, "UPGMA") 
        
def return_gene_block(temp):
    sequences = []
    #An array to hold only the sequences of characters
    edited_sequences = []
    #Iterate through temp and extract 
    for element in temp: 
        arr = element.split()
        #Append just the genome info to the list, converting each string to a list in order to 
        #make the string object mutable. 
        sequences.append(list(arr[6]))

    #The reference genomes for the selected gene is the first in the sequence in this case. 
    for idx, char in enumerate(sequences[0]): 
        if char == "-":
            [x.pop(idx) for x in sequences]
        else: 
            continue
    #Convert list object back to a string object using the .join() method 
    for i, element in enumerate(sequences): 
        sequences[i] = "".join(element)
    return sequences 


def cluster_zero_distances_test(distance_matrix, dictionary): 
    inverse_dict = {v: k for k, v in dictionary.items()}
    flag = True
    cov_cluster = list(dictionary.keys())[0]

    idx = dictionary[cov_cluster]
    for i, col in enumerate(distance_matrix[idx]):
        if col == 0:
            print(col)
        else: 
            #print("species: ", inverse_dict[i])
            #print("species: ", inverse_dict[i])
            [x.pop(i) for x in distance_matrix]
            distance_matrix.pop(i)
            del dictionary[inverse_dict[i]]
            del inverse_dict[i]
            sorted_dict = {k: v for k, v in sorted(dictionary.items(), key = lambda item: item[1])} 
            dict_list = list(sorted_dict.keys())
            for idx, element in enumerate(dict_list): 
                dictionary[element] = idx
                inverse_dict[idx] = element
    print("end cluster function: ", dictionary)

def construct_graph(dictionary): 
    G = gx.Graph()
    G.add_nodes(dictionary.keys())
    for node in dictionary.keys():
        G.increment_sequence_count(node, 1)
    return G

def create_dict(sequences, dictionary):
    idx = 0 
    for seq in sequences: 
        x = seq.split()
        species = x[1]
        dictionary[species] = idx
        idx += 1

def compute_matches(seq1, seq2, count_matrix, dictionary, c_beg, c_end):
    #print("************************sequence begin************************")
    #arr1 = seq1.split()
    #arr2 = seq2.split()
    #print("seq1: ", arr1[1])
    #print("seq2: ", arr2[1])
    #print("c_beg: ", c_beg)
    #print("c_end: ", c_end)
    #sequence1 = arr1[6]
    #sequence2 = arr2[6]
    #_sequence1 = sequence1.replace("-", "")
    _sequence1 = seq1[c_beg:c_end]
    _sequence2 = seq2[c_beg:c_end]
    #print(_sequence1)
    #print(_sequence2)
    #_sequence2 = sequence2.replace("-", "")
    #_sequence2 = _sequence2[c_beg:c_end]
    #print("seq1: ", _sequence1[:100])
    #print("seq2: ", _sequence2[:100])
    #print("************************sequence begin************************")
    #print(sequence1)
    #print("************************sequence close************************")
    #print("************************sequence end************************")
    for i in range(len(_sequence1)): 
        count_matrix[ord(_sequence1[i].upper())][ord(_sequence2[i].upper())] += 1
    #print(count_matrix)
    return count_matrix

def compute_substitution_rate(counts): 
    matches = counts[65][65] + counts[67][67] + counts[71][71] + counts[84][84]
    mismatches = counts[65][67] + counts[65][71] + counts[65][84]
    + counts[67][65] + counts[67][71] + counts[67][84] + counts[71][65]
    + counts[71][67] + counts[71][84] + counts[84][65] + counts[84][67] + counts[84][71]
    #print("matches: ", matches)
    #print("mismatches: ", mismatches)
    if matches + mismatches == 0: 
        return 0
    else:
        return (mismatches / (mismatches + matches))
 
def edit_distance_matrix(distance_matrix, substitution_rate, i, j, dictionary):
    temp_sorted = sorted([i, j])
    distance_matrix[temp_sorted[0]][temp_sorted[1]] = substitution_rate

def print_distance_matrix(results, key, dictionary):
    #print("key: ", key)
    filename = key + "_distances.txt"
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    outfile = open(filename, "w+")
    df = pd.DataFrame(results, columns = [i for i in range(len(results))])
    print("KEY: ", key, '\n\n', file = outfile)
    counter = 0
    print(dictionary, file = outfile)
    print(df, '\n', file = outfile) 
    print(df)

def print_results(results):
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    df = pd.DataFrame(results, columns = [i for i in range(len(results))])
    print(df)

main() 
