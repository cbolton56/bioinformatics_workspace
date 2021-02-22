import os
import pandas as pd
import graph as gx
import phylo_tree as px

FILE_PATH = "/home/turing/mhou/csci652/data/sars-cov2/UCSC/wuhCor1.strainName119way.maf"
ASCII_LENGTH = 128

def main():
    infile = open(FILE_PATH)
    lines = infile.readlines()
    temp = []
    species = []
    #distance_matrix = [[0 for x in range(119)] for x in range(119)]
    dictionary = {}
    for line in lines:
        if line[0] == "s": 
            temp.append(line)
    distance_matrix = [[0 for x in range(len(temp))] for x in range(len(temp))]
    create_dict(temp, dictionary)
    print(dictionary)
    for i in temp: 
        for j in temp:
            if i != j: 
                count_matrix = [[0 for x in range(ASCII_LENGTH)] for x in range(ASCII_LENGTH)]
                compute_matches(i, j, count_matrix, dictionary)
                substitution_rate = compute_substitution_rate(count_matrix)
                edit_distance_matrix(distance_matrix, substitution_rate, i, j, dictionary)

    print_results(distance_matrix, dictionary)
    G = construct_graph(dictionary)
    px.Phylo_Tree(G, distance_matrix, dictionary, "full_genome_tree.txt", "UPGMA")

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

def compute_matches(seq1, seq2, count_matrix, dictionary):
    #print("************************sequence begin************************")
    arr1 = seq1.split()
    arr2 = seq2.split()
    #print("seq1: ", arr1[1])
    #print("seq2: ", arr2[1])
    sequence1 = arr1[6]
    sequence2 = arr2[6]
    #print("seq1: ", sequence1[:100])
    #print("seq2: ", sequence2[:100])
    #print("************************sequence begin************************")
    #print(sequence1)
    #print("************************sequence close************************")
    #print("************************sequence end************************")
    for i in range(len(sequence1)): 
        count_matrix[ord(sequence1[i].upper())][ord(sequence2[i].upper())] += 1
    #print(count_matrix)
    return count_matrix

def compute_substitution_rate(counts): 
    matches = counts[65][65] + counts[67][67] + counts[71][71] + counts[84][84]
    mismatches = counts[65][67] + counts[65][71] + counts[65][84]
    + counts[67][65] + counts[67][71] + counts[67][84] + counts[71][65]
    + counts[71][67] + counts[71][84] + counts[84][65] + counts[84][67] + counts[84][71]
    return (mismatches / (mismatches + matches)) 
 
def edit_distance_matrix(distance_matrix, substitution_rate, i, j, dictionary):
    arr1 = i.split()
    arr2 = j.split()
    species1 = arr1[1]
    species2 = arr2[1]
    idx1 = dictionary[species1]
    idx2 = dictionary[species2]
    temp_sorted = sorted([idx1, idx2])
    distance_matrix[temp_sorted[0]][temp_sorted[1]] += substitution_rate

def print_results(results, dictionary):
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    outfile = open("full_genome_distance_matrix.txt", "w+")
    df = pd.DataFrame(results, columns = [i for i in range(len(results))])
    print(df, file = outfile)
    print(df)
    print(dictionary)

main()
