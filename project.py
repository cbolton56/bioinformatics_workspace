"""
    Project 3
    z1885022 - Courtney Bolton

"""
import sys
import os
import pandas as pd
import graph as gx
import phylo_tree as px
import copy
#The length of the list to store characters when calculating substitution rate, generally is the length of ASCII table.
ASCII_LENGTH = 128

#File path for MAF files
FILE_PATH =  "/home/turing/mhou/csci652/data/guideTreeProject/pw-alignments"


def main():
    """Main routine handles file opening and calls subroutines. 
    
    """
    #A list for storing absolute path of all files in directory
    files = []
    counter = 0
    for path in os.listdir(FILE_PATH): 
        full_path = os.path.join(FILE_PATH, path)
        files.append(full_path)
    
    #Call create_dict( ) to make a dictionary object of all files based on their species name. This will be used to create a lookup table for the index of each species.
    dictionary = create_dict(files)
    #Results is the distance matrix for the substitution rate per species.
    distance_matrix = [[0 for x in range(len(dictionary))] for x in range(len(dictionary))]
    for sourcefile in files: 
    #Open infile in the directory and store all lines in the variable "lines"
        infile = open(sourcefile)
        lines = infile.readlines()

        #A temporary array for storing columns of the file containing genetic information
        temp = []
        #A 2D list for storing results for the pairwise alignment being sequenced.
        counts = [[0 for x in range(ASCII_LENGTH)] for x in range(ASCII_LENGTH)]
    
        #The "s" at the beginning of the lines in the file that contain genetic information is used as an indicator. 
        for line in lines: 
            #If the line contains "s" we store the information in an array and only use the last index of the array, which is the sequence of chars.
            if line[0] == "s":
                #Split on the whitespace 
                arr = line.split()
                #Append the genetic information to our temporary array
                temp.append(arr[6])
                #After gathering 2 sequences of characters (genetic info), we pass the sequences to compute substitution and transition rates.
                if len(temp) == 2:
                    compute_matches(temp[0], temp[1], counts)
                    #Clear the temporary array.
                    temp.clear()
        #Find the substitution rate of the particular pair of species
        substitution_rate = compute_substitution_rate(counts)
        write_to_distance_matrix(distance_matrix, dictionary, substitution_rate, sourcefile)
    #print_results(distance_matrix, dictionary)
    G = construct_graph(dictionary)
    tree1 = px.Phylo_Tree(G, distance_matrix, dictionary,"newest.txt",  "UPGMA")
    #tree2 = px.Phylo_Tree(G, distance_matrix, dictionary, "neighbor_tree.txt",  "NEIGHBOR")
            
def construct_graph(dictionary): 
    G = gx.Graph() 
    G.add_nodes(dictionary.keys())
    for node in dictionary.keys(): 
        G.increment_sequence_count(node, 1)
    return G

def write_to_distance_matrix(results, dictionary, substitution_rate, sourcefile):
    y = sourcefile.split("/")
    z = y[8].split(".")
    species_1 = z[0]
    species_2 = z[1]
    idx_1 = dictionary[species_1]
    idx_2 = dictionary[species_2]
    sorted_list = sorted([idx_1, idx_2])
    results[sorted_list[0]][sorted_list[1]] += substitution_rate

def create_dict(files): 
    """This function creates a dictionary object that stores all the species names and a reference to which index the particular species is being stored at.

    Args: files, a list of all files from which species pairs will be extracted. 
    Returns: @dictionary, a dictionary object where the key is the species name and the value is the index in the 2D results array where the sub rate is 
    stored for that species.
    """
    dictionary = {}
    index = 0
    for file in files:
        y = file.split("/")
        z = y[8].split(".")
        species_1 = z[0]
        species_2 = z[1]

        if species_1 not in dictionary.keys():
            dictionary[species_1] = index 
            index += 1

        if species_2 not in dictionary.keys():
            dictionary[species_2] = index
            index += 1
    return dictionary

def compute_matches(string1, string2, counts):
    """This function computes how many matches per sequence exist between two sequences. 

    Args: Two string objects, @string1 and @string2, which are the sequences of characters extracted from the 
    maf file. 
    @counts: A 2D list of each base pair
    Returns: @counts, a 2D list with updated count values
    """
    #Iterate over the length of the strings (string1 is used here), convert each character into ASCII value and update the value of the base pair.
    for i in range(len(string1)):
        counts[ord(string1[i].upper())][ord(string2[i].upper())] += 1     
    return counts

def compute_substitution_rate(counts):
    """This function computes substitution rate by adding the ASCII value. Since we don't include gaps or indels in the substitution rate the 
    AGCT values are accessed directly and summed together.

    Args: @counts: a 2D list of all of the mismatch and match counts for each ASCII value

    Returns: the substitution rate for the paired genomes.
    """
    #ASCII value based; A = 65, C = 67, etc.
    matches = counts[65][65] + counts[67][67] + counts[71][71] + counts[84][84]
    #Find all mismatches
    mismatches = counts[65][67] + counts[65][71] + counts[65][84] + counts[67][65] + counts[67][71]  + counts[67][84] + counts[71][65] + counts[71][67] + counts[71][84] + counts[84][65] + counts[84][67] + counts[84][71]
    #Sub rate is mismatches over matches plus mismatches
    return (mismatches / (mismatches + matches))

def print_results(results, dictionary):
    df = pd.DataFrame(results, columns= [i for i in range(len(results))])
    print(df)

main()
