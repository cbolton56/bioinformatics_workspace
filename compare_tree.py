import re
import dendropy
import glob as glob
from dendropy.calculate import treecompare
from itertools import combinations


def main():
    #Files are all stored in a files
    files = glob.glob('./*.txt')

    tree_combinations =  list(combinations([i for i in range(len(files))], 2))

    labels = []
    scores = []
    
    for tup in tree_combinations:
        #Name of gene 1
        x = files[tup[0]]
        x = x.replace(".", "")
        x = x.replace("\\", "")
        x = x.replace(".txt", "")
        x = x.replace("txt", "")
        x = x.replace("_", " ")
        x = x.replace("tree", "")
        #Name of gene 2
        y = files[tup[1]]
        y = y.replace(".", "")
        y = y.replace("\\", "")
        y = y.replace(".txt", "")
        y = y.replace("txt", "")
        y = y.replace("_", " ")
        y = y.replace("tree", "")

        label = x + "vs " + y
        print(label)
        infile_1 = open(files[tup[0]])
        infile_2 = open(files[tup[1]])
        lines1 = str(infile_1.readline())
        lines2 = str(infile_2.readline())
        s1 = lines1
        s2 = lines2

        # establish common taxon namespace
        tns = dendropy.TaxonNamespace()
        # ensure all trees loaded use common namespace
        tree1 = dendropy.Tree.get(
                data=s1,
                schema='newick', 
                preserve_underscores=True, 
                suppress_internal_node_taxa=False,
                taxon_namespace=tns)
        tree2 = dendropy.Tree.get(
                data=s2,
                schema='newick',
                preserve_underscores=True,
                suppress_internal_node_taxa = False,
                taxon_namespace=tns)
        ## Unweighted Robinson-Foulds distance
        score = treecompare.symmetric_difference(tree1, tree2)
        print("Comparing tree  ", files[tup[0]], " and ", files[tup[1]],  ": ", score)
        labels.append(label)
        scores.append(score)
    print(scores)
    print(labels)

main()
