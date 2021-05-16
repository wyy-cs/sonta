# Social netwOrks Nodes clusTering plAtform (SONTA)

# Input Files:
SONTA supports two format of input topology files: 1) edge list, and 2) adjacency list. All nodes will be renamed and the corresponding relationship between nodes' new id and original name is stored in a map<string, int> structure.

SONTA supports two format of input cluster ground truth files: 1) each line denotes a node's affilication, and 2) edge line denotes a certain cluster to which a set of nodes belong. Similarly, all clusters will be renamed and the corresponding relationship between clusters' new id and their original name is stored in a map<string, int> structure.

# Tips:


# Commands:
1. g++ sonta.cpp -o sonta -fopenmp -std=c++11
2. ./sonta name_dataset
