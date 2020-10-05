# Social netwOrks Nodes clusTering plAtform (SONTA)

# Input Files:
1. Directed graph file (.graph): each line contains a specifsic node's all neighbors' id. (e.g., the first line in the file is a set of node-1's all neighbors id.). The node id is required to be continuous and begin from number 1.!!! Specifically, all graph is regarded as a directed graph (undirected can be regarded as a directed graph where each linked node pair has two directed edges.)
2. Cluster file (.gt): each line contains all nodes within a specific cluster. (e.g., the first line in the file is a set of nodes belonging to the first cluster.)
3. node feature file (.nf): each line contains a set of feature ids a specific node has. (e.g., the first line in the file is a set of feature ids node-1 has.)
4. Graph topology file (.graph) is essential, its absence will cause the interruption of program; while (.gt) and (.nf) files are not essential.
5. All input files are binary, in which any two adjacent numbers in each line are separated by "Space", instead of "Tab".

# Tips:
1. All one-dimensional arrays' first element's value is the number of all elements with them. (e.g., if we store 10 values in an one-dimensional array A[11], then A[0] is assigned the value 10). Although this operation comes at the cost of a certain amount of space consumption, it will facilitate the indexing of vector elements.
2. For the two-dimensional arrays (e.g., B), the first one-dimensional array within them is set as B[0] = new int[1]; B[0][0] = "the number of one-dimensionaminusl arrays minus 1", other one-dimensional array's first element is also the correspoing number of elements (as Tips 1 shows).
3. Generally, for the three-dimensional arrays (e.g., C), the first element C[0] only has one element C[0][0], storing "the number of two-dimensional arrays minus 1".

# Commands:
1. g++ sonta.cpp -o sonta
2. ./sonta name_dataset
