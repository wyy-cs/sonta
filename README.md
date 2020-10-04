# Social netwOrks Nodes clusTering plAtform (SONTA)

# Input Files:
1. directed graph file (.graph): each line contains node's all neighbors. (For example, the first line in the file is her neighbors set of node-1). The node id is continuous and begin from number 1. !!! Specifically, all graph is regarded as a directed graph (undirected can be regarded as a directed graph)
2. cluster file (.gt): each line contains all nodes in a cluster. (For example, the first line in the file is the nodes belonging to first cluster)
3. node feature file (.nf): each line contains a set of features a specific node has. (For example, the first line in the file is the a set of feature ids node-1 has.)
4. .graph topology file is necessary, while .gt and .nf files are not necessary.

# Tips:
1. All the one-dimensional array's first element's value is its number of all elements. (for example, if we store 10 values in an one-dimensional array A[11], where A[0] = 10). Although this operation will be at a cost of centain memory space, it will facilitate index and calculation.
2. For the two-dimensional array B, similarly the first one-dimensional array within is set as B[0] = new int[1]; B[0][0] = "number of one-dimensional vectors", other one-dimensional array's first element is also the correspoing number of elements (as Tips 1 shows).

# Commands:
1. g++ sonta.cpp -o sonta
2. ./sonta name_dataset
