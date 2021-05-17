# Social netwOrks Nodes clusTering plAtform (SONTA)

# Input graph topological files:
SONTA supports two format of input topology files. All nodes will be renamed and the correspondence between nodes' new id and original name is stored in a map<string, int> structure.

1) edge list.
********
NodeName1 NodeName2
NodeName3 NodeName4
...
...
...
*******

2) adjacency list. (first name is target node)
------------------------
TargetNodeName1 NeighborNodeName1 NeighborNodeName2 NeighborNodeName3 ...... NeighborNodeNameI
TargetNodeName2 NeighborNodeName4
TargetNodeName3 NeighborNodeName5 NeighborNodeName6 NeighborNodeName7 ...... NeighborNodeNameJ
...
...
...
------------------------

#  Input graph cluster ground truth files:
SONTA supports two format of input cluster ground truth files. Similarly, all clusters will be renamed and the correspondence between clusters' new id and their original name is stored in a map<string, int> structure.

1) each line denotes a node's affilication. (first name is target node, others are cluster name)
------------------------
NodeName1 ClusterName1 ClusterName2 ...... ClusterNameI
NodeName2 ClusterName3
NodeName3 ClusterName4 ClusterName5 ...... ClusterNameJ
...
...
...
------------------------

2) edge line denotes a certain cluster to which a set of nodes belong. (first name is target cluster, others are node name)
------------------------
ClusterName1 NodeName1 NodeName2 ...... NodeNameI
ClusterName2 NodeName3
ClusterName3 NodeName4 NodeName5 ...... NodeNameJ
...
...
...
------------------------

# Commands:
1. g++ sonta.cpp -o sonta -fopenmp -std=c++11
2. ./sonta name_dataset
