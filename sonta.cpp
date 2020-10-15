#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <cstdlib>
#include <time.h> 
#include <stdlib.h>
#include <math.h>
#include <map>
#include <algorithm> 
#include <vector>
#include <sstream>
#include <sys/types.h>
#include <iterator>
#include <set>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>
//g++ -fopenmp -std=c++11
#ifdef USE_OPENMP
#include <omp.h>
#endif
#define RndUniDevInt(Range) (rand() % Range) // Randomly generate an integer between 0 and b, [0,b), do not contain b.
#define RndUniDev() (rand() / double(RAND_MAX)) // Ramdomly generate a float between 0 and 1, 0~1.
using namespace std;

typedef vector<int> TIntV;
typedef vector<double> TDblV;
typedef vector<string> TStrV;
typedef vector<vector<int> > TIntVV;
typedef vector<vector<double> > TDblVV;
typedef vector<set<int> > TIntSetV;
typedef vector<set<double> > TDblSetV;
typedef vector<pair<double, int> > TDblIntPrV;
typedef set<int> TIntSet;
typedef set<int>::iterator TIntSetIter;
typedef set<double> TDblSet;


// record the cost time from tod1 to tod2
long long int todiff(struct timeval* tod1, struct timeval* tod2) {
  long long t1, t2;
  t1 = tod1->tv_sec * 1000000 + tod1->tv_usec;
  t2 = tod2->tv_sec * 1000000 + tod2->tv_usec;
  return t1 - t2;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// basic graph class
class GRAPH {
private:
  // dataset name (do not contain the suffix name)
  string DataName;
  /*---------basic topology-----------*/
  // network topology, each set of corresponding node contains its all neighbors
  TIntSetV Neighbor;
  // number of nodes
  int NumNodes;
  // number of edges
  int NumEdges;
  /*--------weighted graph------------*/
  // values on nodes (each node is associated with a single value (e.g., centrality))
  TDblV WgtNodes;
  // values on nodes (each node is associated with a vector,0-1, binary)
  TIntVV BinNodeFea;
  // length of feature of nodes when each node is associated with a vector
  int LenNodeFea;
  /*-------cluster gt information---------*/
  // number of real clusters
  int NumClus; 
  // each set of corresponding node contains the cluster ids she belongs
  TIntSetV ClusInNode;
  // each set of corresponding cluster contains all nodes in this cluster
  TIntSetV ClusInClus;

public:
  // constructor (initialization, file name is needed)
  GRAPH(string dataname) { DataName = dataname; }
  // get data name
  string GetDataName() { return DataName; }
  
  // load a graph (directed) with topology information
  void LoadGraphTpl();
  // get graph topology
  TIntSet GetNeighbor(int NID) { return Neighbor[NID-1]; }
  // get number of nodes
  int GetNumNodes() { return NumNodes; }
  // get number of edges
  int GetNumEdges() { return NumEdges; }
  // get degree of a node
  int GetDeg(int NID) { return Neighbor[NID-1].size(); }
  // calculate the second-order proximity of node pair
  double CalSecOrdProx(int NID1, int NID2, int type);
  // calculate the second-order proximity of node pair (Jaccard Coefficient)
  double CalSecOrdProx(int NID1, int NID2) { return CalSecOrdProx(NID1, NID2, 1); }
  
  // get a node's weight (single value)
  double GetWgtNode(int NID) { return WgtNodes[NID-1]; }
  // load a graph node feature information
  void LoadNodeFea();
  // get length of features
  int GetLenNodeFea() { return LenNodeFea; }
  // get node features (in 0-1)
  TIntV GetBinNodeFea(int NID) { return BinNodeFea[NID-1]; }
  
  // load real cluster information
  void LoadClusGt();
  // get gt clusters listed in cluster
  TIntSetV GetClusInClus() { return ClusInClus; }
  // get a set of nodes in a specific cluster 
  TIntSet GetClusInClus(int CID) { return ClusInClus[CID-1]; }
  // get gt clusters listed in node
  TIntSetV GetClusInNode() { return ClusInNode; }
  // get a set of cluster ids to which a specific node belongs
  TIntSet GetClusInNode(int NID) { return ClusInNode[NID-1]; }
  // get number of clusters
  int GetNumClus() { return NumClus; }
  // calculate the number of common communities of all node pairs
  int CalNumCmnCmty(int NId1, int NId2) { return CalNumJointSets(ClusInNode[NId1-1], ClusInNode[NId2-1]); }
  
  // degree centrality of one node
  int DegCtly(int NID) { return GetDeg(NID); }
  // conductance centrality of one node
  double CondCtly(int NID);
  
  // visualize the graph (output a .gml file which can be fed into software named Cytoscape
  void VisualizeGraph();
  // test the functionality of this Graph class
  void FuncTest();
  
  // neighbor based proximity
  double NghBsdProx(TIntSet FirNode, TIntSet SecNode, int type);
  // calculate the number of elements in the intersection set of two sets
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
};

void GRAPH::LoadGraphTpl() {
  string SGraphFile = DataName + ".graph";
  const char* GraphFile = SGraphFile.c_str();
  if (access(GraphFile, R_OK|W_OK) != 0) {
    printf("No topology (.graph) file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finG;
  finG.open(GraphFile);
  string szLine;
  TIntSet tempset;
  while(getline(finG, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      tempset.insert(atoi((*iter).c_str()));
    }
    Neighbor.push_back(tempset);
    tempset.clear();  
  }
  finG.close();
  finG.clear();
  NumNodes = Neighbor.size();
  NumEdges = 0;
  for (int i = 0; i < Neighbor.size(); i++) { NumEdges += Neighbor[i].size(); }
  // each edge is twicely recorded (undirected graph)
  NumEdges /= 2;
}

// load a graph node feature information
void GRAPH::LoadNodeFea() {
  string SFeaFile = DataName + ".nf";
  const char* FeaFile = SFeaFile.c_str();
  if (access(FeaFile, R_OK|W_OK) != 0) {
    cout << "There is no node feature (.nf) file!!!" << endl;
    return;
  }
  // if there exists a .fea file
  int maxFeaId = 0;
  int tempd;
  ifstream finNF;
  finNF.open(FeaFile);
  string szLine; 
  while(getline(finNF,szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<TStrV >(tData));
    // find the maximum feature id
    for(TStrV::iterator iter = tData.begin(); iter!=tData.end(); ++iter) {
      tempd = atoi((*iter).c_str());
      if (maxFeaId <= tempd) { maxFeaId = tempd; }
    }  
  }
  finNF.close();
  finNF.clear();
  // length of each node's features
  LenNodeFea = maxFeaId;
  cout << "Length of node features: " << LenNodeFea << endl;
  BinNodeFea.resize(NumNodes);
  int i, j;
  for (i = 0; i < NumNodes; i++) { 
    for (j = 0; j < LenNodeFea; j++) { 
    BinNodeFea[i].push_back(0); 
  }
  }
  i = 0;
  finNF.open(FeaFile);
  while(getline(finNF,szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<TStrV >(tData));
    // reopen and write to fea_nodes
    // if this node has a feature, then the value of corresponding position is 1, 0 otherwise.
    for(TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) {
      BinNodeFea[i][atoi((*iter).c_str()) - 1] = 1;
    }
    i++;  
  }
  finNF.close();
  finNF.clear();
}

// load real cluster information
void GRAPH::LoadClusGt() {
  string SClusFile = DataName + ".gt";
  const char* ClusFile = SClusFile.c_str();
  if (access(ClusFile, R_OK|W_OK) != 0) { printf("No cluster (.gt) file exists!!!\n"); return; }
  ifstream finGT;
  finGT.open(ClusFile);
  string szLine;
  TIntSet tempset;
  while(getline(finGT, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) {
      tempset.insert(atoi((*iter).c_str()));
    }
    ClusInClus.push_back(tempset);
    tempset.clear();
  }
  finGT.close();
  finGT.clear();
  NumClus = ClusInClus.size();
  ClusInNode.resize(NumNodes);
  for (int i = 1; i <= ClusInClus.size(); i++) {
    for (TIntSetIter it_set = ClusInClus[i - 1].begin(); it_set != ClusInClus[i - 1].end(); it_set++) {
      ClusInNode[*it_set - 1].insert(i);
    }
  }
}

// calculate the second-order proximity
double GRAPH::CalSecOrdProx(int NID1, int NID2, int type) {
  if (type == 1) { printf("Second-order proximity: Jaccard Coefficient.\n"); }
  else if (type == 2) { printf("Second-order proximity: Salton Index.\n"); }
  else if (type == 3) { printf("Second-order proximity: Sorensen Index.\n");  }
  else if (type == 4) { printf("Second-order pairroximity: Hub Promoted Index.\n"); }
  else if (type == 5) { printf("Second-order pairroximity: Hub Depressed Index.\n"); }
  else if (type == 6) { printf("Second-order pairroximity: Leicht-Holme-Newman Index.\n"); }
  double Proximity = NghBsdProx(Neighbor[NID1-1], Neighbor[NID2-1], type);
  return Proximity;
}

// visualize the graph (output a .gml file which can be fed into software named Cytoscape
void GRAPH::VisualizeGraph() {
  string SGmlFile = DataName + ".gml";
  const char* GmlFile = SGmlFile.c_str();
  ofstream foutG;
  foutG.open(GmlFile);
  foutG << "graph" << endl;
  foutG << "[" << endl;
  foutG << "  directed 0" << endl;
  int i, j;
  TIntSetIter it_set;
  for (i = 0; i < NumNodes; i++) {
    foutG << "  node" << endl;
    foutG << "  [" << endl;
    foutG << "    id " << i + 1 << endl;
    foutG << "    label " << "\"" << i + 1 << "\"" << endl;
    // non-overlapping clusters
    // if one node does not belong to any cluster, the flag equals 0;
    it_set = ClusInNode[i].begin();
    foutG << "    cluster " << "\"" << *it_set << "\"" << endl;
    foutG << "  ]" << endl;
  }
  
  for (i = 0; i < NumNodes; i++) {
  for (it_set = Neighbor[i].begin(); it_set != Neighbor[i].end(); it_set++) {
      if (i < *it_set) {
        foutG << "  edge" << endl;
        foutG << "  [" << endl;
        foutG << "    source " << i + 1 << endl;
        foutG << "    target " << *it_set << endl;
        foutG << "  ]" << endl;
      }
    }
  }
  foutG << "]" << endl;
  foutG.close();
}

// test the functionality of this Graph class
void GRAPH::FuncTest() {
  if (NumNodes > 50) {
    printf("Data scale is too big not suitable for full-print, please specilize it!!!!\n");
    return;
  }
  int i, j;
  TIntSetIter it_set;
  if (Neighbor.empty()) { 
    printf("No topology information!!!\n");
    return;  
  } else {
    printf("============== Number of Nodes: %d, Number of edges: %d =================\n", NumNodes, NumEdges);
    for (i = 0; i < Neighbor.size(); i++) {
      printf("node-%d: ", i);
      for (it_set = Neighbor[i].begin(); it_set != Neighbor[i].end(); it_set++) {
        printf("%d ", *it_set);
      }
      printf("\n");
  }
  }
  
  if (BinNodeFea.empty()) { 
    printf("No node feature information!!!\n"); 
  } else {
    printf("================= Length of node feature: %d =====================\n", LenNodeFea);
    for (i = 0; i < NumNodes; i++) {
      printf("node-%d: ", i);
      for (j = 0; j < LenNodeFea; j++) {
        printf("%d ", BinNodeFea[i][j]);
      }
      printf("\n");
    }
  }
  
  if (ClusInClus.empty()) { 
    printf("No ground-truth information!!!\n"); 
  } else {
    printf("============ Number of clusters: %d ==================\n", NumClus);
    for (i = 0; i < ClusInClus.size(); i++) {
      printf("cluster-%d: ", i + 1);
      for (it_set = ClusInClus[i].begin(); it_set != ClusInClus[i].end(); it_set++) {
        printf("%d ", *it_set);
      }
      printf("\n");
    }
  }

  VisualizeGraph();
}

// conductance centrality (Gleich, et al. KDD12)
double GRAPH::CondCtly(int NID) {
  double ImptNode = 0.0; // small value indicates an important node.
  TIntSet NghNId(Neighbor[NID - 1]);
  NghNId.insert(NID); // did in SNAP
  int LenNghs = NghNId.size();
  if (LenNghs < 5) {
    ImptNode = 1.0;
	return ImptNode;
  }
  int Edges2 = 2 * NumEdges;
  int Vol = 0,  Cut = 0;
  for (TIntSetIter it_set = NghNId.begin(); it_set != NghNId.end(); it_set++) {
    for (TIntSetIter it_set1 = Neighbor[*it_set - 1].begin(); it_set1 != Neighbor[*it_set - 1].end(); it_set1++) {
      // whether her neighbor's neighbor is also her neighbor.
      if (NghNId.find(*it_set1) == NghNId.end()) { Cut += 1; }
    }
    // Vol store the summation of degree of all nodes inside this set 
    Vol += Neighbor[*it_set - 1].size();
  }
  // get conductance
  if (Vol != Edges2) {
    if (2 * Vol > Edges2) { ImptNode = Cut / double (Edges2 - Vol); }
    else if (Vol == 0) { ImptNode = 0.0; }
    else { ImptNode = Cut / double(Vol); }
  } else {
    if (Vol == Edges2) { ImptNode = 1.0; }
  }
  return ImptNode;
}

// neighbor based proximity
double GRAPH::NghBsdProx(TIntSet FirNode, TIntSet SecNode, int type) {
  int NumJoint = CalNumJointSets(FirNode, SecNode);
  int LenFirNode = FirNode.size();
  int LenSecNode = SecNode.size();
  int NumUnion = LenFirNode + LenSecNode - NumJoint;
  NumJoint += 2;
  NumUnion += 2;
  LenFirNode++;
  LenSecNode++;
  // Jaccard Coefficient
  if (type == 1) { return (double)NumJoint / NumUnion; }
  // Salton Index
  else if (type == 2) { return (double)NumJoint / (sqrt(LenFirNode * LenSecNode)); }
  // Sorensen Index
  else if (type == 3) { return 2.0 * NumJoint / (LenFirNode + LenSecNode); }
  // Hub Promoted Index
  else if (type == 4) {
    if (LenFirNode < LenSecNode) { return 2.0 * NumJoint / LenFirNode; }
    else { return 2.0 * NumJoint / LenSecNode; }
  }
  // Hub Depressed Index
  else if (type == 5) {
    if (LenFirNode < LenSecNode) { return 2.0 * NumJoint / LenSecNode; }
    else { return 2.0 * NumJoint / LenFirNode; }
  }
  // Leicht-Holme-Newman Index
  else if (type == 6) { return 2.0 * NumJoint / LenFirNode / LenSecNode; }
  return 0.0;
}


/////////////////////////////////////////////
//Class of performance test on nodes clustering
/*Tips:
1. The cluster structure of calculating AvgF1, NMI, and ARI is the second Type (each line denotes a set of nodes in a common community)
2. The cluster structure of calculating Omega Index is the first Type (each line denotes a set of clusters a specific node belongs)
3. The number of nodes is needed when calculate NMI, Omega Index, and ARI.
*/
class ClusterTest {
public: 
  // no ground-truth information 
  // Modularity of non-overlapping cluster strucuture. (Newman, et al. PRE, 2004)
  // Our1 and Our2 are same cluster structure with different format.
  // each set in Our1 denotes a cluster, containing a set of nodes who belong it.
  // each set in Our2 is a set of clusters to which a specific node belongs.
  double CalNonOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges);
  // Modularity of overlapping cluster structure (Shen, et al. Physica A, 2009)
  double CalOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges);
  // Tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
  double CalTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor);
  // Adjusted tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
  // penalizing very small and very large clusters and produces well-balanced solutions
  double CalAdjTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor);
  
public:
  // has ground-truth (.gt file) information
  // 1. Average of f1-score (AvgF1)
  double CalF1(TIntSet SetA, TIntSet SetB);
  double CalMaxF1(TIntSet SetA, TIntSetV TargetClus);
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalAvgF1(TIntSetV Our, TIntSetV GT);
  // 2. Normalized mutual information (NMI)
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalNMI(TIntSetV Our, TIntSetV GT, int num_nodes);
  // 3. Omega Index
  // each line in ClusInNode_Our and ClusInNode_GT is a set of clusters to which a specific node belongs
  double CalOmegaIndex(TIntSetV ClusInNode_Our, TIntSetV ClusInNode_GT, int num_nodes);
  // 4. Adjusted Rand Index (ARI)
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalARI(TIntSetV Our, TIntSetV GT, int num_nodes);
  
  // calculate the number of elements in the intersection set of two sets
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
};

// calculate modularity of non-overlapping cluster structure. (Newman, PRE, 2004.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// if input is a overlapping cluster structure, this function will only take the first cluster for each node as her FuzzyMembership.
// num_edges is the number of all edges in the graph.
double ClusterTest::CalNonOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
  int i = 0, k = 0;
  // number of all nodes in the graph
  int num_nodes = neighbor.size();
  TDblV theta;
  theta.resize(Our1.size());
  TIntSetIter it_set;
  for (k = 0; k < theta.size(); k++) {
    theta[k] = 0.0;
    for (it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      theta[k] += neighbor[*it_set - 1].size();
    }
    theta[k] /= 2.0 * num_edges;
  }
  // modularity of single node
  TDblV singleModul;
  singleModul.resize(num_nodes);
  double tempValue = 0.0, SumModul = 0.0;
  for (i = 0; i < num_nodes; i++) {
    tempValue = 0.0;
    singleModul[i] = 0.0;
    // determine whether node-i does not belong to any cluster.
    if (Our2[i].empty()) { 
      singleModul[i] = 0.0; 
    } else if (neighbor[i].empty() == 0) {
      it_set = Our2[i].begin();
      tempValue = CalNumJointSets(neighbor[i], Our1[*it_set - 1]) / (double)neighbor[i].size();
      singleModul[i] = (tempValue - theta[*it_set - 1]) * neighbor[i].size() / 2.0 / num_edges;
    }
    // summation
    SumModul += singleModul[i];
  }
  return SumModul;
}

// calculating modularity of overlapping and hierarchical cluster structure.(Shen, et al. Physica A, 2009.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// num_edges is the number of all edges in the graph.
double ClusterTest::CalOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
  int i = 0, j = 0, k = 0;
  int id_i, id_j;
  int num_clusters = Our1.size();
  TDblV singleModul;
  singleModul.resize(num_clusters);
  double SumModul = 0.0;
  for (k = 0; k < num_clusters; k++) {
    if (Our1[k].size() < 2) { singleModul[k] = 0.0; continue; }
    singleModul[k] = 0.0;  
    TIntSet VerifySet;
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      id_i = *it_set;
      VerifySet.insert(id_i);
      if (Our2[id_i - 1].empty()) { continue; }
      for (TIntSetIter it_set1 = Our1[k].begin(); it_set1 != Our1[k].end(); it_set1++) {
        if (VerifySet.find(*it_set1) != VerifySet.end()) { continue; }
        else {
          id_j = *it_set1;
          if (Our2[id_j - 1].empty()) { continue; }
          if (neighbor[id_i - 1].find(id_j) == neighbor[id_i - 1].end()) {
            singleModul[k] = -neighbor[id_i - 1].size() * neighbor[id_j - 1].size() / 2.0 / (double)num_edges / Our2[id_i - 1].size() / Our2[id_j - 1].size();
          } else {
            singleModul[k] = (1.0 - neighbor[id_i - 1].size() * neighbor[id_j - 1].size() / 2.0 / (double)num_edges) / Our2[id_i - 1].size() / Our2[id_j - 1].size();
          }
        }
      }
    }
    VerifySet.clear();
    SumModul += singleModul[k];
  }
  SumModul = 0.5 * SumModul / num_edges;
  return SumModul;
}

// calculating tightness of cluster structure (Bu, et al. Information Fusion, 2017)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
double ClusterTest::CalTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
  int id_i, id_j;
  int flag = 0;
  double FirTerm = 0.0, SecTerm = 0.0;
  int num_clusters = Our1.size();
  int num_nodes = neighbor.size();
  double SumTgt = 0.0;
  TDblV SinglTgt;
  SinglTgt.resize(num_clusters);
  TIntV numInEdges;
  numInEdges.resize(num_clusters);
  TIntV numOutEdges;
  numOutEdges.resize(num_clusters);
  // a cluster k
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      id_i = *it_set;
      // a neighbor of node id_i
      for (TIntSetIter it_set1 = neighbor[id_i - 1].begin(); it_set1 != neighbor[id_i - 1].end(); it_set1++) {
        id_j = *it_set1;
        TIntSet ngtClust(Our2[id_j - 1]);
        // j's FuzzyMemberships, if she belongs cluster k, thus numInEdges[k] plus 1.
        // it can be applied to non-overlapping as well as overlapping clusters
        flag = 0;
        for (TIntSetIter it_set1 = ngtClust.begin(); it_set1 != ngtClust.end(); it_set1++) {
          if ((*it_set1) == (k+1)) { numInEdges[k]++; }
          else if (flag != 1) {
            numOutEdges[k]++;
            // guarantee numOutEdges[k] is counted once for node-j,
            // when node-j belongs to several outter clusters besides cluster-k 
            flag = 1;
          }
        }
      }
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
    if (Our1[k].empty() == 0 && Our1[k].size() != num_nodes) {
      FirTerm = 2.0 * numInEdges[k] / Our1[k].size() / (double)Our1[k].size();
      SecTerm = numOutEdges[k] / Our1[k].size() / (double)(num_nodes - Our1[k].size());
      SinglTgt[k] = FirTerm - SecTerm;
    }
    SumTgt += SinglTgt[k];
  }
  return SumTgt;
}

// calculating adjusted tightness of cluster structure (Bu, et al. Information Fusion, 2017)
// penalizing very small and very large clusters and produces well-balanced solutions
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
double ClusterTest::CalAdjTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
  int id_i, id_j;
  int flag = 0;
  double FirTerm = 0.0;
  int num_clusters = Our1.size();
  int num_nodes = neighbor.size();
  double SumTgt = 0.0;
  TDblV SinglTgt;
  SinglTgt.resize(num_clusters);
  TIntV numInEdges;
  numInEdges.resize(num_clusters);
  TIntV numOutEdges;
  numOutEdges.resize(num_clusters); 
  // a cluster k
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      id_i = *it_set;
      // a neighbor of node id_i
      for (TIntSetIter it_set1 = neighbor[id_i - 1].begin(); it_set1 != neighbor[id_i - 1].end(); it_set1++) {
        id_j = *it_set1;
        TIntSet ngtClust(Our2[id_j - 1]);
        // j's FuzzyMemberships, if she belongs cluster k, thus numInEdges[k] plus 1.
        // it can be applied to non-overlapping as well as overlapping clusters
        flag = 0;
        for (TIntSetIter it_set1 = ngtClust.begin(); it_set1 != ngtClust.end(); it_set1++) {
          if ((*it_set1) == (k+1)) { numInEdges[k]++; }
          else if (flag != 1) {
            numOutEdges[k]++;
            // guarantee numOutEdges[k] is counted once for node-j,
            // when node-j belongs to several outter clusters besides cluster-k 
            flag = 1;
          }
        }
      }
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
	if (Our1[k].empty() == 0) {
      FirTerm = 2.0 * (num_nodes - Our1[k].size()) * numInEdges[k] /  (double)Our1[k].size();
      SinglTgt[k] = FirTerm - numOutEdges[k];
	}
    SumTgt += SinglTgt[k];
  }
  return SumTgt;
}

// calculate F1
double ClusterTest::CalF1(TIntSet SetA, TIntSet SetB) {
  int K_A = SetA.size();
  int K_B = SetB.size();
  int NumJoint = CalNumJointSets(SetA, SetB);
  double precision = 0.0, recall = 0.0;
  double f1 = 0.0;
  if (K_A != 0 && K_B != 0) {
    precision = (double)NumJoint / K_A;
    recall = (double)NumJoint / K_B;
  }
  if ((precision + recall) != 0 ) {
    f1 = 2.0 * precision * recall / (precision + recall);
    return f1;
  } else { return f1; }
}

// calculate maximum of F1
double ClusterTest::CalMaxF1(TIntSet SetA, TIntSetV TargetClus) {
  int K_A = SetA.size();
  int Size_T = TargetClus.size();
  double maxf1 = 0.0;
  for (int q = 0; q < Size_T; q++) {
    // double ClusterTest::CalF1(int* A, int* B)
    double tempf1 = CalF1(SetA, TargetClus[q]);
    if (tempf1 >= maxf1) { maxf1 = tempf1; }
  }
  return maxf1;
}

// calculating the average of f1-score
double ClusterTest::CalAvgF1(TIntSetV Our, TIntSetV GT) {
  int K_1 = Our.size();
  int K_2 = GT.size();
  TDblV MaxF1_1;
  MaxF1_1.resize(K_1);
  TDblV MaxF1_2;
  MaxF1_2.resize(K_2);
  int p = 0, q = 0;
  double sumMAXF1_1 = 0.0, sumMAXF1_2 = 0.0;
  for (p = 0; p < K_1; p++) {
    MaxF1_1[p] = CalMaxF1(Our[p], GT);
    sumMAXF1_1 += MaxF1_1[p];
  }
  for (q = 0; q < K_2; q++) {
    MaxF1_2[q] = CalMaxF1(GT[q], Our);
    sumMAXF1_2 += MaxF1_2[q];
  }
  double AvgF1 = sumMAXF1_1 / 2.0 / K_1 + sumMAXF1_2 / 2.0 / K_2;
  return AvgF1;
}

// calculating the normalized mutual information (NMI)
// [S. Fortunato. Community detection in graphs. Physics Reports, 486(3-5):75 – 174, 2010.]
double ClusterTest::CalNMI(TIntSetV Our, TIntSetV GT, int num_nodes) {
  int LenA = Our.size();
  int LenB = GT.size();
  int i = 1, j = 1;
  TDblVV number_joint;
  number_joint.resize(LenA);
  for (i = 0; i < LenA; i++) {
    number_joint[i].resize(LenB);
    for (j = 0; j < LenB; j++) {
      number_joint[i][j] = (double)CalNumJointSets(Our[i], GT[j]);
    }
  }
  double NMI = 0.0;
  TDblV NMI_A;
  NMI_A.resize(LenA);
  TDblV NMI_B;
  NMI_B.resize(LenB);
  for (i = 0; i < LenA; i++) { NMI_A[i] = 0.0; }
  for (j = 0; j < LenB; j++) { NMI_B[j] = 0.0; }
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) { NMI_A[i] += number_joint[i][j]; }
  }
  for (j = 0; j < LenB; j++) {
    for (i = 0; i < LenA; i++) { NMI_B[j] += number_joint[i][j]; }
  }
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) {
      if (number_joint[i][j] != 0.0) {
        NMI += number_joint[i][j] * log(number_joint[i][j] * num_nodes / NMI_A[i] / NMI_B[j]) / log(2);
      }
    }
  }
  NMI = NMI * (-1 * 2.0);
  double sum_NMI_A = 0.0;
  double sum_NMI_B = 0.0;
  for(i = 0; i < LenA; i++) {
    if (NMI_A[i] != 0.0) {
      sum_NMI_A += NMI_A[i] * log(NMI_A[i] / num_nodes) / log(2);
    }
  }
  for(j = 0; j < LenB; j++) {
    if(NMI_B[j] != 0) { sum_NMI_B += NMI_B[j] * log(NMI_B[j] / num_nodes) / log(2); }    
  }
  NMI = NMI / (sum_NMI_A + sum_NMI_B);
  return NMI;
}

// calculating the value of Omega index
// [S. Gregory. Fuzzy overlapping communities in networks. J. of Stat. Mech.: Theory and Experiment, 2011.]
// each set in ClusInNode_Our and ClusInNode_GT denotes a set of clusters to which a specific node belongs.
double ClusterTest::CalOmegaIndex(TIntSetV ClusInNode_Our, TIntSetV ClusInNode_GT, int num_nodes) {
  int SumTemp = 0;
  double results = 0.0;
  for (int i = 0; i < num_nodes; i++) {
    for (int j = i + 1; j < num_nodes; j++) {
      if (CalNumJointSets(ClusInNode_Our[i], ClusInNode_Our[j]) == CalNumJointSets(ClusInNode_GT[i], ClusInNode_GT[j])) { SumTemp++; }
    }
  }
  results = (double)SumTemp / num_nodes / num_nodes;
  return results;
}

// calculating Adjusted Rand Index (ARI)
// each set in Our1 and GT denotes a cluster, containing a set of nodes who belong it.
double ClusterTest::CalARI(TIntSetV Our, TIntSetV GT, int num_nodes) {
  int LenA = Our.size();
  int LenB = GT.size();
  int i, j;
  double numerator = 0.0, denominator = 0.0;
  double TempVal = 0.0;
  double numerator_1 = 0.0, numerator_2 = 0.0;
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) {
	  TempVal = CalNumJointSets(Our[i], GT[j]);
	  numerator_1 = numerator_1 + TempVal * (TempVal - 1.0) / 2.0;
	}
  }		
  double TempVal_1 = 0.0;
  for (i = 0; i < LenA; i++) {
    TempVal_1 += Our[i].size() * (Our[i].size() - 1.0) / 2.0;
  }
  double TempVal_2 = 0.0;
  for (j = 0; j < LenB; j++) {
    TempVal_2 += GT[j].size() * (GT[j].size() - 1.0) / 2.0;
  }
  numerator_2 = (TempVal_1 + TempVal_2) * 2.0 / num_nodes / (num_nodes - 1);
  denominator = (TempVal_1 + TempVal_2) / 2.0 - numerator_2;
  double results = (numerator_1 - numerator_2) / denominator;
  return results;
}

///////////////////////////////////
///////////////////////////////////////////////
class SIMGT {
private:
  GRAPH G;
  int NumComs; // number of communities
  TDblVV F; // membership matrix (size: NumNodes * NumClus)
  TDblV SumF; // community size (Sum_i F_ic for community c)
  
public:
  SIMGT(GRAPH graph): G(graph.GetDataName()) {}
  void Initialization(); // load data
  void NeighborComInit(); // initialize F with best neighborhood communities
  void PeformSIMGT(); // perform SIMGT
  
  void inline UpdateCom(int NID, int CID, double Val) {
    // update F and SumF
    SumF[CID-1] -= F[NID-1][CID-1];
    F[NID-1][CID-1] = Val;
    SumF[CID-1] += Val;
  }
};

void SIMGT::Initialization() {
  G.LoadGraphTpl(); // load .graph file (pure directed graph)
  printf("Number of Nodes: %d, Number of edges: %d\n", G.GetNumNodes(), G.GetNumEdges());
  G.LoadClusGt(); // load .gt file (cluster ground-truth)
  printf("Number of clusters: %d\n", G.GetNumClus());
}

void SIMGT::NeighborComInit() {
  // initialize F with best neighborhood communities (Gleich et.al. KDD'12)
  int i, j;
  F.resize(G.GetNumNodes());
  NumComs = G.GetNumClus();
  for (j = 0; j < NumComs; j++) { SumF.push_back(0.0); }
  for (i = 0; i < F.size(); i++) {
	for (j = 0; j < NumComs; j++) { F[i].push_back(0.0); }
  }
  TDblIntPrV NIdPhiV;
  TIntSet InvalidNIDS;
  // compute conductance of neighborhood community
  for (int NId = 1; NId <= F.size(); NId++) {
	double Phi = G.CondCtly(NId);
    NIdPhiV.emplace_back(make_pair(Phi, NId));
  }
  sort(NIdPhiV.begin(), NIdPhiV.end());
  printf("conductance computation completed\n");
  // choose nodes with local minimum in conductance
  int CurCID = 1;
  for (int ui = 0; ui < NIdPhiV.size(); ui++) {
    int UID = NIdPhiV[ui].second;
    if (InvalidNIDS.find(UID) != InvalidNIDS.end()) { continue; }
    // add the node and its neighbors to the current community
	UpdateCom(UID, CurCID, 1.0);
	TIntSet Ngh(G.GetNeighbor(UID)); // G.GetNeighbor(UID).begin() is a bug. (abandoned)
	for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
	  UpdateCom(*it_set, CurCID, 1.0);
	  // exclude its neighbors from the next considerations
	  InvalidNIDS.insert(*it_set);
    }
	CurCID++;
	if (CurCID > NumComs) { break; }
  }
  if (NumComs >= CurCID) {
    printf("%d communities needed to fill randomly\n", NumComs - CurCID + 1);
  }
  //assign a member to zero-member community (if any)
  for (int CID = 1; CID <= SumF.size(); CID++) {
    if (SumF[CID-1] == 0.0) {
      int ComSz = 10;
      for (int u = 0; u < ComSz; u++) {
        int UID = RndUniDevInt(G.GetNumNodes()) + 1;
        UpdateCom(UID, CID, RndUniDev());
      }
    }
  }
}

void SIMGT::PeformSIMGT() {
  Initialization();
  NeighborComInit();
  for (int i = 0; i < F.size(); i++) {
	  for (int j = 0; j < G.GetNumClus(); j++) {
		  cout << F[i][j] << " ";
	  }
	  cout << endl;
  }
}


int main(int argc, char **argv)
{
  string head = argv[1]; // name of dataset
  //int TypeProx = atoi(argv[2]);
  //float stepsize = atof(argv[3]);
  //int maxiter = atoi(argv[4]);
  //float delta = atof(argv[5]);
  
  GRAPH myGraph(head); // initialize an object
  struct timeval tod1, tod2; // record time
// Tools: SIMGT
  cout << "========= " << "Node clustering by SIMGT." << "================== " << endl;
  gettimeofday(&tod1, NULL);
  SIMGT simgt(myGraph);
  simgt.PeformSIMGT();
  gettimeofday(&tod2, NULL);
  cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
  cout << endl;
 /*
//Tools: Performance Test on nodes clustering
  cout << "========= " << "Performance test on nodes clustering." << "================== "<< endl;
  ClusterTest mytest;
  gettimeofday(&tod1, NULL);
  // no gt
  cout << "##### no ground-truth information:" << endl;
  cout << "Modularity (no overlapping) = " << mytest.CalNonOverlapModul(mygraph.GetClusInClus(), mygraph.GetClusInNode(), mygraph.GetNeighbor(), mygraph.GetNumEdges()) << endl;
  cout << "Modularity (overlapping) = " << mytest.CalOverlapModul(mygraph.GetClusInClus(), mygraph.GetClusInNode(), mygraph.GetNeighbor(), mygraph.GetNumEdges()) << endl;
  cout << "Tightness = " << mytest.CalTgt(mygraph.GetClusInClus(), mygraph.GetClusInNode(), mygraph.GetNeighbor()) << endl;
  cout << "Adjusted Tightness = " << mytest.CalAdjTgt(mygraph.GetClusInClus(), mygraph.GetClusInNode(), mygraph.GetNeighbor()) << endl;
  cout << endl;
  // with gt
  cout << "##### with ground-truth information:" << endl;
  cout << "AvgF1 = " << mytest.CalAvgF1(mygraph.GetClusInClus(), mygraph.GetClusInClus()) << endl;
  cout << "NMI = " << mytest.CalNMI(mygraph.GetClusInClus(), mygraph.GetClusInClus(), mygraph.GetNumNodes()) << endl;
  cout << "Omega Index = " << mytest.CalOmegaIndex(mygraph.GetClusInNode(), mygraph.GetClusInNode(), mygraph.GetNumNodes()) << endl;
  cout << "ARI = " << mytest.CalARI(mygraph.GetClusInClus(), mygraph.GetClusInClus(), mygraph.GetNumNodes()) << endl;
  gettimeofday(&tod2, NULL);
  cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
  cout << endl;*/

  return 0;
}
