// g++ sonta.cpp -o sonta -fopenmp -std=c++11 
// ./sonta DataName
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cstdlib>
#include <time.h> 
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#define RndUniDevInt(Range) (rand() % Range) // Randomly generate an integer between 0 and b, [0,b), do not contain b.
#define RndUniDev() (rand() / double(RAND_MAX)) // Ramdomly generate a float between 0 and 1, 0~1.
using namespace std;

typedef vector<int> TIntV;
typedef vector<double> TDblV;
typedef vector<string> TStrV;
typedef vector<bool> TBoolV;
typedef vector<vector<int> > TIntVV;
typedef vector<vector<double> > TDblVV;
typedef vector<set<int> > TIntSetV;
typedef vector<set<double> > TDblSetV;
typedef vector<pair<double, int> > TDblIntPrV;
typedef set<int> TIntSet;
typedef set<int>::iterator TIntSetIter;
typedef set<double> TDblSet;
typedef map<int, int> TIntIntMap;
typedef queue<int> TIntQueqe;
typedef stack<int> TIntStack;

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
  TIntSetV GetNeighbor() { return Neighbor; }
  // get number of nodes
  int GetNumNodes() { return NumNodes; }
  // get number of edges
  int GetNumEdges() { return NumEdges; }
  // get degree of a node
  int GetDeg(int NID) { return Neighbor[NID-1].size(); }
  // get degree of all nodes
  TIntV GetDegV() { 
    TIntV DegAllNodes;
    for (int NID = 1; NID <= GetNumNodes(); NID++) { DegAllNodes.push_back(GetDeg(NID)); } 
    return DegAllNodes;
  }
  
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
  
  // test the functionality of this Graph class
  void FuncTest();
  
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
}


//////////////////////////////////
///////////////////////////////////////////////
// class of evaluating a graph's basic characteristics
// constructed in in 2021-01-16 by Yuyao Wang
class EVALGRAPH {
public:
  TIntV GetMaxDeg(GRAPH G); // return index from 1 and maxvalue 
  TIntV GetMinDeg(GRAPH G); // return index from 1 and minvalue 
  int GetAvgDeg(GRAPH G); // return average degree of a graph
  TIntIntMap GetDegreeDistribution(GRAPH G); // key is number of degree value, value is number of nodes
  bool GetConnectivityBFS(TIntSetV Neighbor, int TargetNode); // BFS for determining the connectivity of the graph
  bool GetConnectivityDFS(TIntSetV Neighbor, int TargetNode); // DFS for determining the connectivity of the graph
  
  void FuncTest(GRAPH G); // for function test
}; 

TIntV EVALGRAPH::GetMaxDeg(GRAPH G) {
  TIntV DegAllNodes = G.GetDegV();
  vector<int>::iterator MaxPos = max_element(DegAllNodes.begin(), DegAllNodes.end());
  TIntV MaxPosValue;
  MaxPosValue.push_back(distance(DegAllNodes.begin(), MaxPos) + 1); // position index, begin from 1
  MaxPosValue.push_back(*MaxPos); // value
  cout << "max degree, node-id: " << MaxPosValue[0] << ", value= " << MaxPosValue[1] << endl;
  return MaxPosValue;
}

TIntV EVALGRAPH::GetMinDeg(GRAPH G) {
  TIntV DegAllNodes = G.GetDegV();
  vector<int>::iterator MinPos = min_element(DegAllNodes.begin(), DegAllNodes.end());
  TIntV MinPosValue;
  MinPosValue.push_back(distance(DegAllNodes.begin(), MinPos) + 1); // position index, begin from 1
  MinPosValue.push_back(*MinPos); // value
  cout << "min degree, node-id: " << MinPosValue[0] << ", value= " << MinPosValue[1] << endl;
  return MinPosValue;
}

int EVALGRAPH::GetAvgDeg(GRAPH G) {
  int SumDeg = 0;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) { SumDeg += G.GetDeg(NID); }
  SumDeg /= G.GetNumNodes();
  cout << "average degree: " << SumDeg << endl;
  return SumDeg;
}

TIntIntMap EVALGRAPH::GetDegreeDistribution(GRAPH G) {
  TIntIntMap DegDistirb;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) {
    int NDeg = G.GetDeg(NID);
    //cout << "Node-id: " << NID << ", Degree: " << NDeg << endl; 
    TIntIntMap::iterator iter = DegDistirb.find(NDeg);
    if (iter != DegDistirb.end()) {
      DegDistirb[NDeg] = iter->second + 1; // if the key exists, corresponding value plus 1
    } else {
      DegDistirb.insert(pair<int, int>(NDeg, 1)); // add a new key if the key does not exist, set the value as 1
    }
  }
  /*cout << endl;
  for (TIntIntMap::iterator iter = DegDistirb.begin(); iter != DegDistirb.end(); iter++) { 
    cout << "Degree: " << iter->first << ", Number: " << iter->second << endl; 
  }*/
  return DegDistirb;
}

bool EVALGRAPH::GetConnectivityBFS(TIntSetV Neighbor, int TargetNode) {
  // BFS
  TBoolV Flag(Neighbor.size(), false);
  queue<int> Q1;
  Q1.push(TargetNode);
  while (!Q1.empty()) {
    int PushNode = Q1.front();
    Q1.pop();
    if (Flag[PushNode-1] == 0) { 
      Flag[PushNode-1] = 1;
      TIntSet TempNgh = Neighbor[PushNode-1];
      for (TIntSetIter iter = TempNgh.begin(); iter != TempNgh.end(); ++iter) {
        Q1.push(*iter);
      }
    }
  }
  int NumIterated = 0;
  for (int j = 0; j < Flag.size(); j++) { NumIterated += Flag[j]; }
  if (NumIterated==Flag.size()) {
    cout << "This is a connected graph!!! (BFS)" << endl;
    return 1;
  } else { 
    cout << "This is not a connected graph!!! (BFS)" << endl;
    return 0;
  }
}

bool EVALGRAPH::GetConnectivityDFS(TIntSetV Neighbor, int TargetNode) {
  // DFS
  TBoolV Flag(Neighbor.size(), false);
  stack<int> S1;
  S1.push(TargetNode);
  while (!S1.empty()) {
    int PushNode = S1.top();
    S1.pop();
    if (Flag[PushNode-1] == 0) { 
      Flag[PushNode-1] = 1;
      TIntSet TempNgh = Neighbor[PushNode-1];
      for (TIntSetIter iter = TempNgh.begin(); iter != TempNgh.end(); ++iter) {
        S1.push(*iter);
      }
    }
  }
  int NumIterated = 0;
  for (int j = 0; j < Flag.size(); j++) { NumIterated += Flag[j]; }
  if (NumIterated==Flag.size()) {
    cout << "This is a connected graph!!! (DFS)" << endl;
    return 1;
  } else {
    cout << "This is not a connected graph!!! (DFS)" << endl;
    return 0; 
  }
}

void EVALGRAPH::FuncTest(GRAPH G) {
  cout << "===========graph statistics==============" << endl;
  GetMaxDeg(G);
  GetMinDeg(G);
  GetAvgDeg(G);
  GetConnectivityBFS(G.GetNeighbor(), 1);
  GetConnectivityDFS(G.GetNeighbor(), 1);
  cout << "=========================================" << endl;
}


//////////////////////////////////
///////////////////////////////////////////////
// class of evaluating nodes' centrality
// constructed in in 2021-01-14 by Yuyao Wang
class CENTRALITY {
public:
  // conductance centrality of one node
  double CondOneNode(GRAPH G, int NID);
  // conductance centrality of all nodes in a graph
  TDblV CondAllNodes(GRAPH G);
};

double CENTRALITY::CondOneNode(GRAPH G, int NID) {
// conductance centrality (Gleich, et al. KDD12)
// last modified in 2021-01-14 by Yuyao Wang
  double ImptNode = 0.0; // small value indicates an important node.
  TIntSetV Neighbor = G.GetNeighbor();
  TIntSet NghNId(Neighbor[NID - 1]);
  NghNId.insert(NID); // did in SNAP
  int LenNghs = NghNId.size();
  if (LenNghs < 5) {
    ImptNode = 1.0;
    return ImptNode;
  }
  int Edges2 = 2 * G.GetNumEdges();
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

TDblV CENTRALITY::CondAllNodes(GRAPH G) {
// conductance centrality of all nodes in a graph
// small value indicates an important node.
  TDblV Conductance;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) {
    Conductance.push_back(CondOneNode(G, NID));
  }
  return Conductance;
}


//////////////////////////////////
///////////////////////////////////////////////
// class of file output
// constructed in 2021-01-14 by Yuyao Wang
class IOCLASS {
public:
  // input a TwoDimVV (stored in TIntVV).
  void LoadIntVV(string DataName, string Suffix, TIntVV& TwoDimVV);
  // input a TwoDimVV (stored in TDblVV).
  void InputDblVV(string DataName, string Suffix, TDblVV& TwoDimVV);
  
public:
  // output a vector (stored in T[X]V).
  template <class T> void OutputVec(string FileName, string Suffix, vector<T> OneDimVec);
  // output a set (stored in T[X]Set).
  template <class T> void OutputSet(string FileName, string Suffix, set<T> OneDimSet);
  // output a TwoDimVV (stored in T[X]VV).
  template <class T> void OutputVV(string FileName, string Suffix, vector<vector<T> > TwoDimVV);
  // output a TwoDimSetV (stored in T[X]SetV).
  template <class T> void OutputSetV(string FileName, string Suffix, vector<set<T> > TwoDimSetV);
  // output a map (stored in T[X][X]Map)
  template <class T> void OutputMap(string FileName, string Suffix, map<T, T> FMap);
  // output a .gml file which can be fed into Cytoscape software for visualizing a graph.
  void OutputGML(string DataName, TIntSetV Neighbor, TIntSetV ClusInNode);
};

void IOCLASS::LoadIntVV(string DataName, string Suffix, TIntVV& TwoDimVV) {
  string FileName = DataName + "." + Suffix;
  const char* CharFileName = FileName.c_str();
  if (access(CharFileName, R_OK|W_OK) != 0) {
    printf("No IntVV file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finIntVV;
  finIntVV.open(CharFileName);
  string szLine;
  TIntV TempV;
  while(getline(finIntVV, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      TempV.push_back(atoi((*iter).c_str()));
    }
    TwoDimVV.push_back(TempV);
    TempV.clear();  
  }
  finIntVV.close();
  finIntVV.clear();
}

void IOCLASS::InputDblVV(string DataName, string Suffix, TDblVV& TwoDimVV) {
  string FileName = DataName + "." + Suffix;
  const char* CharFileName = FileName.c_str();
  if (access(CharFileName, R_OK|W_OK) != 0) {
    printf("No DblVV file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finDblVV;
  finDblVV.open(CharFileName);
  string szLine;
  TDblV TempV;
  while(getline(finDblVV, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      TempV.push_back(atof((*iter).c_str()));
    }
    TwoDimVV.push_back(TempV);
    TempV.clear();  
  }
  finDblVV.close();
  finDblVV.clear();
}

template <class T> 
void IOCLASS::OutputVec(string FileName, string Suffix, vector<T> OneDimVec) {
// output a vector.
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutV;
  foutV.open(File);
  for (int ind1 = 0; ind1 < OneDimVec.size(); ind1++) { foutV << OneDimVec[ind1] << endl; }
  foutV.close();
  foutV.clear();
}

template <class T> 
void IOCLASS::OutputSet(string FileName, string Suffix, set<T> OneDimSet) {
// output a vector.
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutSet;
  foutSet.open(File);
  typename std::set<T>::iterator TempIter;
  for (TempIter = OneDimSet.begin(); TempIter != OneDimSet.end(); TempIter++) {
    foutSet << *TempIter << " ";
  }
  foutSet.close();
  foutSet.clear();
}

template <class T> 
void IOCLASS::OutputVV(string FileName, string Suffix, vector<vector<T> > TwoDimVV) {
// output a TwoDimVV (stored in T[X]VV)
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutVV;
  foutVV.open(File);
  for (int ind1 = 0; ind1 < TwoDimVV.size(); ind1++) {
    for (int ind2 = 0; ind2 < TwoDimVV[ind1].size(); ind2++) { foutVV << TwoDimVV[ind1][ind2] << " "; }
    foutVV << endl;
  }
  foutVV.close();
  foutVV.clear();
}

template <class T> 
void IOCLASS::OutputSetV(string FileName, string Suffix, vector<set<T> > TwoDimSetV) {
// output a TwoDimSetV (stored in T[X]SetV)
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutSetV;
  foutSetV.open(File);
  typename std::set<T>::iterator TempIter;
  for (int ind1 = 0; ind1 < TwoDimSetV.size(); ind1++) {
    for (TempIter = TwoDimSetV[ind1].begin(); TempIter != TwoDimSetV[ind1].end(); TempIter++) {
      foutSetV << *TempIter << " ";
    }
  foutSetV << endl;
  }
  foutSetV.close();
  foutSetV.clear();
}

template <class T>
void IOCLASS::OutputMap(string FileName, string Suffix, map<T, T> FMap) {
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutMap;
  foutMap.open(File);
  typename std::map<T, T>::iterator TempIter;
  for (TempIter = FMap.begin(); TempIter != FMap.end(); TempIter++) { 
    foutMap << TempIter->first << " " << TempIter->second << endl; 
  }
  foutMap.close();
  foutMap.clear();
}

void IOCLASS::OutputGML(string FileName, TIntSetV Neighbor, TIntSetV ClusInNode) {
// visualize the graph (output a .gml file which can be fed into software named Cytoscape
  string SGmlFile = FileName + ".gml";
  const char* GmlFile = SGmlFile.c_str();
  ofstream foutG;
  foutG.open(GmlFile);
  foutG << "graph" << endl;
  foutG << "[" << endl;
  foutG << "  directed 0" << endl;
  TIntSetIter it_set;
  int i;
  for (i = 0; i < Neighbor.size(); i++) {
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
  
  for (i = 0; i < Neighbor.size(); i++) {
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
  foutG.clear();
}


//////////////////////////////////
///////////////////////////////////////////////
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
  TIntSet inline CalJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    return JointSet;
  }
};

double ClusterTest::CalNonOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
// last modified in 2021-01-13 by Yuyao Wang.
// calculate modularity of non-overlapping cluster structure. (Newman, PRE, 2004.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// if input is a overlapping cluster structure, this function will only take the first cluster for each node as her FuzzyMembership.
// num_edges is the number of all edges in the graph.
// the value of modularity lies between -1 and 1.
  if (neighbor.size() != Our2.size()) { 
    cout << "inappropriate inputs!!!" << endl;
    exit(1); 
  }
  TDblV theta;
  theta.resize(Our1.size());
  TIntSetIter it_set;
  for (int CID = 0; CID < theta.size(); CID++) {
    theta[CID] = 0.0;
    // total degree (in-degree + out-degree) of kth community, for the pre-computation of \sum_j d_j in the equation
    // here we sum all nodes including the target node who is removed sometimes in some papers.
    for (it_set = Our1[CID].begin(); it_set != Our1[CID].end(); it_set++) {
      theta[CID] += neighbor[*it_set - 1].size();
    }
    theta[CID] /= 2.0 * num_edges;
  }
  // modularity of single node
  TDblV singleModul;
  singleModul.resize(neighbor.size());
#pragma omp parallel for
  for (int NID = 1; NID <= neighbor.size(); NID++) {
    double tempValue = 0.0;
    singleModul[NID-1] = 0.0;
    // determine whether NID does not belong to any cluster.
    if (Our2[NID-1].empty()) { 
      singleModul[NID-1] = 0.0; 
    } else if (neighbor[NID-1].empty() == 0) {
      it_set = Our2[NID-1].begin(); // adopt the first value for non-overlapping communities (if there are multiple values)
      tempValue = CalNumJointSets(neighbor[NID-1], Our1[*it_set - 1]); // first term of equation
      singleModul[NID-1] = (tempValue - theta[*it_set - 1] * neighbor[NID-1].size()) / 2.0 / num_edges;
    }    
  }
  double SumModul = 0.0;
  // total value of modularity
  for (int NID = 1; NID <= neighbor.size(); NID++) { SumModul += singleModul[NID-1]; }
  return SumModul;
}

double ClusterTest::CalOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating modularity of overlapping and hierarchical cluster structure.(Huawei Shen, et al. Eqn.(2) in Physica A, 2009.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// num_edges is the number of all edges in the graph.
// here we tranform the defintion of community-oriented to node-oriented.
// In particular, we divide 1/O_iO_j.
  if (neighbor.size() != Our2.size()) { 
    cout << "inappropriate inputs!!!" << endl;
    exit(1); 
  }
  // modularity of single node
  TDblV singleModul;
  singleModul.resize(neighbor.size());
//#pragma omp parallel for
  for (int NID = 1; NID <= neighbor.size(); NID++) {
    double tempValue1 = 0.0, tempValue2 = 0.0;
    singleModul[NID-1] = 0.0;
    // determine whether NID does not belong to any cluster.
    if (Our2[NID-1].empty()) { 
      singleModul[NID-1] = 0.0; 
    } else if (neighbor[NID-1].empty() == 0) {
      for (TIntSetIter it_set = neighbor[NID-1].begin(); it_set != neighbor[NID-1].end(); it_set++) {
        if ((Our2[*it_set - 1].empty() == 0) && (CalNumJointSets(Our2[NID-1], Our2[*it_set - 1]) > 0)) {
          tempValue1 += (1.0 / Our2[*it_set - 1].size());
          tempValue2 += (neighbor[NID-1].size() * neighbor[*it_set - 1].size() / 2.0 / num_edges / Our2[*it_set - 1].size());
        }
      }
      singleModul[NID-1] = (tempValue1 - tempValue2) / Our2[NID-1].size();
    }    
  }
  double SumModul = 0.0;
  // total value of modularity
  for (int NID = 1; NID <= neighbor.size(); NID++) { SumModul += singleModul[NID-1]; }
  SumModul /= (2.0 * num_edges);
  return SumModul;
}

double ClusterTest::CalTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating tightness of cluster structure (Zhan Bu, et al. Information Fusion, 2017)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
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
#pragma omp parallel for
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      //id_i = *it_set;
      int TempNumJointSets = CalNumJointSets(neighbor[*it_set - 1], Our1[k]);
      numInEdges[k] += TempNumJointSets;
      numOutEdges[k] += (neighbor[*it_set - 1].size() - TempNumJointSets);
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
  }
  for (int k = 0; k < num_clusters; k++) {
    if (Our1[k].empty() == 0 && Our1[k].size() != num_nodes) {
      FirTerm = 2.0 * numInEdges[k] / Our1[k].size() / (double)Our1[k].size();
      SecTerm = numOutEdges[k] / Our1[k].size() / (double)(num_nodes - Our1[k].size());
      SinglTgt[k] = FirTerm - SecTerm;
    }
    SumTgt += SinglTgt[k];
  }
  return SumTgt;
}

double ClusterTest::CalAdjTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating adjusted tightness of cluster structure (Zhan Bu, et al. Information Fusion, 2017)
// penalizing very small and very large clusters and produces well-balanced solutions
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
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
#pragma omp parallel for
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      //id_i = *it_set;
      int TempNumJointSets = CalNumJointSets(neighbor[*it_set - 1], Our1[k]);
      numInEdges[k] += TempNumJointSets;
      numOutEdges[k] += (neighbor[*it_set - 1].size() - TempNumJointSets);
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
  }
  for (int k = 0; k < num_clusters; k++) {
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
  CENTRALITY myCentrality;
  for (int NId = 1; NId <= F.size(); NId++) {
	double Phi = myCentrality.CondOneNode(G, NId);
    NIdPhiV.push_back(make_pair(Phi, NId));
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


//////////////////////////////////
///////////////////////////////////////////////
class EGTCD {
private:
  GRAPH G;
  int NumComs; // number of communities
  //TIntSetV Neighbor1; // topology of first tier graph
  
  int TypeSecOrdProximity; // type of calculation of second-order proximity
  double ThresGenSecOrdGraph; // [0,1], threshold for generating second order graph.
  
  TIntSetV Neighbor2; // topology of second-order graph
  int NumEdges2; // number of edges of second-order graph
  TIntVV NumJointNgh; // store the number of joint neighbors
  
public:
  EGTCD(GRAPH graph): G(graph.GetDataName()), TypeSecOrdProximity(1), ThresGenSecOrdGraph(0.0) {}
  void Initialization(); // load data
  double NghBsdProx(int NID1, int NID2, int type);
  void GenerateSecondOrderGraph(TIntSetV& Neighbor2, int TypeSecOrdProximity, double ThresGenSecOrdGraph);
  void PlotMotivation(); // PlotMotivation
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
};

void EGTCD::Initialization() {
  G.LoadGraphTpl(); // load .graph file (pure directed graph)
  printf("Number of Nodes: %d, Number of edges: %d\n", G.GetNumNodes(), G.GetNumEdges());
  G.LoadClusGt(); // load .gt file (cluster ground-truth)
  printf("Number of clusters: %d\n", G.GetNumClus());
}

// neighbor based proximity
double EGTCD::NghBsdProx(int NID1, int NID2, int type) {
  TIntSet FirNode = G.GetNeighbor(NID1);
  TIntSet SecNode = G.GetNeighbor(NID2);
  int LenFirNode = FirNode.size();
  int LenSecNode = SecNode.size();
  if (LenFirNode == 0 || LenSecNode == 0) { return 0.0; }
  int NumJoint = NumJointNgh[NID1-1][NID2-NID1-1]; // pre-computed
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

void EGTCD::GenerateSecondOrderGraph(TIntSetV& Neighbor2, int TypeSecOrdProximity, double ThresGenSecOrdGraph) {
  Neighbor2.resize(G.GetNumNodes()); // core dumped if the space is not pre-allocated.
  NumEdges2 = 0;  
  // allocate space for SecOrdProx
  TDblVV SecOrdProx(G.GetNumNodes());
  for(int NID = 1; NID <= SecOrdProx.size(); NID++) { SecOrdProx[NID-1].resize(SecOrdProx.size()-NID); }
  
#pragma omp parallel for
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      SecOrdProx[NID1-1][NID2-NID1-1] = NghBsdProx(NID1, NID2, TypeSecOrdProximity);
    }
  }
  
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      if(SecOrdProx[NID1-1][NID2-NID1-1] > ThresGenSecOrdGraph) {
        // construct a link when two nodes share at least one neighbor.
        Neighbor2[NID1-1].insert(NID2); // this opearation can not be parallized
        Neighbor2[NID2-1].insert(NID1);
      }
    }
  }
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) { NumEdges2 += Neighbor2[NID1-1].size(); }
  NumEdges2 = NumEdges2 / 2;
  // cout << "Number of Edges in Second-Tier Graph: " << NumEdges2 << endl;
}

void EGTCD::PlotMotivation() {
  // for motivation plot
  Initialization();
  
  // allocate space for NumJointNgh.
  NumJointNgh.resize(G.GetNumNodes());
  for(int NID = 1; NID <= NumJointNgh.size(); NID++) { NumJointNgh[NID-1].resize(NumJointNgh.size()-NID); }
  // calculate NumJointNgh
#pragma omp parallel for
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    //cout << "for parallization test: " << NID1 << endl;
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      NumJointNgh[NID1-1][NID2-NID1-1] = CalNumJointSets(G.GetNeighbor(NID1), G.GetNeighbor(NID2));
    }
  }
  cout << "Calculation of Number of Joint Neighbors is fininshed!!!" << endl;
  
  TDblVV StoreMetric;
  StoreMetric.resize(1000);
  ClusterTest mytest; 
  // Modularity-1st-order (no overlapping)
  double Modul1 = mytest.CalNonOverlapModul(G.GetClusInClus(), G.GetClusInNode(), G.GetNeighbor(), G.GetNumEdges());
  // Tightness-1st-order
  double AdjTgt1 = mytest.CalAdjTgt(G.GetClusInClus(), G.GetClusInNode(), G.GetNeighbor());
  cout << "Modul1= " << Modul1 << endl;
  cout << "AdjTgt1= " << AdjTgt1 << endl;
  StoreMetric[0].push_back(Modul1);
  StoreMetric[0].push_back(AdjTgt1);
  
  for (TypeSecOrdProximity = 1; TypeSecOrdProximity <= 6; TypeSecOrdProximity++) {
    if (TypeSecOrdProximity == 1) { printf("Second-order proximity: Jaccard Coefficient.\n"); }
      else if (TypeSecOrdProximity == 2) { printf("Second-order proximity: Salton Index.\n"); }
      else if (TypeSecOrdProximity == 3) { printf("Second-order proximity: Sorensen Index.\n");  }
      else if (TypeSecOrdProximity == 4) { printf("Second-order proximity: Hub Promoted Index.\n"); }
      else if (TypeSecOrdProximity == 5) { printf("Second-order proximity: Hub Depressed Index.\n"); }
      else if (TypeSecOrdProximity == 6) { printf("Second-order proximity: Leicht-Holme-Newman Index.\n"); }
    for (int IndThresGenSecOrdGraph = 0; IndThresGenSecOrdGraph < 100; IndThresGenSecOrdGraph++) {
      ThresGenSecOrdGraph = IndThresGenSecOrdGraph / 100.0;
      TIntSetV TempNeighbor2;
      GenerateSecondOrderGraph(TempNeighbor2, TypeSecOrdProximity, ThresGenSecOrdGraph);
      Neighbor2 = TempNeighbor2;
      int index = (TypeSecOrdProximity-1)*100 + IndThresGenSecOrdGraph;
      // Modularity-2nd-order (no overlapping)
      StoreMetric[index+1].push_back(mytest.CalNonOverlapModul(G.GetClusInClus(), G.GetClusInNode(), TempNeighbor2, NumEdges2));
      // Tightness-2nd-order
      StoreMetric[index+1].push_back(mytest.CalAdjTgt(G.GetClusInClus(), G.GetClusInNode(), TempNeighbor2));
      //cout << index+1 << ": " << StoreMetric[index][0] << " " << StoreMetric[index][1]<< endl;
      cout << StoreMetric[index][0] << " ";
    }
    cout << endl;
  }
  
  IOCLASS myOutput;
  myOutput.OutputVV(G.GetDataName(), "MetricResults", StoreMetric);
  string DataName = G.GetDataName() + "2";
  myOutput.OutputGML(DataName, Neighbor2, G.GetClusInNode());
}


////////////////////
//////////////////////////////
class FUNCTIONTEST {
private:
  GRAPH G;
public:
  FUNCTIONTEST(GRAPH graph): G(graph.GetDataName()) {}
  void PerformTest();
};

void FUNCTIONTEST::PerformTest() {
  G.LoadGraphTpl(); // load .graph file (pure directed graph)
  printf("Number of Nodes: %d, Number of edges: %d\n", G.GetNumNodes(), G.GetNumEdges());
  EVALGRAPH EvalGraph;
  EvalGraph.FuncTest(G);
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
  gettimeofday(&tod1, NULL);
  FUNCTIONTEST FuncTest(myGraph);
  FuncTest.PerformTest();
  gettimeofday(&tod2, NULL);
  cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
  cout << endl;

  return 0;
}
