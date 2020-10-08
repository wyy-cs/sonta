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
#define FALSE 0

using namespace std;

// generate a random vector whose summation of all elements equals 1.
double* randArrSumOne(int NumEntry) {
	//NumEntry is the number of elements
	srand((unsigned)time(0) * rand());
	double* Array = new double[NumEntry + 1];
	Array[0] = NumEntry;
	double SumEntry = 0.0;
	int i, TempRand;
	for (i = 1; i <= NumEntry; i++) {
		TempRand = rand();// randomly generate a random value 0 - 0x7FFF(i.e., 0 -- RAND_MAX)
		Array[i] = (double)TempRand / RAND_MAX;//map the random value to the internal (0, 1)
		SumEntry += Array[i];
	}
	for (i = 1; i <= NumEntry; i++) { Array[i] /= SumEntry; }
	return Array;
}

// record the cost time from tod1 to tod2
long long int todiff(struct timeval* tod1, struct timeval* tod2) {
    long long t1, t2;
    t1 = tod1->tv_sec * 1000000 + tod1->tv_usec;
    t2 = tod2->tv_sec * 1000000 + tod2->tv_usec;
    return t1 - t2;
}

// concatenate two string (in char*) into one string (in char*)
char* strcatre(char* s1, char* s2) {
	int i, len1, len2;
	len1 = strlen(s1);
	len2 = strlen(s2);
	char* s3 = new char[len1 + len2 + 1];
	while(*s1 != '\0') { *(s3++) = *(s1++); }
	while(*s2 != '\0') { *(s3++) = *(s2++); }
	*(s3++) = '\0';
	for(i = 0; i < len1 + len2 + 1; i++) { s3--; }
	return s3;
}

// determine whether a value exists in an array. 
// If it exists, return the index of the first element equal to this value;
// return 0, otherwise.
template<typename Type>
int find(Type* set1, Type object) {
	if (set1[0] == 0) { return 0; }
	int i = 0;
	int ind = 0;
	for (i = 1; i <= set1[0]; i++) {
		if (object == set1[i]) { ind = i; break; }
	}
	return ind;
}

// size of intersection set between two sets
// the elements of each set must be stored in ascending order!!!!
template <typename Type>
int NumIntersection(const Type* set1, const Type* set2) {
	int size_1 = set1[0];
	int size_2 = set2[0];
	int i = 1;
	int j = 1;
	int joint_number = 0;
	while(1) {
		if (i > set1[0] || j > set2[0]) { break; }
		else if (set1[i] > set2[j]) { j++; }
		else if (set1[i] < set2[j]) { i++; }
		else {
			i++;
			j++;
			joint_number++;
		}
	}
	return joint_number;
}

// size of union set between two sets
// the elements of each set must be stored in ascending order!!!!
template <typename Type>
int NumUnion(const Type* set1, const Type* set2) {
	int size_1 = set1[0];
	int size_2 = set2[0];
	int i = 1;
	int j = 1;
	int joint_number = 0;
	while(1) {
		if (i > set1[0] || j > set2[0]) { break; }
		else if (set1[i] > set2[j]) { j++; }
		else if (set1[i] < set2[j]) { i++; }
		else {
			i++;
			j++;
			joint_number++;
		}
	}
	int union_number = size_1 + size_2 - joint_number;
	return union_number;
}

// basic graph class
class Graph {
public:
	// dataset name (do not contain the suffix name)
	char* name_dataset;
	// fixed for calculating second-proximity
	int TypeTopolSimil;
	
	/*----------------------------basic topology information of pure graph--------------------------*/
	// number of nodes
	int num_nodes;
	// number of edges
	int num_edges;
	// network topology, each one-dimension array of corresponding node contains its all neighbors,
	// and the first position (for example, neighbor[i][0] for node i) stores the number of its neighbors
	int** neighbor;
	// second-order proximity
	double** proximity;

	/*----------------------------real cluster information (ground truth)----------------------------*/
	// each one-dimension array of corresponding node contains the cluster ids she belongs,
	// and the first position (for example, cluster[i][0] for node i) stores the number of clusters she belongs
	int** cluster;
	// each one-dimension array of corresponding cluster contains all nodes in this cluster,
	// and the first position (for example, clusters[k][0] for cluster k) stores all nodes contained in this cluster
	int** clusters;
	// number of real clusters
	int num_clusters; 
	// number of common communities of node pairs
	int** NumCommCommu;
	// the correlation between the number of common clusters and average of proximity.
	// first value correlation[0] is not the number of elements, 
	// here is the correlation between (common clusters equals 0) and Proximity.
	double* correlation;
	
	/*------------------weighted graph (have not used/defined, 2020-10-4)------------*/
	// values on nodes (each node is associated with a single value (e.g., centrality))
	double* wgt_nodes;
	// values on nodes (each node is associated with a vector)
	double** fea_nodes;
	// length of feature of nodes when each node is associated with a vector
	int len_fea_nodes;
	
	// values on edges (each edge is associated with a single value (e.g., weight))
	double** wgt_edges;
	// values on edges (each edge is associated with a vector)
	double*** fea_edges;
	// length of feature of edges when each edge is associated with a vector
	int len_fea_edges;
	
public:
	// constructor (initialization, file name is needed)
	Graph(char* dataname) { name_dataset = dataname; }
	// load a graph with topology information
	void ReadPureDirGraph();
	// load a graph with topology and node feature information
	void ReadNodeFeature();
	// calculate second-order proximity of one node pair
	double NeighborBasedSimilarity(int* Node_1, int* Node_2, int type);
	// calculate second-order proximity of all node pairs in a graph
	void CalSecOrdProx(int type);
	
	// load real cluster information
	void ReadGTclusters();
	// calculate the number of common communities of node pairs
	void CalculateNumCommCommu();
	// calculate the correlation
	void CalCorrelation();
	// output the correlation
	void OutputCorrelation();
	
	// visualize the real cluster strucutre
	void VisualizeGraph();
};

void Graph::ReadPureDirGraph() {
	ifstream fin1;
	int i = 0;
	int j = 0;
	int n = 0;
	int m = 0;
	string szLine;

	char* graphtail = (char*)".graph";
	char* graphfile = strcatre(name_dataset,graphtail);
	if (access(graphfile, R_OK|W_OK) != 0) {
		cout << "There is no graph topology (.graph) file, please check it!!!" << endl;
		exit(0);
	}
	
	fin1.open(graphfile);
	while(getline(fin1,szLine)) { i++; }
	fin1.close();
	fin1.clear();
	num_nodes = i;
	
	neighbor = new int*[num_nodes + 1];
	neighbor[0] = new int[1];
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	neighbor[0][0] = num_nodes;
	
	num_edges = 0;
	i = 1;
	
	fin1.open(graphfile);
	while(getline (fin1,szLine)) {
		vector<string> tData;
		istringstream iss(szLine);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tData));
		neighbor[i] = new int[tData.size() + 1];
		neighbor[i][0] = tData.size();
		j = 1;
		set<int> tempset;
		set<int>::iterator it_set;
		// rewrite into set<int>, sort the neighbor set, 
		// avoiding bug in the process of calculating union/intersection set
		for(vector<string>::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
			tempset.insert(atoi((*iter).c_str()));
		}
		for (it_set = tempset.begin(); it_set != tempset.end(); it_set++) {
			num_edges++;
			neighbor[i][j] = *it_set;	
			j++;
		}
		tempset.clear();
		i++;	
	}
	fin1.close();
	fin1.clear();
	// each edge is twicely recorded (undirected graph)
	num_edges = num_edges / 2;
	cout << "Number of Nodes: " << num_nodes << endl;
	cout << "Number of Edges: " << num_edges << endl;
}

void Graph::ReadNodeFeature() {
	ifstream fin2;
	int i = 0;
	int j = 0;
	int n = 0;
	int m = 0;
	int TGraph = 0;
	string szLine;
	
	char* featuretail = (char*)".nf";
	char* featurefile = strcatre(name_dataset, featuretail);
	if (access(featurefile, R_OK|W_OK) != 0) {
		cout << "There is no node feature (.nf) file!!!" << endl;
		return;
	}
	
	// if there exists a .fea file
	int maxFeaId = 0;
	fin2.open(featurefile);
	while(getline(fin2,szLine)) {
		vector<string> tData;
		istringstream iss(szLine);
		copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<vector<string> >(tData));
		// find the maximum feature id
		for(vector<string>::iterator iter = tData.begin(); iter!=tData.end(); ++iter) {
			int tempd = atoi((*iter).c_str());
			if (maxFeaId <= tempd) { maxFeaId = tempd; }
		}	
	}
	fin2.close();
	fin2.clear();
	// length of each node's features
	len_fea_nodes = maxFeaId;
	cout << "Length of node features: " << len_fea_nodes << endl;
	
	fea_nodes = new double*[num_nodes + 1];
	fea_nodes[0] = new double[1];
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	fea_nodes[0][0] = num_nodes;
	for (i = 1; i <= num_nodes; i++) {
		fea_nodes[i] = new double[len_fea_nodes + 1];
		fea_nodes[i][0] = len_fea_nodes;
		for (j = 1; j <= len_fea_nodes; j++) { fea_nodes[i][j] = 0.0; }
	}
	
	i = 1;
	fin2.open(featurefile);
	while(getline(fin2,szLine)) {
		vector<string> tData;
		istringstream iss(szLine);
		copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<vector<string> >(tData));
		// reopen and write to fea_nodes
		// if this node has a feature, then the value of corresponding position is 1, 0 otherwise.
		for(vector<string>::iterator iter = tData.begin(); iter != tData.end(); ++iter) {
			fea_nodes[i][atoi((*iter).c_str())] = 1;
		}
		i++;	
	}
	fin2.close();
	fin2.clear();
}

double Graph::NeighborBasedSimilarity(int* Node_1, int* Node_2, int type) {
	double joint_number = 0.0;
	joint_number = (double)NumIntersection(Node_1, Node_2);
	double size_a = Node_1[0];
	double size_b = Node_2[0];
	double union_number = size_a + size_b - joint_number;
	joint_number += 2;
	union_number += 2;
	size_a++;
	size_b++;
	// Jaccard Coefficient
	if (type == 1) { return joint_number / union_number; }
	// Salton Index
	else if (type == 2) { return joint_number / (sqrt(size_a * size_b)); }
	// Sorensen Index
	else if (type == 3) { return 2 * joint_number / (size_a + size_b); }
	// Hub Promoted Index
	else if (type == 4) {
		if (size_a < size_b) { return 2 * joint_number / size_a; }
		else { return 2 * joint_number / size_b; }
	}
	// Hub Depressed Index
	else if (type == 5) {
		if (size_a < size_b) { return 2 * joint_number / size_b; }
		else { return 2 * joint_number / size_a; }
	}
	// Leicht-Holme-Newman Index
	else if (type == 6) { return 2 * joint_number / size_a / size_b; }
	return 0.0;
}

// calculate the second-order proximity
void Graph::CalSecOrdProx(int type) {
	TypeTopolSimil = type; // initialize class FuzzyMember variable in class Graph
	if (type == 1) { cout << "Second-order proximity: Jaccard Coefficient" << endl; }
	else if (type == 2) { cout << "Second-order proximity: Salton Index" << endl; }
	else if (type == 3) { cout << "Second-order proximity: Sorensen Index" << endl;	}
	else if (type == 4) { cout << "Second-order pairroximity: Hub Promoted Index" << endl; }
	else if (type == 5) { cout << "Second-order pairroximity: Hub Depressed Index" << endl; }
	else if (type == 6) { cout << "Second-order pairroximity: Leicht-Holme-Newman Index" << endl; }
	int i = 1;
	int j = 1;
	int k = 1;
	proximity = new double*[num_nodes + 1];
	proximity[0] = new double[1];
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	proximity[0][0] = num_nodes;
	
	for (i = 1; i <= num_nodes; i++) {
		proximity[i] = new double[num_nodes - i + 1];
		proximity[i][0] = num_nodes - i;
		for (j = i+1; j <= num_nodes; j++) {
			// double Graph::NeighborBasedSimilarity(int* Node_1, int* Node_2, int type)
			proximity[i][j-i] = NeighborBasedSimilarity(neighbor[i], neighbor[j], TypeTopolSimil);
		}
		printf("\rProximity Calculation Completed Progress: %.2lf%%(%d)", i / double(num_nodes) * 100, i);
		fflush(stdout);
	}
	cout << endl;
}

// load real cluster information
void Graph::ReadGTclusters() {
	ifstream fin;
	
	char* clustertail = (char*) ".gt";
	char* clusterfile = strcatre(name_dataset, clustertail);
	
	string szLine;
	int i = 0;
	int j = 0;
	int temp = 0;
	int q = 0;
	set<int>::iterator it_set;
	
	if (access(clusterfile, R_OK|W_OK) != 0) {
		cout << "There is no cluster (.gt) file!!!" << endl;
		return;
	}
	fin.open(clusterfile);
	
	while(getline(fin,szLine)) { i++; }
	fin.close();
	fin.clear();
	
	num_clusters = i;
	cout << "Number of clusters: " << num_clusters << endl;
	
	clusters = new int*[num_clusters + 1];
	clusters[0] = new int[1];
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	clusters[0][0] = num_clusters;
	
	// store the cluster information, each array of corresponding node has num_clusters elements,
	// the corresponding value is 1 if this node belongs to the corresponding cluster
	int** state = new int*[num_nodes + 1];
	state[0] = new int[1]; 
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	state[0][0] = num_nodes;
	
	for (i = 1; i <= num_nodes; i++) {
		state[i] = new int[num_clusters + 1];
		for (j = 0; j <= num_clusters; j++) { state[i][j] = 0; }
	}
	
	fin.open(clusterfile);
	int p = 1;
	while(getline (fin,szLine)) {
		vector<string> tData;
		istringstream iss(szLine);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tData));
		set<int> tempset;
		for(vector<string>::iterator iter = tData.begin(); iter!=tData.end(); ++iter) {
			tempset.insert(atoi((*iter).c_str()));
		}
		int tempK = tempset.size();
		clusters[p] = new int[tempK + 1];
		clusters[p][0] = tempK;
		j=1;
		for (it_set=tempset.begin(); it_set != tempset.end(); it_set++) {
			clusters[p][j] = *it_set;
			state[clusters[p][j]][p] = 1;
			j++;
		}
		p++;
		tempset.clear();
	}
	
	cluster = new int*[num_nodes + 1];
	cluster[0] = new int[1]; 
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	cluster[0][0] = num_nodes;
	
	for (i = 1; i <= num_nodes; i++) {
		temp = 0;
		for (j = 1; j <= num_clusters; j++) { temp += state[i][j]; }
		cluster[i] = new int[temp + 1];
		cluster[i][0] = temp;
		q = 1;
		for (j = 1; j <= num_clusters; j++) {
			if (state[i][j]!=0) {
				cluster[i][q] = j; 
				q++;
			}
		}
	}
	
	// free memory space
	// state
	for(j = 0; j <= num_nodes; j++) { delete []state[j]; }
	delete []state;
}

void Graph::CalculateNumCommCommu() {
	int i = 1;
	int j = 1;
	
	NumCommCommu = new int*[num_nodes + 1];
	NumCommCommu[0] = new int[1]; 
	//do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	NumCommCommu[0][0] = num_nodes;
	
	for (i = 1; i <= num_nodes; i++) {
		NumCommCommu[i] = new int[num_nodes - i + 1];
		NumCommCommu[i][0] = num_nodes - i;
		for (j = i+1; j <= num_nodes; j++) {
			NumCommCommu[i][j-i] = NumIntersection(cluster[i], cluster[j]);
		}
		printf("\r#Common Communities Calculation Completed Progress: %.2lf%%(%d)", i / double(num_nodes) * 100, i);
		fflush(stdout);
	}
	cout << endl;
}

void Graph::CalCorrelation() {
	CalculateNumCommCommu();
	int i = 1;
	int j = 1;
	correlation = new double[num_clusters + 1];
	correlation[0] = 0.0;
	
	// store the number of node pairs of specific number of common clusters, 
	// for subsequent calculation of average
	int* NumNodePairs;
	NumNodePairs = new int[num_clusters + 1];
	NumNodePairs[0] = num_clusters;	
	// initialization
	for (i = 0; i <= num_clusters; i++) {
		correlation[i] = 0;
		NumNodePairs[i] = 0; 
	}	
	// summation
	for (i = 1; i <= num_nodes; i++) {
		for (j = i + 1; j <= num_nodes; j++) {
			NumNodePairs[NumCommCommu[i][j-i]]++;
			correlation[NumCommCommu[i][j-i]] += proximity[i][j-i];
		}
	}	
	// average
	for (i = 0; i <= num_clusters; i++) {
		correlation[i] /= NumNodePairs[i]; 
		if(isnan(correlation[i]) == 0) {
			cout << i << " common clusters: " << correlation[i] << endl;
		}
	}
}

void Graph::OutputCorrelation() {
	// Output the correlation to a file 
	char* graphtail = (char*)".correlation";
	char* graphfile = strcatre(name_dataset, graphtail);
	ofstream ffout;
	ffout.open(graphfile);
	ffout << "Type of Proximity: ";
	if (TypeTopolSimil == 1) { ffout << "Jaccard Coefficient" << endl; }
	else if (TypeTopolSimil == 2) { ffout << "Salton Index" << endl; }
	else if (TypeTopolSimil == 3) { ffout << "Sorensen Index" << endl; }
	else if (TypeTopolSimil == 4) { ffout << "Hub Promoted Index" << endl; }
	else if (TypeTopolSimil == 5) { ffout << "Hub Depressed Index" << endl; }
	else if (TypeTopolSimil == 6) { ffout << "Leicht-Holme-Newman Index" << endl; }
	ffout << "NumCommCommu   " << " Average_proximity_of_all_corresponding_node_pairs" << endl;
	for (int i = 0; i <= num_clusters; i++) {
		if(isnan(correlation[i]) == 0) { ffout << i << " " << correlation[i] << endl; }
	}
	ffout.close();
}

// visualize the graph (output a .gml file which can be fed into software named Cytoscape
void Graph::VisualizeGraph() {
	char* graphtail = (char*)".gml";
	char* graphfile = strcatre(name_dataset, graphtail);
	ofstream fout;
	fout.open(graphfile);
	fout << "graph" << endl;
	fout << "[" << endl;
	fout << "	directed 0" << endl;
	int i, j;
	for(i = 1; i<= num_nodes; i++) {
		fout << "	node" << endl;
		fout << "	[" << endl;
		fout << "		id " << i << endl;
		fout << "		label " << "\"" << i << "\"" << endl;
		// non-overlapping clusters
		fout << "		cluster " << "\"" << cluster[i][1] << "\"" << endl;
		fout << "	]" << endl;
	}
	
	for(i = 1; i <= num_nodes; i++) {
		for(j = 1; j <= neighbor[i][0]; j++) {
			if(i < neighbor[i][j]) {
				fout << "	edge" << endl;
				fout << "	[" << endl;
				fout << "		source " << i << endl;
				fout << "		target " << neighbor[i][j] << endl;
				fout << "	]" << endl;
			}
		}
	}
	fout << "]" << endl;
	fout.close();
}


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
	// each one-dimensional array in Our1 denotes a cluster, containing a set of nodes who belong it.
	// each one-dimensional array in Our2 is a set of clusters to which a specific node belongs.
	double CalNonOverlapModul(int** Our1, int** Our2, int** neighbor, int num_edges);
	// Modularity of overlapping cluster structure (Shen, et al. Physica A, 2009)
	double CalOverlapModul(int** Our1, int** Our2, int** neighbor, int num_edges);
	// Tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
	double CalTgt(int** Our1, int** Our2, int** neighbor);
	// Adjusted tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
	// penalizing very small and very large clusters and produces well-balanced solutions
	double CalAdjTgt(int** Our1, int** Our2, int** neighbor);
	
public:
	// has ground-truth (.gt file) information
	// 1. Average of f1-score (AvgF1)
	double CalculateF1(int* A, int* B);
	double CalculateMaxF1(int* A, int** targetset);
	double CalAvgF1(int** Our, int** GT);
	// 2. Normalized mutual information (NMI)
	double CalNMI(int** Our, int** GT, int num_nodes);
	// 3. Omega Index
	// each line in cluster_1 and cluster_2 is a set of clusters to which a specific node belongs.
	double CalOmegaIndex(int** cluster_1, int** cluster_2, int num_nodes);
	// 4. Adjusted Rand Index (ARI)
	double CalARI(int** Our, int** GT, int num_nodes);
};

// calculating modularity of non-overlapping cluster structure. (Newman, PRE, 2004.)
double ClusterTest::CalNonOverlapModul(int** Our1, int** Our2, int** neighbor, int num_edges) {
	// each one-dimensional array in Our1 denotes a cluster, containing a set of nodes who belongs it.
	// each one-dimensional array in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
	// if the cluster structure is non-overlapping, each one-dimensional array in Our2 only has two values, the first value is 1 which is useless.
	// if input is a overlapping cluster structure, this function will only take the first cluster for each node as her FuzzyMembership.
	// num_edges is the number of all edges in the graph.
	int i = 0;
	int k = 0;
	// number of all nodes in the graph
	int num_nodes = neighbor[0][0];
	
	double* theta = new double[Our1[0][0] + 1];
	theta[0] = Our1[0][0];
	for (k = 1; k <= theta[0]; k++) {
		theta[k] = 0;
		for (i = 1; i <= Our1[k][0]; i++) { theta[k] += neighbor[Our1[k][i]][0]; }
		theta[k] /= 2 * num_edges;
	}
	
	// modularity of single node
	double* singleModul = new double[num_nodes + 1];
	singleModul[0] = num_nodes;
	double tempValue = 0.0;
	double SumModul = 0.0;
	for (i = 1; i <= num_nodes; i++) {
		tempValue = 0.0;
		singleModul[i] = 0.0;
		// determine whether node-i does not belong to any cluster.
		if (Our2[i][0] == 0) { 
			singleModul[i] = 0; 
		} else if (neighbor[i][0] != 0) {
			tempValue = NumIntersection(neighbor[i], Our1[Our2[i][1]]) / (double)neighbor[i][0];
			singleModul[i] = (tempValue - theta[Our2[i][1]]) * neighbor[i][0] / 2.0 /num_edges;
		}
		// summation
		SumModul += singleModul[i];
	}
	return SumModul;
}

// calculating modularity of overlapping and hierarchical cluster structure.(Shen, et al. Physica A, 2009.)
double ClusterTest::CalOverlapModul(int** Our1, int** Our2, int** neighbor, int num_edges) {
// each one-dimensional array in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each one-dimensional array in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
	int i = 0;
	int j = 0;
	int k = 0;
	int id_i, id_j;
	int num_clusters = Our1[0][0];
	double* singleModul = new double[num_clusters + 1];
	singleModul[0] = num_clusters;
	double SumModul = 0.0;
	
	for (k = 1; k <= num_clusters; k++) {
		//cout << k << ": for test!!!" << endl;
		if (Our1[k][0] < 2) { singleModul[k] = 0.0; continue; }
		singleModul[k] = 0.0;
		for (i = 1; i < Our1[k][0]; i++) {
			id_i = Our1[k][i];
			if (Our2[id_i][0] == 0) { continue; }
			for (j = i+1; j <= Our1[k][0]; j++) {
				id_j = Our1[k][j];
				if (Our2[id_j][0] == 0) { continue; }
				if (find(neighbor[id_i], id_j) == 0) {
					singleModul[k] = -neighbor[id_i][0] * neighbor[id_j][0] / 2.0 / (double)num_edges / Our2[id_i][0] / Our2[id_j][0];
				} else { 
					singleModul[k] = (1.0 - neighbor[id_i][0] * neighbor[id_j][0] / 2.0 / (double)num_edges) / Our2[id_i][0] / Our2[id_j][0];
				}
			}
		}
		SumModul += singleModul[k];
	}
	SumModul = 0.5 * SumModul / num_edges;
	return SumModul;
}

// calculating tightness of cluster structure (Bu, et al. Information Fusion, 2017)
double ClusterTest::CalTgt(int** Our1, int** Our2, int** neighbor) {
	int i = 0;
	int j = 0;
	int k = 0;
	int q = 0;
	int flag = 0;
	double FirTerm = 0.0;
	double SecTerm = 0.0;
	int num_clusters = Our1[0][0];
	int num_nodes = neighbor[0][0];
	double SumTgt = 0.0;
	
	double* SinglTgt = new double[num_clusters + 1];
	SinglTgt[0] = num_clusters;
	int* numInEdges = new int[num_clusters + 1];
	numInEdges[0] = num_clusters;
	int* numOutEdges = new int[num_clusters + 1];
	numOutEdges[0] = num_clusters;
	
	// a cluster k
	for (k = 1; k <= num_clusters; k++) {
		SinglTgt[k] = 0.0;
		numInEdges[k] = 0;
		numOutEdges[k] = 0;
		// a node i in the cluster-k
		for (i = 1; i <= Our1[k][0]; i++) {
			// a neighbor of node i, i.e., Our1[k][i]
			for (j = 1; j <= neighbor[Our1[k][i]][0]; j++) {
				int* ngtClust = Our2[neighbor[Our1[k][i]][j]];
				// j's FuzzyMemberships, if she belongs cluster k, thus numInEdges[k] plus 1.
				// it can be applied to non-overlapping as well as overlapping clusters
				flag = 0;
				for (q = 1; q <= ngtClust[0]; q++) {
					if (ngtClust[q] == k) { 
						numInEdges[k]++; 
					} else if (flag != 1) { 
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
		if (Our1[k][0] != 0 && Our1[k][0] != num_nodes) {
			FirTerm = 2.0 * numInEdges[k] / Our1[k][0] / (double)Our1[k][0];
			SecTerm = numOutEdges[k] / Our1[k][0] / (double)(num_nodes - Our1[k][0]);
			SinglTgt[k] = FirTerm - SecTerm;
		}
		SumTgt += SinglTgt[k];
	}
	return SumTgt;
}

// calculating adjusted tightness of cluster structure (Bu, et al. Information Fusion, 2017)
// penalizing very small and very large clusters and produces well-balanced solutions
double ClusterTest::CalAdjTgt(int** Our1, int** Our2, int** neighbor) {
	int i = 0;
	int j = 0;
	int k = 0;
	int q = 0;
	int flag = 0;
	double FirTerm = 0.0;
	int num_clusters = Our1[0][0];
	int num_nodes = neighbor[0][0];
	double SumTgt = 0.0;
	
	double* SinglTgt = new double[num_clusters + 1];
	SinglTgt[0] = num_clusters;
	int* numInEdges = new int[num_clusters + 1];
	numInEdges[0] = num_clusters;
	int* numOutEdges = new int[num_clusters + 1];
	numOutEdges[0] = num_clusters;
	
	// a cluster k
	for (k = 1; k <= num_clusters; k++) {
		SinglTgt[k] = 0.0;
		numInEdges[k] = 0;
		numOutEdges[k] = 0;
		// a node i in the cluster-k
		for (i = 1; i <= Our1[k][0]; i++) {
			// a neighbor of node i, i.e., Our1[k][i]
			for (j = 1; j <= neighbor[Our1[k][i]][0]; j++) {
				int* ngtClust = Our2[neighbor[Our1[k][i]][j]];
				// j's FuzzyMemberships, if she belongs cluster k, thus numInEdges[k] plus 1.
				// it can be applied to non-overlapping as well as overlapping clusters
				flag = 0;
				for (q = 1; q <= ngtClust[0]; q++) {
					if (ngtClust[q] == k) { 
						numInEdges[k]++; 
					} else if (flag != 1) { 
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
		if (Our1[k][0] != 0) {
			FirTerm = 2.0 * (num_nodes - Our1[k][0]) * numInEdges[k] /  (double)Our1[k][0];
			SinglTgt[k] = FirTerm - numOutEdges[k];
		}
		SumTgt += SinglTgt[k];
	}
	return SumTgt;
}

// calculating F1
double ClusterTest::CalculateF1(int* A, int* B) {
	int K_A = A[0];
	int K_B = B[0];
	int nmuber_joint = NumIntersection(A, B);
	double precision = 0;
	double recall = 0;
	double f1 = 0;
	if (K_A != 0 && K_B != 0) {
		precision = nmuber_joint / K_A;
		recall = nmuber_joint / K_B;
	}
	if ((precision + recall) != 0 ) {
		f1 = 2 * precision * recall / (precision + recall);
		return f1;
	} else { return f1; }
}

// calculating maximum of F1
double ClusterTest::CalculateMaxF1(int* A, int** targetset) {
	int K_A = A[0];
	int Size_T = targetset[0][0];
	double maxf1 = 0.0;
	int q = 1;
	for (q = 1; q <= Size_T; q++) {
		// double ClusterTest::CalculateF1(int* A, int* B)
		double tempf1 = CalculateF1(A, targetset[q]);
		if (tempf1 >= maxf1) { maxf1 = tempf1; }
	}
	return maxf1;
}

// calculating the average of f1-score
double ClusterTest::CalAvgF1(int** Our, int** GT) {
	int K_1 = Our[0][0];
	int K_2 = GT[0][0];
	double* MaxF1_1 = new double[K_1 + 1];
	MaxF1_1[0] = (double)K_1;
	double* MaxF1_2 = new double[K_2 + 1];
	MaxF1_2[0] = (double)K_2;
	int p = 0;
	int q = 0;
	double sumMAXF1_1 = 0.0;
	double sumMAXF1_2 = 0.0;
	for (p = 1; p <= K_1; p++) {
		// double ClusterTest::CalculateMaxF1(int* A, int** targetset)
		MaxF1_1[p] = CalculateMaxF1(Our[p], GT);
		sumMAXF1_1 += MaxF1_1[p];
	}
	for (q = 1; q <= K_2; q++) {
		MaxF1_2[q] = CalculateMaxF1(GT[q], Our);
		sumMAXF1_2 += MaxF1_2[q];
	}	
	double AvgF1 = sumMAXF1_1 / 2 / K_1 + sumMAXF1_2 / 2 / K_2;
	return AvgF1;
}

// calculating the normalized mutual information (NMI)
// [S. Fortunato. Community detection in graphs. Physics Reports, 486(3-5):75 – 174, 2010.]
double ClusterTest::CalNMI(int** Our, int** GT, int num_nodes) {
	int size_a = Our[0][0];
	int size_b = GT[0][0];
	int i = 1;
	int j = 1;
	
	double** number_joint = new double*[size_a + 1];
	number_joint[0] = new double[1]; 
	// do not ignore it !!! otherwise, "core dumped" happens when free memory space.
	number_joint[0][0] = (double)size_a;
	
	for (i = 1; i <= size_a; i++) {
		number_joint[i] = new double[size_b + 1];
		for (j = 1; j <= size_b; j++) {
			number_joint[i][j] = (double)NumIntersection(Our[i], GT[j]);
		}
	}

	double NMI = 0.0;
	double* NMI_A = new double[size_a + 1];
	double* NMI_B = new double[size_b + 1];
	
	for (i = 1; i <= size_a; i++) { NMI_A[i] = 0.0; }
	NMI_A[0] = size_a;
	for (j = 1; j <= size_b; j++) { NMI_B[j] = 0.0; }
	NMI_B[0] = size_b;
	
	for (i = 1; i <= size_a; i++) {
		for (j = 1; j <= size_b; j++) { NMI_A[i] += number_joint[i][j]; }
	}
	for (j = 1; j <= size_b; j++) {
		for (i = 1; i <= size_a; i++) { NMI_B[j] += number_joint[i][j]; }
	}
	for (i = 1; i <= size_a; i++) {
		for (j = 1; j <= size_b; j++) {
			if (number_joint[i][j] != 0) {
				NMI += number_joint[i][j] * log(number_joint[i][j] * num_nodes / NMI_A[i] / NMI_B[j]) / log(2);
			}
		}
	}
	NMI = NMI * (-1 * 2);

	double sum_NMI_A = 0.0;
	double sum_NMI_B = 0.0;
	for(i = 1; i <= size_a; i++) {
		if (NMI_A[i] != 0) {
			sum_NMI_A += NMI_A[i] * log(NMI_A[i] / num_nodes) / log(2);
		}
	}
	for(j = 1; j <= size_b; j++) {
		if(NMI_B[j] != 0) { sum_NMI_B += NMI_B[j] * log(NMI_B[j] / num_nodes) / log(2); }		
	}
	NMI = NMI / (sum_NMI_A + sum_NMI_B);
	return NMI;
}

// calculating the value of Omega index
// [S. Gregory. Fuzzy overlapping communities in networks. J. of Stat. Mech.: Theory and Experiment, 2011.]
double ClusterTest::CalOmegaIndex(int** cluster_1, int** cluster_2, int num_nodes) {
	// each line in cluster_1 and cluster_2 denotes a set of clusters to which a specific node belongs.
	int i = 1;
	int j = 1;
	int sum_temp = 0;
	double results = 0.0;
	for (i = 1; i <= num_nodes; i++) {
		for (j = i + 1; j <= num_nodes; j++) {
			if (NumIntersection(cluster_1[i], cluster_1[j]) == NumIntersection(cluster_2[i], cluster_2[j])) { sum_temp++; }
		}
	}
	results = (double)sum_temp / num_nodes / num_nodes;
	return results;
}

// calculating Adjusted Rand Index (ARI)
double ClusterTest::CalARI(int** Our, int** GT, int num_nodes) {
	int size_a = Our[0][0];
	int size_b = GT[0][0];
	int i = 1;
	int j = 1;
	double numerator = 0.0;
	double denominator = 0.0;
	double temp_value = 0.0;
	double numerator_1 = 0.0;
	double numerator_2 = 0.0;
	
	for (i = 1; i <= size_a; i++) {
		for (j = 1; j <= size_b; j++) {
			temp_value = NumIntersection(Our[i], GT[j]);
			numerator_1 = numerator_1 + temp_value * (temp_value - 1) / 2.0;
		}
	}
	
	double temp_value_1 = 0;
	for (i = 1; i <= size_a; i++) {
		temp_value_1 += Our[i][0] * (Our[i][0] - 1) / 2.0;
	}
	double temp_value_2 = 0;
	for (j = 1; j <= size_b; j++) {
		temp_value_2 += GT[j][0] * (GT[j][0] - 1) / 2.0;
	}	
	
	numerator_2 = (temp_value_1 + temp_value_2) * 2.0 / num_nodes / (num_nodes - 1);
	denominator = (temp_value_1 + temp_value_2) / 2.0 - numerator_2;
	double results = (numerator_1 - numerator_2) / denominator;
	return results;
}


// class of BIGCLAM (Jaewon Yang, Jure Leskovec. WSDM2013)
class BIGCLAM {
public:
	// all nodes' FuzzyMembership vector
	double** FuzzyMemb;
	// cluster size
	double* SizeClus;
	// step size (learn rate)
	double SizeStep;
	// gradient of FuzzyMemb
	double** Grad;
	// temporal FuzzyMemb for verify
	double** TempFuzzyMemb;
	// total likelihood
	double TotalLikelihood;
	// Delta for determine the FuzzyMembership according to FuzzyMemb
	double Delta;
	// Memb in clusters (similar with mygraph.clusters)
	int** MembClus;
	// Memb in nodes (similar with mygraph.cluster)
	int** MembNode;
	// maximum number of iteration
	int MaxIter;
	
public:
	// constructor (initialization, file name is needed)
	BIGCLAM(double sizestep, int maxiter, double delta): SizeStep(sizestep), MaxIter(maxiter), Delta(delta) { }
	// randomly initialize FuzzyMemb, TempFuzzyMemb, Grad, also the SizeClus is derived
	void IniVar(Graph mygraph);
	// calculate the gradient
	void CalGrad(Graph mygraph);
	// update FuzzyMemb
	void CalTempFuzzyMemb(Graph mygraph);
	// calculate total likelihood
	double CalLikelihood(Graph mygraph);
	// calculate temporal total likelihood
	double CalTempLikelihood(Graph mygraph);
	// generate crisp clusters according FuzzyMemb and Delta
	void GenClus(Graph mygraph);
	// perform BIGCLAM
	void performBIGCLAM(Graph mygraph);
};

// randomly initialize FuzzyMemb, TempFuzzyMemb, Grad; SizeClus is also derived.
void BIGCLAM::IniVar(Graph mygraph) {
	int i, k;
	int num_nodes = mygraph.num_nodes;
	int num_clusters = mygraph.num_clusters;
	FuzzyMemb = new double*[num_nodes + 1];
	FuzzyMemb[0] = new double[1];
	FuzzyMemb[0][0] = num_nodes;
	for (int i = 1; i <= num_nodes; i++) {
		FuzzyMemb[i] = new double[num_clusters + 1];
		FuzzyMemb[i] = randArrSumOne(num_clusters);
	}
	// derive the size of each clusters
	SizeClus = new double[num_clusters + 1];
	SizeClus[0] = num_clusters;
	for (k = 1; k <= num_clusters; k++) {
		SizeClus[k] = 0.0;
		for (i = 1; i <= num_nodes; i++) { SizeClus[k] += FuzzyMemb[i][k]; }
	}
	// initialize TempFuzzyMemb, begin by making same with FuzzyMemb
	// space must be allocated and cannot be directly assigned as a whole part
	TempFuzzyMemb =  new double*[num_nodes + 1];
	TempFuzzyMemb[0] = new double[1];
	TempFuzzyMemb[0][0] = num_nodes;
	for (i = 1; i <= num_nodes; i++) {
		TempFuzzyMemb[i] = new double[num_clusters + 1];
		TempFuzzyMemb[i][0] = num_clusters;
		for (k = 1; k <= num_clusters; k++) { TempFuzzyMemb[i][k] = FuzzyMemb[i][k]; }
	}
	// initialize the Grad
	Grad =  new double*[num_nodes + 1];
	Grad[0] = new double[1];
	Grad[0][0] = num_nodes;
	for (i = 1; i <= num_nodes; i++) {
		Grad[i] = new double[num_clusters + 1];
		Grad[i][0] = num_clusters;
		for (k = 1; k <= num_clusters; k++) { Grad[i][k] = 0.0; }
	}
	// initialize Delta
	// Delta = sqrt(-log(1-pow(10,-8)));
	// allocate space for MembClus
	MembClus = new int*[num_clusters + 1];
	MembClus[0] = new int[1];
	MembClus[0][0] = num_clusters;
	// allocate space for MembNode
	MembNode = new int*[num_nodes + 1];
	MembNode[0] = new int[1];
	MembNode[0][0] = num_nodes;
}

// calculate the gradient
void BIGCLAM::CalGrad(Graph mygraph) {
	int i, j; // for node index
	int k, q; // for cluster index
	double TempValue = 0.0;
	double FirTerm = 0.0;
	double SecTerm = 0.0;
	// calculate gradient
	for (i = 1; i <= mygraph.num_nodes; i++) {
		for (q = 1; q <= mygraph.num_clusters; q++) {
			FirTerm = 0.0;
			for (j = 1; j <= mygraph.neighbor[i][0]; j++) {
				TempValue = 0.0;
				for (k = 1; k <= mygraph.num_clusters; k++) { TempValue += FuzzyMemb[i][k] * FuzzyMemb[mygraph.neighbor[i][j]][k]; }
				if (TempValue == 0) { 
					FirTerm += FuzzyMemb[mygraph.neighbor[i][j]][q] / pow(10, -7); 
				} else {
					FirTerm += FuzzyMemb[mygraph.neighbor[i][j]][q] * exp(-TempValue) / (1 - exp(-TempValue));
				}				
			}
			TempValue = 0.0;
			SecTerm = 0.0;
			for (j = 1; j <= mygraph.neighbor[i][0]; j++) { TempValue += FuzzyMemb[mygraph.neighbor[i][j]][q]; }
			SecTerm = SizeClus[q] - FuzzyMemb[i][q] - TempValue;
			Grad[i][q] = FirTerm - SecTerm;
		}
		//printf("\r#Gradient Calculation Completed Progress: %.2lf%%(%d)", i / double(mygraph.num_nodes) * 100, i);
		//fflush(stdout);
	}
}

// calculate the TempFuzzyMemb
void BIGCLAM::CalTempFuzzyMemb(Graph mygraph) {
	int i, k;
	double TempValue = 0.0;
	for (i = 1; i <= mygraph.num_nodes; i++) {
		for (k = 1; k <= mygraph.num_clusters; k++) {
			TempValue = 0.0;
			TempValue = FuzzyMemb[i][k] + SizeStep * (Grad[i][k]);
			if (TempValue < 0) {
				TempFuzzyMemb[i][k] = pow(10, -8);
			} else { TempFuzzyMemb[i][k] = TempValue; }
		}
	}
	//normalization
	for (i = 1; i <= mygraph.num_nodes; i++) {
		double temp = 0.0;
		for (k = 1; k <= mygraph.num_clusters; k++) { temp += TempFuzzyMemb[i][k]; }
		for (k = 1; k <= mygraph.num_clusters; k++) { TempFuzzyMemb[i][k] /= temp; }
	}
}

// calculate total likelihood
double BIGCLAM::CalLikelihood(Graph mygraph) {
	int num_nodes = mygraph.num_nodes;
	int num_clusters = mygraph.num_clusters;
	int** neighbor = mygraph.neighbor;
	int i, j; // for node index
	int k; // for cluster index
	double TempValue = 0.0;
	double FirTerm = 0.0;
	double SecTerm = 0.0;
	double TotalLikelihood = 0.0;
	double* TempVec = new double[num_clusters + 1];
	TempVec[0] = num_clusters;
	for (i = 1; i <= num_nodes; i++) {
		FirTerm = 0.0;
		for (j = 1; j <= neighbor[i][0]; j++) { 
			TempValue = 0.0;
			for (k = 1; k <= num_clusters; k++) { TempValue += FuzzyMemb[i][k] * FuzzyMemb[neighbor[i][j]][k]; }
			FirTerm += log(1 - exp(-TempValue));
		}
		SecTerm = 0.0;
		for (k = 1; k <= num_clusters; k++) {
			for (j = 1; j <= neighbor[i][0]; j++) { TempValue += FuzzyMemb[neighbor[i][j]][k]; }
			TempVec[k] = SizeClus[k] - FuzzyMemb[i][k] - TempValue;
			SecTerm += FuzzyMemb[i][k] * TempVec[k];
		}
		//cout << FirTerm << " ";
		//cout << SecTerm << endl;
		TotalLikelihood += FirTerm - SecTerm;
		printf("\r#likelihood Calculation Completed Progress: %.2lf%%(%d)", i / double(num_nodes) * 100, i);
		fflush(stdout);
	}
	cout << endl;
	cout << "Total likelihood: " << TotalLikelihood << endl;
	return TotalLikelihood;
}

// calculate temporal total likelihood
double BIGCLAM::CalTempLikelihood(Graph mygraph) {
	int num_nodes = mygraph.num_nodes;
	int num_clusters = mygraph.num_clusters;
	int** neighbor = mygraph.neighbor;
	int i, j; // for node index
	int k; // for cluster index
	double TempValue = 0.0;
	double FirTerm = 0.0;
	double SecTerm = 0.0;
	double newTotalLikelihood = 0.0;
	double* TempVec = new double[num_clusters + 1];
	TempVec[0] = num_clusters;
	for (i = 1; i <= num_nodes; i++) {
		FirTerm = 0.0;
		for (j = 1; j <= neighbor[i][0]; j++) { 
			TempValue = 0.0;
			for (k = 1; k <= num_clusters; k++) { TempValue += TempFuzzyMemb[i][k] * TempFuzzyMemb[neighbor[i][j]][k]; }
			FirTerm += log(1 - exp(-TempValue));
		}
		SecTerm = 0.0;
		for (k = 1; k <= num_clusters; k++) { 
			for (j = 1; j <= neighbor[i][0]; j++) { TempValue += TempFuzzyMemb[neighbor[i][j]][k]; }
			TempVec[k] = SizeClus[k] - TempFuzzyMemb[i][k] - TempValue;
			SecTerm += TempFuzzyMemb[i][k] * TempVec[k];
		}
		//cout << FirTerm << " ";
		//cout << SecTerm << endl;
		newTotalLikelihood += FirTerm - SecTerm;
		//printf("\r#CalLikelihood: %.2lf%%(%d)", i / double(num_nodes) * 100, i);
		fflush(stdout);
	}
	return newTotalLikelihood;
}

// generate cluster
void BIGCLAM::GenClus(Graph mygraph) {
	int temp = 0;
	int i = 0;
	int k = 0;
	int c = 0;
	for (i = 1; i <= mygraph.num_nodes; i++) {
		temp = 0;
		// determine the space size for space allocation
		for (k = 1; k <= mygraph.num_clusters; k++) {
			if (FuzzyMemb[i][k] > Delta) { temp++; }
		}
		// allocate space
		MembNode[i] = new int[temp + 1];
		MembNode[i][0] = temp; 
		c = 1;
		for (k = 1; k <= mygraph.num_clusters; k++) {
			if (FuzzyMemb[i][k] > Delta) { 
				MembNode[i][c] = k; 
				c++;
			}
		}
	}
	
	for (k = 1; k <= mygraph.num_clusters; k++) {
		temp = 0;
		// determine the space size for space allocation
		for (i = 1; i <= mygraph.num_nodes; i++) {
			if (FuzzyMemb[i][k] > Delta) { temp++; }
		}
		// allocate space
		MembClus[k] = new int[temp + 1];
		MembClus[k][0] = temp; 
		c = 1;
		for (i = 1; i <= mygraph.num_nodes; i++) {
			if (FuzzyMemb[i][k] > Delta) { 
				MembClus[k][c] = i; 
				c++;
			}
		}
	}
}	

// perform BIGCLAM
void BIGCLAM::performBIGCLAM(Graph mygraph) {
	int i, k;
	IniVar(mygraph);
	double oldLikelihood = 0.0;
	double newLikelihood = CalLikelihood(mygraph);
	double change = 0.00001;
	
	int iter = 1;
	while (change >= 0.00001 && iter <= MaxIter) {
		cout << iter << "-th iteration: ";
		oldLikelihood = newLikelihood;
		FuzzyMemb = TempFuzzyMemb;
		for (k = 1; k <= mygraph.num_clusters; k++) {
			SizeClus[k] = 0.0;
			for (i = 1; i <= mygraph.num_nodes; i++) { SizeClus[k] += FuzzyMemb[i][k]; }
		}
		CalGrad(mygraph);
		CalTempFuzzyMemb(mygraph);
		newLikelihood = CalTempLikelihood(mygraph);
		cout << "newLikelihood: " << newLikelihood << endl;
		change = abs((newLikelihood - oldLikelihood) / oldLikelihood); 
		
		iter++;
	}
	GenClus(mygraph);
	
	/*for (i = 1; i <= mygraph.num_nodes; i++) {
		for (k = 1; k <= mygraph.num_clusters; k++) {
			// cout << FuzzyMemb[i][k] << ":(" << Grad[i][k] << "," << TempFuzzyMemb[i][k] << ") ";
			if (FuzzyMemb[i][k] > Delta) { 
				cout << 1 << " ";
			} else { cout << 0 << " "; }
		}
		cout << endl;
	}
	for (k = 1; k <= mygraph.num_clusters; k++) {
		for (i = 1; i <= mygraph.num_nodes; i++) {
			// cout << FuzzyMemb[i][k] << ":(" << Grad[i][k] << "," << TempFuzzyMemb[i][k] << ") ";
			if (FuzzyMemb[i][k] > Delta) { cout << i << " "; }
		}
		cout << endl;
	}
	for (i = 1; i <= mygraph.num_nodes; i++) {
		cout << "(node-" << i << ") ";
		for (k = 1; k <= MembNode[i][0]; k++) {
			cout << MembNode[i][k] << " ";
		}
		cout << endl;
	}*/
}
	

int main(int argc, char **argv)
{
	char* head = argv[1]; // name of dataset
	// type of 2nd-order proximity: 1. Jaccard Coefficient; 2. Salton Index; 3. Sorensen Index; 
	// 4. Hub Promoted Index; 5. Hub Depressed Index; 6. Leicht-Holme-Newman Index
	int TypeProx = atoi(argv[2]);
	float stepsize = atof(argv[3]);
	int maxiter = atoi(argv[4]);
	float delta = atof(argv[5]);
	
	Graph mygraph(head); // initialize an object
	struct timeval tod1, tod2; // record time
	
// Phase: Loading data
	cout << "========= " << "Load graph: " << head << " ================== "<< endl;
	gettimeofday(&tod1, NULL);		
	mygraph.ReadPureDirGraph(); // load .graph file (pure directed graph)
	//mygraph.ReadNodeFeature(); // load .nf file (node feature)
	mygraph.ReadGTclusters(); // load .gt file (cluster ground-truth)
	//VisualizeGraph(mygraph); //output .gml file
	gettimeofday(&tod2, NULL);
	cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
	cout << endl;

//Tools: Calculate second-order proximity and correlation between common clusters of node pairs
	// 1. Jaccard Coefficient; 2. Salton Index; 3. Sorensen Index; 
	// 4. Hub Promoted Index; 5. Hub Depressed Index; 6. Leicht-Holme-Newman Index
	cout << "========= " << "Calculate 2nd proximity" << " ================== "<< endl;
	gettimeofday(&tod1, NULL);
	mygraph.CalSecOrdProx(TypeProx); //CalSecOrdProx(int type)
	//mygraph.CalCorrelation();
	//mygraph.OutputCorrelation();
	gettimeofday(&tod2, NULL);
	cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
	cout << endl;
		
// Tools: BIGCLAM
	cout << "========= " << "Node clustering by BIGCLAM." << "================== " << endl;
	BIGCLAM myBIGCLAM(stepsize, maxiter, delta);
	gettimeofday(&tod1, NULL);
	myBIGCLAM.performBIGCLAM(mygraph);
	gettimeofday(&tod2, NULL);
	cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
	cout << endl;
	
//Tools: Performance Test on nodes clustering
	cout << "========= " << "Performance test on nodes clustering." << "================== "<< endl;
	ClusterTest mytest;
	gettimeofday(&tod1, NULL);
	// no gt
	cout << "##### no ground-truth information:" << endl;
	//cout << "Modularity (no overlapping) = " << mytest.CalNonOverlapModul(mygraph.clusters, mygraph.cluster, mygraph.neighbor, mygraph.num_edges) << endl;
	cout << "Modularity (overlapping) = " << mytest.CalOverlapModul(myBIGCLAM.MembClus, myBIGCLAM.MembNode, mygraph.neighbor, mygraph.num_edges) << endl;
	cout << "Tightness = " << mytest.CalTgt(myBIGCLAM.MembClus, myBIGCLAM.MembNode, mygraph.neighbor) << endl;
	cout << "Adjusted Tightness = " << mytest.CalAdjTgt(myBIGCLAM.MembClus, myBIGCLAM.MembNode, mygraph.neighbor) << endl;
	cout << endl;
	// with gt
	cout << "##### with ground-truth information:" << endl;
	cout << "AvgF1 = " << mytest.CalAvgF1(myBIGCLAM.MembClus, mygraph.clusters) << endl;
	cout << "NMI = " << mytest.CalNMI(myBIGCLAM.MembClus, mygraph.clusters, mygraph.num_nodes) << endl;
	cout << "Omega Index = " << mytest.CalOmegaIndex(myBIGCLAM.MembNode, mygraph.cluster, mygraph.num_nodes) << endl;
	cout << "ARI = " << mytest.CalARI(myBIGCLAM.MembClus, mygraph.clusters, mygraph.num_nodes) << endl;
	gettimeofday(&tod2, NULL);
	cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
	
	return 0;
}
