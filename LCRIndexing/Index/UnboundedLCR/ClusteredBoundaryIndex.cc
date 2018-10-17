#include <cstdlib>
#include <set>

using namespace std;

class ClusteredBoundaryIndex {
	public: 
		void construct(Graph);

	private:
		vector<int> vertexToClusterIndex;
		vector<set<int>> clusterToVertexIndex;
		int vertexCount;
		int edgeCount;
		int clusterCount;
		vector<set<int>> clusterLabels;
		vector<set<int>> clusterBoundaryNodes;
}

void ClusteredBoundaryIndex::construct(Graph* mg) {
	this->graph = mg;
	this->vertexCount = this->graph.getNumberOfVertices();
	this->edgeCount = this->graph.getNumberOfEdges();
	

	constStartTime = getCurrentTimeInMilliSec();
    
    createClustering(method);

    buildIndex();

    constEndTime = getCurrentTimeInMilliSec();
}

void ClusteredBoundaryIndex::createClustering(int method) {
	switch (method) {
		case 1: randomClustering();
				break;
	}
}

void randomClustering() {
	// arbitrary num of clusters
	this->clusterCount = this->graph.getNumberOfVertices() / 4;
	int mergingEdgesCount = this->graph.getNumberOfVertices() - numOfClusters;

	// initialize cluster index to be unique
	for (int i = 0; i < this->vertexCount; i++) {
		this->vertexToClusterIndex.push_back(i);
		set<int> temp;
		temp.insert(i);
		this->clusterToVertexIndex.push_back(temp);
	}


	set<int> mergedEdges;

	for (int i = 0; i < mergingEdgesCount; i++) {
		bool hasMerged = false;
		while (!hasMerged):
			edgeNum = randn() % this->vertexCount;
			if (mergedEdges.find(edgeNum) == mergedEdges.end()) {
				Edge selectedEdge = this->graph->EdgeSet.at(edgeNum)
				int sourceCID = selectedEdge.first;
				int targetCID = selectedEdge.first.first;
				mergeCluster(sourceCID, targetCID);
				hasMerged = true;
			}
	}

	// relabel the clusters by consecutive index 0, 1, 2, ...
	reLabelClusters()
}

void mergeCluster(int sourceCID, int targetCID) {
	vector<int> targetVIDs = this->clusterToVertexIndex.at(targetCID);
	
	for (vector<int>::iterator targetVID = targetVIDs.begin(); targetVID != targetVIDs.end(); targetVID ++) {
		// modify vertexToCluster
		this->vertexToClusterIndex.at(targetVID) = sourceCID;
		// modify ClusterToVertex
		this->clusterToVertexIndex.at(sourceCID).insert(targetVID);
	}

	// clear targetCID
	this->clusterToVertexIndex.at(targetCID).clear();
}

void reLabelClusters() {
	int headCID = 0;
	int tailCID = this->vertexCount;
	vector<set<int>>::iterator head = this->clusterToVertexIndex.begin();
	vector<set<int>>::iterator tail = this->clusterToVertexIndex.end();
	while (headCID < tailCID) {
		// find first empty cluster from head
		while (!head.empty()) { ++head; ++headCID; }
		// find first non-empty cluster from tail
		while (tail.empty()) { --tail; --tailCID; }
		// set all cID of v in tail to headCID
		for (set<int>::iterator vID = tail.begin(); vID != tail.end(); vID ++) {
			this->vertexToClusterIndex.at(vID) = headCID
		}
		// set head to tail vIDs 
		this->clusterToVertexIndex.at(head) = this->clusterToVertexIndex.at(tail);
		// clear tail vIDs
		this->clusterToVertexIndex.at(tail.clear())
		
		++ head;
		-- tail;
	}

}

void buildIndex() {
	getAllLabelsAndBoundaryNodesPerCluster()
	buildRBI()
	buildRRBI()
	buildRCI()
}

void getAllLabelsAndBoundaryNodesPerCluster() {
	for (int cID = 0; cID < this->clusterCount; cID ++) {
		set<int> vIDs = this->clusterToVertexIndex.at(cID);
		set<int> labels;
		set<int> boundaryNodes;
		for (set<int>::iterator vID = vIDs.begin(); vID != vIDs.end(); vID ++) {
			SmallEdgeSet outNeighbours;
			this->graph->getOutNeighbours(vID, outNeighbours);
			for (SmallEdgeSet::iterator outNeighbour = outNeighbours.begin(); outNeighbour != outNeighbours.end(); outNeighbour ++) {
				// insert label to label set
				if (vIDs.find(outNeightbour.first) != vIDs.end()) {
					labels.insert(outNeighbour.second);
				}
				// insert boundary nodes, since it has edges that connects to other clusters
				else {
					boundaryNodes.insert(vID);
				}
			}
		}
		this->clusterLabels.push_back(labels);
		this->clusterBoundaryNodes.push_back(boundaryNodes);
	}
}


// reacheable boundary index
void buildRBI() {
	for (int cID = 0; cID < this->clusterCount; cID ++) {
		set<int> vIDs = this->clusterToVertexIndex.at(cID);
		set<int> bIDs = this->clusterBoundaryNodes.at(cID);
		for (set<int>::iterator vID = vIDs.begin(); vID != vIDs.end(); ++vID) {
				findMinLabelSet(vID, bIDs);
		}
	}
}
