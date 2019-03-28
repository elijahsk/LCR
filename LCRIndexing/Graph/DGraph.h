#include "Graph.h"

#include <set>
#include <stack>
#include <boost/dynamic_bitset.hpp>

using namespace graphns;
using namespace boost;

#ifndef DGRAPH_H
#define DGRAPH_H

/**
* This class contains a normal graph, just two lists with In and Out edges
* per node. It is a directed graph.
*/
class DGraph : public Graph {
private:
    SmallEdgeSets outE, inE;
    double constStartTime, constEndTime;
    vector<int> weight;
    bool allowMultipleEdges;

public:
    DGraph(string fileName);
    DGraph(EdgeSet* edgeSet);
    DGraph(EdgeSet* edgeSet, int pN, int pL);
    DGraph(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges);
    DGraph(SmallEdgeSets inE, SmallEdgeSets outE, int pN, int pL, int pM, vector<int> weight);
    void construct(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges);
    ~DGraph();

    void buildGraph(graphns::EdgeSet* edges);
    EdgeSet* loadEdgeFile(string fileName);
    void loadEdgeStats(graphns::EdgeSet* edgeSet);

    int getGraphSizeInBytes();
    double getGraphConstructionTime();

    void getOutNeighbours(graphns::VertexID w, SmallEdgeSet& outNeighbours);
    void getInNeighbours(graphns::VertexID w, SmallEdgeSet& inNeighbours);
    void getAllNeighbours(graphns::VertexID w, SmallEdgeSet& allNeighbours);

    // prints stats of the graph
    std::string getStats();
    std::string toString();

    int getNumberOfVertices();
    int getNumberOfLabels();
    int getNumberOfEdges();
    void getInE(SmallEdgeSets& tempInE);
    void getOutE(SmallEdgeSets& tempOutE);
    void getInE(VertexID v, SmallEdgeSet& tempInE);
    void getOutE(VertexID v, SmallEdgeSet& tempOutE);

    void addNode();
    void removeNode(graphns::VertexID w);
    void addEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel);
    void addMultiEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelSet newLabelSet);
    void removeEdge(graphns::VertexID v, graphns::VertexID w);
    void removeInEdges(graphns::VertexID v);
    void removeOutEdges(graphns::VertexID v);
    void changeLabel(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel);
    bool hasMultiEdge(graphns::VertexID v , graphns::VertexID w, graphns::LabelSet ls);
    void setWeight(vector<int> w);
    void setWeight(VertexID v, int w);
    void getWeight(VertexID v, int& tempWeight);
    void getWeights(vector<int>& weights);

    bool findInsertablePosition(graphns::VertexID w, graphns::SmallEdgeSet& ses, int& pos);
    void insertEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel, SmallEdgeSet& ses);
    LabelID getLabelID(graphns::VertexID v , graphns::VertexID w);
    long getCountPerLabel(graphns::LabelID l);

    bool hasEdge(VertexID v, VertexID w);

    double computeAverageDiameter();
    int computerNumberOfTriangles();
    int computeNumberOfConnectedTriplets();
    double computeClusterCoefficient();

    void initializeUnionFind(vector<VertexID>& parent);
    int find(vector<VertexID> parent, VertexID v);
    void connect(vector<VertexID>& parent, VertexID v, VertexID w);
    void randomClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID);
    
    void getNodesWithNoInDegree(DGraph* tempGraph, vector<VertexID>& nodesWithNoInDegree);
    void findLongestChain(DGraph* tempGraph, VertexID v, vector<VertexID>& chain);
    void findLongestChain(DGraph* tempGraph, vector<VertexID> nodesWithNoInDegree, vector<VertexID>& chain);
    void segmentChain(vector<VertexID> chain, int radius, vector<int> weights, vector<vector<VertexID>>& segments);
    void addSegments(int head, int tail, vector<VertexID> chain, vector<vector<VertexID>>& segments, int radius);
    void growSegment(DGraph* tempGraph, VertexID cID, vector<VertexID>& startVertices, int maxClusterSize, vector<vector<VertexID>>& clusters, vector<int>& vToCID);
    void modifyGraph(DGraph* tempGraph, VertexID v, VertexID w, vector<vector<VertexID>>& clusters, vector<int>& vToCID);
    void modifyGraph(DGraph* tempGraph, int segmentCount, vector<vector<VertexID>>& clusters, vector<int>& vToCID);
    void newClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID, int radius, int maxClusterSize, int minClusterSize);


    void tarjan(vector< vector<VertexID> >& SCCs);
    void tarjanStrongConnect(int v, int& index, stack<VertexID>& q, vector< int >& indexPerNode,
                             vector< int >& lowlinkPerNode, vector< bool >& onStack, vector< vector<VertexID> >& SCCs);
    int computeDiameter();

    void getStatus();

};
#endif
