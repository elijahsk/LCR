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
class DGraph : public Graph
{
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
        DGraph(SmallEdgeSets inE， SmallEdgeSets outE, int pN, int pL, int pM, vector<int> weight);
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
        SmallEdgeSets getInE();
        SmallEdgeSets getOutE();

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
        int getWeight(VertexID v);

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
        void newClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID);
        vector<VertexID> findLongestChain();
        vector<VertexID> findLongestChain(VertexID v);
        void tarjan(vector< vector<VertexID> >& SCCs);
        void tarjanStrongConnect(int v, int& index, stack<VertexID>& q, vector< int >& indexPerNode,
            vector< int >& lowlinkPerNode, vector< bool >& onStack, vector< vector<VertexID> >& SCCs);
        int computeDiameter();

};
#endif
