#include "Index.h"
#include "../../Graph/DGraph.h"

#include <boost/dynamic_bitset.hpp>
#include <map>
#include <unordered_map>
#include <unordered_set>

using namespace graphns;
using namespace indexns;
using namespace boost;
using namespace std;

#ifndef NEW_H
#define NEW_H

    namespace newns
    {
        typedef struct
        {
            LabelSet ls;
            VertexID x;
            int dist;
        }
        NewIndexBitEntry;

        typedef vector< Triplet > BitEntries;

        struct compBitEntries
        {
           bool operator()( const NewIndexBitEntry& lhs,  const VertexID& rhs ) const
           {
               return lhs.x < rhs;
           }

           bool operator()( const VertexID& lhs, const NewIndexBitEntry& rhs ) const
           {
               return lhs > rhs.x;
           }
        };

        struct NewIndexPQBitEntries
        {
            bool operator()(NewIndexBitEntry const & t1, NewIndexBitEntry const & t2)
            {
                // return "true" if "p1" is ordered before "p2", for example:
                return t1.dist > t2.dist;
            }
        };
    }

    class NewIndex : public Index
    {

        public:
            NewIndex(Graph* mg);
            ~NewIndex();
            void initializeLocalIndexes();
            void buildIndex();
            void labeledBFSPerCluster(int cID, Graph* graph);
            void labeledBFSPerVertex(int cID, VertexID v, Graph* graph, vector<bool>& indexed);
            void labeledBFSAcrossClusters(int cID, vector<vector<VertexID>> clusters, vector<bool>& BNIndexed);
            void getRBI(int cID, Graph* graph, vector<vector<VertexID>> clusters);
            void getRRBI(int cID, Graph* graph, vector<vector<VertexID>> clusters);
            void getRRCI(int cID);
			void printNewRBI();
            void printRBI();
            void printRRBI();
            void printRRCI();
            void printTIn(int cID);

            bool query(VertexID source, VertexID target, LabelSet ls);
            bool queryShell(VertexID source, VertexID target, LabelSet ls);

            bool tryInsert(VertexID w, VertexID v, LabelSet ls);
            unsigned long getIndexSizeInBytes();

            void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);

        private:
            vector<Graph*> subGraphs;
            vector<int> vToSubGraphID;
            vector<int> vToCID;
            map<VertexID, int> positionInC;
            vector<bool> isBoundaryNode;
            vector<vector<VertexID>> boundaryNodesPerCluster;
            vector<VertexID> boundaryNodes;
            Graph* D;
            vector<vector<VertexID>> clusters;

            // each node -> a list of boundary nodes and respective minimal labelset
			// boundary nodes that a node can reach within a cluster
            vector<vector<pair<VertexID, vector<LabelSet>>>> RBI;
			// boundary nodes that can reach the current node within a cluster
            vector<vector<pair<VertexID, vector<LabelSet>>>> RRBI;
			// neighbouring boundary nodes that a node can reach outside a cluster
			vector<vector<pair<VertexID, LabelSet>>> ROBI;
			
			// each node -> a mapping from labelset to a list of BNs that the node can reach
			vector<unordered_map<LabelSet, unordered_set<VertexID>>> newRBI;

            // each cluster -> a list of clusters and minimal labelset
            vector<vector<pair<int, vector<LabelSet>>>> RRCI;
    };

#endif
