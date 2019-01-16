#include "NewIndex.h"

#include <queue>
#include <iostream>
#include <unordered_set>
#include <map>

using namespace newns;
using namespace indexns;


NewIndex::NewIndex(Graph* mg)
{
    this->graph = mg;
    this->isBlockedMode = false;
    this->indexType = IndexType::New;
    buildIndex();

    cout << "NewIndex-index size(byte)=" << getIndexSizeInBytes() << ", time(s)=" << getIndexConstructionTimeInSec() << endl;
};

NewIndex::~NewIndex()
{

};

unsigned long NewIndex::getIndexSizeInBytes()
{
    return getIndexSizeInBytesM();
};

void NewIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach)
{

};

void NewIndex::initializeLocalIndexes() {
    int N = graph->getNumberOfVertices();

    for (int i = 0; i < N; i++) {
        vector<pair<VertexID, vector<LabelSet>>> v1, v2;
        vector<pair<int, vector<LabelSet>>> v3;
        RBI.push_back(v1);
        RRBI.push_back(v2);
        RRCI.push_back(v3);
    }
}

void NewIndex::buildIndex()
{
    constStartTime = getCurrentTimeInMilliSec();

    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    this->didComplete = false;

    // construct tIn or cIn

    initializeIndex();
    initializeLocalIndexes();

    // first we divide graph G into a set of random clusters
    // each cluster -> a list of vertexID
    this->clusters = vector< vector < VertexID > >();
    // each vertex -> a cluster IDs
    this->vToCID = vector<int>(N, -1);
    this->isBoundaryNode = vector<bool>(N, false);
    // this->graph->tarjan(SCCs);
    this->graph->randomClustering(clusters, vToCID);
    
    // create a subgraph for each cluster containing only the right edges
    subGraphs = vector< Graph* >();
    vector< int > countPerCluster = vector< int >(clusters.size(), 0);

    // Initialize graph for each cluster
    for(int i = 0; i < clusters.size(); i++)
    {
        EdgeSet* es = new EdgeSet();
        DGraph* ng = new DGraph(es,clusters[i].size(),L);
        subGraphs.push_back( ng );
        vector<VertexID> v;
        boundaryNodesPerCluster.push_back(v);
    }

    // position of nodes in each cluster starting from 0 to cluster.size() - 1
    for (int i = 0; i < clusters.size(); i++) {
        vector<VertexID> cluster = clusters.at(i);
        for (int j = 0; j < cluster.size(); j++) {
            positionInC[cluster[j]] = j;
        }
    }

    for(int i = 0; i < N; i++)
    {
        int iCID = vToCID[i];
        //cout << "i=" << i << " ,vToSCCID[i]=" << iID << ", vToSubGraphID[i]=" << iSID << endl;

        SmallEdgeSet ses;
        graph->getOutNeighbours(i, ses);
        for(int j = 0; j < ses.size(); j++)
        {
            int jCID = vToCID[ses[j].first];
            // If in the same cluster, add edge to subgraph
            if( iCID == jCID )
            {
                // int jSID = vToSubGraphID[ ses[j].first ];
                //cout << "- ses[j].first =" << ses[j].first << ", jID=" << iID << ", jSID=" << jSID << endl;
                subGraphs[iCID]->addEdge( positionInC.at(i), positionInC.at(ses[j].first), labelSetToLabelID(ses[j].second) );
            }
            else
            {
                // DAGedges.push_back( make_pair(i, make_pair( ses[j].first, ses[j].second )) );
                if (isBoundaryNode[i] == false) {
                    boundaryNodesPerCluster.at(iCID).push_back(positionInC.at(i));    
                    isBoundaryNode[i] = true;
                }
                if (isBoundaryNode[ses[j].first] == false) {
                    boundaryNodesPerCluster.at(jCID).push_back(positionInC.at(ses[j].first));    
                    isBoundaryNode[ses[j].first] = true;    
                }
            }
        }
    }

    // subGraphs[4]->addEdge(1, 2, labelSetToLabelID(4));
    // cout << subGraphs[4]->toString() << endl;
    cout << "Boundary Nodes: ";
    for (int i = 0; i < N; i++) {
        if (isBoundaryNode[i]) {
            cout << i << " ";
        }
        cout << endl;
    }

    cout << "Step 1 (clustering): " << print_digits( getCurrentTimeInMilliSec()-constStartTime, 2 ) << endl;

    // build an index for each cluster first
    for(int i = 0; i < subGraphs.size(); i++) {
        labeledBFSPerCluster(i, subGraphs[i]);
        getRBI(i, subGraphs[i], clusters);
        getRRBI(i, subGraphs[i], clusters);
    }

    cout << "Step 2 (RBI & RRBI indices built): " << print_digits( getCurrentTimeInMilliSec()-constStartTime, 2 ) << endl;
    printRBI();
    printRRBI();

    initializeIndex();
    vector<bool> BNIndexed = vector<bool>(N, false);
    for (int i = 0; i < clusters.size(); i++) {
        labeledBFSAcrossClusters(i, clusters, BNIndexed);
        getRRCI(i);
    }

    printRRCI();

    cout << "Step 3 (RRCI indices built): " << print_digits( getCurrentTimeInMilliSec()-constStartTime, 2 ) << endl;

    this->didComplete = true;
    constEndTime = getCurrentTimeInMilliSec();
    totalConstTime = constEndTime - constStartTime;

};

void NewIndex::labeledBFSPerCluster(int cID, Graph* sG)
{

    int N = sG->getNumberOfVertices();
    cout << "buildIndex sG->N=" << N << " ,cID=" << cID << endl;
    vector<bool> indexed = vector<bool>(N, false);

    for (int i = 0; i < N; i++) {
        tIn.at(i).clear();
    }

    for(int i = 0; i < N; i++)
    {
        cout <<"At vertex: " << i << endl;
        labeledBFSPerVertex(cID, i, sG, indexed);
    }
    // if (cID == 4) {
    //     SmallEdgeSet ses;
    //     sG->getOutNeighbours(1, ses);
    //     for (int i = 0; i < ses.size(); i++) {
    //         cout<<ses[i].first << " " << ses[i].second << endl;
    //     }
    //     // for (int x = 0; x < N; x++) {
    //     //     int wCount = tIn.at(x).size();
    //     //     cout << x << ": ";
    //     //     for (int j = 0; j < wCount; j++) {
    //     //         int lsCount = tIn.at(x).at(j).second.size();
    //     //         cout << tIn.at(x).at(j).first << "(";
    //     //         for (int k = 0; k < lsCount; k++) {
    //     //             cout << tIn.at(x).at(j).second.at(k) << " ";
    //     //         }
    //     //         cout << ")";
    //     //     }
    //     //     cout << endl;
    //     // }
    // }
};

void NewIndex::labeledBFSPerVertex(int cID, VertexID v, Graph* sG, vector<bool>& indexed) {
    priority_queue< NewIndexBitEntry, vector<NewIndexBitEntry>, NewIndexPQBitEntries > q;
    NewIndexBitEntry t;
    t.x = v;
    t.ls = 0;
    t.dist = 0;

    q.push(t);
    while( q.empty() == false ) {
        // cout<<"Queue size: "<<q.size()<<endl;
        NewIndexBitEntry tr = q.top();
        VertexID v1 = tr.x;
        LabelSet ls1 = tr.ls;
        q.pop();
        // cout<<"Next element: "<<v1<<" "<<ls1<<endl;
        if (v != v1) {
            if (tryInsert(v, v1, ls1) == false) {
                continue;
            }
            if (indexed[v1] == true) {
                cout << "using index of " << v1 << endl;
                for (int i = 0; i < tIn[v1].size(); i++) {
                    pair<VertexID, vector<LabelSet>> p = tIn[v1].at(i);
                    vector<LabelSet> lss = p.second;
                    for (int j = 0; j < lss.size(); j++) {
                        cout << v << " " << p.first << " " << joinLabelSets(ls1, lss.at(j)) << endl;
                        tryInsert(v, p.first, joinLabelSets(ls1, lss.at(j)));
                    }
                }
                continue;
            }
        }

        SmallEdgeSet ses;
        sG->getOutNeighbours(v1, ses);

        for(int i = 0; i < ses.size(); i++)
        {
            VertexID v2 = ses[i].first;
            LabelSet ls2 = ses[i].second;
            // cout<<"Add to queue: "<<v2 <<" "<<ls2<<endl;
            LabelSet ls3 = joinLabelSets(ls1, ls2);

            if( v2 == v1 )
            {
                continue;
            }

            int dist = tr.dist;
            if( ls3 != ls1 || ls3 != ls2 )
            {
                dist += 1; // labels are added one by one
            }

            NewIndexBitEntry tr2;
            tr2.x = v2;
            tr2.ls = ls3;
            tr2.dist = dist;

            q.push( tr2 );
        }

    }
    indexed[v] = true;
}

bool NewIndex::tryInsert(VertexID s, VertexID v, LabelSet ls)
{
    if( s == v ) {
        return true;
    }
    
    bool b2 = tryInsertLabelSetToIndex(ls, s, v);
    //cout << "tryInsert: w=" << w << ",v=" << v << ",ls=" << labelSetToString(ls) << ",b2=" << b2 << endl;
    return b2;
}

void NewIndex::printRBI() {
    cout << "RBI" <<endl;
    int l1 = RBI.size();
    for (int i = 0; i < l1; i++) {
        cout << "Vertex " << i << ": ";
        int l2 = RBI[i].size();
        for (int j = 0; j < l2; j++) {
            cout << " ( " << RBI[i][j].first << " ";
            int l3 = RBI[i][j].second.size();
            for (int k = 0; k < l3; k++) {
                cout << (RBI[i][j].second)[k] << ", ";
            } 
            cout << ") ";
        }
        cout << endl;
    }
}

void NewIndex::printRRBI() {
    cout << "RRBI" << endl;
    int l1 = RRBI.size();
    for (int i = 0; i < l1; i++) {
        cout << "Vertex " << i << ": ";
        int l2 = RRBI[i].size();
        for (int j = 0; j < l2; j++) {
            cout << " ( " << RRBI[i][j].first << " ";
            int l3 = RRBI[i][j].second.size();
            for (int k = 0; k < l3; k++) {
                cout << (RRBI[i][j].second)[k] << ", ";
            } 
            cout << ") ";
        }
        cout << endl;
    }
}

void NewIndex::printRRCI() {
    cout << "RRCI" << endl;
    int l1 = RRCI.size();
    for (int i = 0; i < l1; i++) {
        cout << "Vertex " << i << ": ";
        int l2 = RRCI[i].size();
        for (int j = 0; j < l2; j++) {
            cout << " ( " << RRCI[i][j].first << " ";
            int l3 = RRCI[i][j].second.size();
            for (int k = 0; k < l3; k++) {
                cout << (RRCI[i][j].second)[k] << ", ";
            } 
            cout << ") ";
        }
        cout << endl;
    }
}

void NewIndex::printTIn(int cID) {
    cout << "tIn of " << cID << endl;
    int l1 = tIn.size();
    for (int i = 0; i < l1; i++) {
        cout << "Vertex " << i << ": ";
        int l2 = tIn[i].size();
        for (int j = 0; j < l2; j++) {
            cout << " ( " << tIn[i][j].first << " ";
            int l3 = tIn[i][j].second.size();
            for (int k = 0; k < l3; k++) {
                cout << (tIn[i][j].second)[k] << ", ";
            } 
            cout << ") ";
        }
        cout << endl;
    }
}

void NewIndex::getRBI(int cID, Graph* sG, vector<vector<VertexID>> clusters) {
    int N = sG->getNumberOfVertices();
    printTIn(cID);
    for (int i = 0; i < N; i++) {
        vector<pair<VertexID, vector<LabelSet>>> RBPerVertex;
        vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(i);
        for (int j = 0; j < closure.size(); j++) {
            if (isBoundaryNode[clusters.at(cID).at(closure.at(j).first)]) {
                RBPerVertex.push_back(closure.at(j));
            }
        }
        RBI.at(clusters.at(cID).at(i)) = RBPerVertex;
    }
}

void NewIndex::getRRBI(int cID, Graph* sG, vector<vector<VertexID>> clusters) {
    int N = sG->getNumberOfVertices();
    vector<VertexID> boundaryNodes = this->boundaryNodesPerCluster.at(cID);
    for (int i = 0; i < boundaryNodes.size(); i++) {
        // cout << i << " " << boundaryNodes[i];
        vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(boundaryNodes[i]);
        for (int j = 0; j < closure.size(); j++) {
            vector<LabelSet> RRBInstanceLabelSets = closure.at(j).second;
            VertexID vPositionInC = closure.at(j).first;
            RRBI.at(clusters.at(cID).at(vPositionInC)).push_back(make_pair(boundaryNodes.at(i), RRBInstanceLabelSets));
        }
    }
}

void NewIndex::labeledBFSAcrossClusters(int cID, vector<vector<VertexID>> clusters, vector<bool>& BNIndexed) {
    int N = graph->getNumberOfVertices();
    vector<VertexID> boundaryNodes = boundaryNodesPerCluster.at(cID);
    for (int i = 0; i < boundaryNodes.size(); i++) {
        int globalVID = clusters.at(cID).at(boundaryNodes.at(i));
        priority_queue< NewIndexBitEntry, vector<NewIndexBitEntry>, NewIndexPQBitEntries > q;
        NewIndexBitEntry t;
        t.x = globalVID;
        t.ls = 0;
        t.dist = 0;

        q.push(t);
        // cout << "Vertex :" << globalVID << endl;
        while (q.empty() == false) {
            NewIndexBitEntry tr = q.top();
            VertexID v1 = tr.x;
            LabelSet ls1 = tr.ls;
            q.pop();
            // cout << v1 << " " << ls1 << endl;
            if (v1 != globalVID) {
                if (tryInsert(globalVID, v1, ls1) == false) {
                    continue;
                }
                if (BNIndexed[v1] == true) {
                    // cout << "indexed" << v1 << " " << tIn[v1].size() << endl;
                    for (int k = 0; k < tIn[v1].size(); k++) {
                        pair<VertexID, vector<LabelSet>> p = tIn[v1].at(k);
                        vector<LabelSet> lss = p.second;
                        for (int j = 0; j < lss.size(); j++) {
                            // cout << p.first << " " << lss.at(j) << endl;
                            tryInsert(globalVID, p.first, joinLabelSets(ls1, lss.at(j)));
                        }
                        
                    }
                    continue;
                }
            }

            SmallEdgeSet ses;
            graph->getOutNeighbours(v1, ses);
            for (int i = 0; i < ses.size(); i++) {
                VertexID v2 = ses[i].first;
                LabelSet ls2 = ses[i].second;
                LabelSet ls3 = joinLabelSets(ls1, ls2);

                if (vToCID.at(v2) != vToCID.at(v1)) {
                    if (isBoundaryNode[v2] == true) {
                        int dist = tr.dist;
                        if( ls3 != ls1 || ls3 != ls2 )
                        {
                            dist += 1; // labels are added one by one
                        }

                        NewIndexBitEntry tr2;
                        tr2.x = v2;
                        tr2.ls = ls3;
                        tr2.dist = dist;

                        q.push( tr2 );
                    }
                } 
            }

            // No RBI paths should be added within source cluster
            if (vToCID.at(v1) == cID) {
                continue;
            }
            vector<pair<VertexID, vector<LabelSet>>> RBIv1 = RBI.at(v1);
            for (int i = 0; i < RBIv1.size(); i++) {
                VertexID v2 = clusters.at(vToCID.at(v1)).at(RBIv1.at(i).first);
                if (isBoundaryNode[v2] == true) {
                    vector<LabelSet> lss = RBIv1.at(i).second;
                    for (int j = 0; j < lss.size(); j++) {
                        LabelSet ls2 = lss.at(j);
                        LabelSet ls3 = joinLabelSets(ls1, ls2);
                        int dist = tr.dist;
                        if( ls3 != ls1 || ls3 != ls2 )
                        {
                            dist += 1; // labels are added one by one
                        }

                        NewIndexBitEntry tr2;
                        tr2.x = v2;
                        tr2.ls = ls3;
                        tr2.dist = dist;

                        q.push( tr2 );
                    }
                }
            }
        }
        BNIndexed[globalVID] = true;
    }
}

void NewIndex::getRRCI(int cID) {
    vector<VertexID> boundaryNodes = boundaryNodesPerCluster.at(cID);
    for (int i = 0; i < boundaryNodes.size(); i++) {
        int globalVID = clusters.at(cID).at(boundaryNodes.at(i));
        vector<pair<VertexID, vector<LabelSet>>> closure = tIn[globalVID];
        for (int j = 0; j < closure.size(); j++) {
            int v1 = closure.at(j).first;
            vector<LabelSet> lss = closure.at(j).second;
            int cIDj = vToCID[v1];
            RRCI.at(cIDj).push_back(make_pair(cID, lss));
        }
    }
}

bool NewIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    //cout << "NewIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();
    bool b = queryShell(source, target, ls);
    queryEndTime = getCurrentTimeInMilliSec();
    //cout << "NewIndex::query answer =" << b << endl;
    return b;
}

bool NewIndex::queryShell(VertexID source, VertexID target, LabelSet ls)
{
    if(source == target)
        return true;
    if( ls == 0 )
        return false;

    if (vToCID[source] == vToCID[target]) {
        unordered_set<int> possibleBoundaryNodes;
        vector<pair<VertexID, vector<LabelSet>>> RBIs = RBI.at(source);
        for (int i = 0; i < RBIs.size(); i++) {
            vector<LabelSet> BNLabelSets = RBIs.at(i).second;
            for (int j = 0; j < BNLabelSets.size(); j++) {
                if (joinLabelSets(ls, BNLabelSets.at(j)) == ls) {
                    possibleBoundaryNodes.insert((int)RBIs.at(i).first);
                    break;
                }
            } 
        }
        vector<pair<VertexID, vector<LabelSet>>> RRBIt = RRBI.at(target);
        for (int i = 0; i < RRBIt.size(); i++) {
            if (possibleBoundaryNodes.find(RRBIt.at(i).first) != possibleBoundaryNodes.end()) {
                vector<LabelSet> BNLabelSets = RRBIt.at(i).second;
                for (int j = 0; j < BNLabelSets.size(); j++) {
                    if (joinLabelSets(ls, BNLabelSets.at(j)) == ls) {
                        return true;
                    }
                } 
            }
        
        queue< VertexID > q;
        q.push(source);

        int N = graph->getNumberOfVertices();
        dynamic_bitset<> marked = dynamic_bitset<>(N);

        while( q.empty() == false )
        {
            VertexID x = q.front();
            q.pop();

            if( x == target )
                return true;

            if( marked[x] == 1 ){
                continue;
            }
            marked[x] = 1;

            SmallEdgeSet ses;
            graph->getOutNeighbours(x, ses);
            for(int i = 0; i < ses.size(); i++) {
                if( isLabelSubset(ses[i].second, ls) == true ) {
                    q.push( ses[i].first );
                }
            }
        }

        return false;  
        }
    }

    unordered_set<int> X;
    vector<pair<int, vector<LabelSet>>> RRCIt = RRCI.at(vToCID[target]);
    for (int i = 0; i < RRCIt.size(); i++) {
        // cout << RRCIt.at(i).first << " ";
        vector<LabelSet> lss = RRCIt.at(i).second;
        for (int j = 0; j < lss.size(); j++) {
            // cout << lss.at(j) << " " << (joinLabelSets(ls, lss.at(j)) == ls) << " xxx";
            if (joinLabelSets(ls, lss.at(j)) == ls) {
                X.insert(RRCIt.at(i).first);
                break;
            }
        } 
        // cout << endl;
    }
    X.insert(vToCID[target]);

    if (X.find(vToCID[source]) == X.end()) {
        return false;
    }

    unordered_set<int> BS;
    vector<pair<VertexID, vector<LabelSet>>> RBIs = RBI.at(source);
    for (int i = 0; i < RBIs.size(); i++) {
        vector<LabelSet> BNLabelSets = RBIs.at(i).second;
        for (int j = 0; j < BNLabelSets.size(); j++) {
            if (joinLabelSets(ls, BNLabelSets.at(j)) == ls) {
                BS.insert((int)RBIs.at(i).first);
                break;
            }
        } 
    }

    if (isBoundaryNode[source]) {
        BS.insert(source);
    }

    if (BS.empty()) {
        return false;
    }

    unordered_set<int> BT;
    vector<pair<VertexID, vector<LabelSet>>> RRBIt = RRBI.at(target);
    for (int i = 0; i < RRBIt.size(); i++) {
        vector<LabelSet> BNLabelSets = RRBIt.at(i).second;
        for (int j = 0; j < BNLabelSets.size(); j++) {
            if (joinLabelSets(ls, BNLabelSets.at(j)) == ls) {
                BT.insert((int)RRBIt.at(i).first);
                break;
            }
        } 
    }

    if (isBoundaryNode[target]) {
        BT.insert(target);
    }

    if (BT.empty()) {
        return false;
    }

    // cout << "BT Size: " << BT.size() << endl;

    int N = graph->getNumberOfVertices();
    map<VertexID, int> BNToID;
    vector<bool> visited;
    int count = 0;
    for (int i = 0; i < N; i++) {
        if (isBoundaryNode[i]) {
            visited.push_back(false);
            BNToID[i] = count;
            count ++;
        }
    }
    // for (int i = 0; i < visited.size(); i++) {
    //     cout << visited[i] << " ";
    // }
    // cout << endl;

    while (!BS.empty()) {
        unordered_set<int> BS1;
        for (unordered_set<int>::iterator i = BS.begin(); i != BS.end(); i++) {
            // cout << "v in BS: " << *i << endl;
            vector<pair<VertexID, vector<LabelSet>>> RBIi = RBI.at(*i);
            // cout << "RBIi size: " << RBIi.size() << endl; 
            for (int j = 0; j < RBIi.size(); j++) {
                int globalVID = clusters.at(vToCID[*i]).at(RBIi.at(j).first);
                // cout << "globalVID: " << globalVID << endl;
                if (visited[BNToID[globalVID]] == false) { 
                    vector<LabelSet> lss = RBIi.at(j).second;
                    for (int k = 0; k < lss.size(); k++) {
                        if (joinLabelSets(lss.at(k), ls) == ls) {
                            BS1.insert(globalVID);
                            // cout << "Inserted!" << endl;
                            visited[BNToID[globalVID]] = true;
                            break;
                        }
                    }
                }
            }

            SmallEdgeSet ses;
            graph->getOutNeighbours(*i, ses);
            for (int j = 0; j < ses.size(); j++) {
                VertexID v2 = ses[j].first;
                LabelSet ls2 = ses[j].second;
                // cout << "v2: " << v2 << " ls2 :" << ls2 << endl;
                // cout << "Cv: " << vToCID[*i] << " Cv2: "<< vToCID[v2] << endl;
                // cout << "ls: " << joinLabelSets(ls2, ls) << endl;
                // cout << "visited? " << visited[BNToID[v2]] << endl;
                // for (unordered_set<int>::iterator k = X.begin(); k != X.end(); k++) {
                //     cout << *k << "__";
                // }
                // cout << endl;
                // cout << (X.find(vToCID[v2]) != X.end()) << endl;
                // cout << (BT.find(v2) != BT.end()) << endl;
                if ((vToCID[*i] != vToCID[v2]) && (joinLabelSets(ls2, ls) == ls)) {
                    if ((visited[BNToID[v2]] == false) && (X.find(vToCID[v2]) != X.end())) {
                        if (BT.find(v2) != BT.end()) {
                            return true;
                        }
                        BS1.insert((int)v2);
                        visited[BNToID[v2]] = true;
                    }
                } 
            }
        }

        BS = BS1;
    }

    return false;
};
