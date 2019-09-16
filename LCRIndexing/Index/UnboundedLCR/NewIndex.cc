#include "NewIndex.h"

#include <queue>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>

using namespace newns;
using namespace indexns;


NewIndex::NewIndex(Graph* mg) {
    this->graph = mg;
    this->isBlockedMode = false;
    this->indexType = IndexType::New;
    buildIndex();

    cout << "NewIndex-index size(byte)=" << getIndexSizeInBytes() << ", time(s)=" << getIndexConstructionTimeInSec() << endl;
};

NewIndex::~NewIndex() {

};

unsigned long NewIndex::getIndexSizeInBytes() {
    unsigned long size = 0;
    int N = graph->getNumberOfVertices();
    int L = sizeof(LabelSet);
    int emptyVectorSize = 3 * sizeof(int); // estimated 3 times 32-bits
    unsigned long RBISize = 0;
    unsigned long RRBISize = 0;
    unsigned long RRCISize = 0;

    cout << "RBI RRBI" << endl;
    for (int i = 0; i < N; i++) {
        vector<pair<VertexID, vector<LabelSet>>> RBIi = RBI[i];
        for (int j = 0; j < RBIi.size(); j++) {
            vector<LabelSet> lss = RBIi[j].second;
            if (lss.size() > 0) {
                RBISize += sizeof(VertexID) + lss.size() * L + emptyVectorSize;
            }
        }

        vector<pair<VertexID, vector<LabelSet>>> RRBIi = RRBI[i];
        for (int j = 0; j < RRBIi.size(); j++) {
            vector<LabelSet> lss = RRBIi[j].second;
            if (lss.size() > 0) {
                RRBISize += sizeof(VertexID) + lss.size() * L + emptyVectorSize;
            }
        }
    }

    cout << "RRCI" << endl;
    for (int i = 0; i < clusters.size(); i++) {
        vector<pair<int, vector<LabelSet>>> RRCIi = RRCI[i];
        for (int j = 0; j < RRCIi.size(); j++) {
            vector<LabelSet> lss = RRCIi[j].second;
            if (lss.size() > 0) {
                RRCISize += sizeof(VertexID) + lss.size() * L + emptyVectorSize;
            }
        }
    }

    cout << "RBI Size: " << RBISize << endl;
    cout << "RRBI Size: " << RRBISize << endl;
    cout << "RRCI Size: " << RRCISize << endl;

    size += RBISize + RRBISize + RRCISize + graph->getGraphSizeInBytes();

    return size;
};

void NewIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach) {

};

void NewIndex::initializeLocalIndexes() {
    int N = graph->getNumberOfVertices();

    for (int i = 0; i != N; i++) {
		vector<pair<VertexID, vector<LabelSet>>> v1, v2;
		vector<pair<VertexID, LabelSet>> v3;
        vector<pair<int, vector<LabelSet>>> v4;
		unordered_map<LabelSet, unordered_set<VertexID>> m1;
        RBI.push_back(v1);
        RRBI.push_back(v2);
		ROBI.push_back(v3);
        RRCI.push_back(v4);
		newRBI.push_back(m1);
    }
}

void NewIndex::buildIndex() {
    constStartTime = getCurrentTimeInMilliSec();

    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    this->didComplete = false;

    // construct tIn or cIn

    initializeIndex();
    initializeLocalIndexes();

    // first we divide graph G into a set of random clusters
    // each cluster -> a list of vertexID
    this->clusters;
    // each vertex -> a cluster IDs
    this->vToCID = vector<int>(N, -1);
    this->isBoundaryNode = vector<bool>(N, false);
    // this->graph->newClustering(clusters, vToCID, 15, N / 10, N / 30);
    this->graph->minBoundaryNodesClustering(clusters, vToCID, max(N / 20, 10));


     for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
         vector<VertexID> cluster = clusters[i];
         cout << "Cluster " << i << ": ";
         for (int j = 0, sizeJ = cluster.size(); j != sizeJ; ++j) {
             cout << cluster[j] << " ";
         }
         cout << endl;
     }

    // this->graph->randomClustering(clusters, vToCID);

    // create a subgraph for each cluster containing only the right edges
    subGraphs = vector< Graph* >();
    vector<int> countPerCluster = vector<int>(clusters.size(), 0);

    // Initialize graph for each cluster
    for (int i = 0; i < clusters.size(); i++) {
        EdgeSet* es = new EdgeSet();
        DGraph* ng = new DGraph(es, clusters[i].size(), L);
        subGraphs.push_back( ng );
        vector<VertexID> v;
        boundaryNodesPerCluster.push_back(v);
    }

    // position of nodes in each cluster starting from 0 to cluster.size() - 1
    for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
        vector<VertexID> cluster = clusters.at(i);
        for (int j = 0, sizeJ = cluster.size(); j != sizeJ; ++j) {
            positionInC[cluster[j]] = j;
        }
    }

    for (int i = 0; i < N; i++) {
        int iCID = vToCID[i];
        //cout << "i=" << i << " ,vToSCCID[i]=" << iID << ", vToSubGraphID[i]=" << iSID << endl;

        SmallEdgeSet ses;
        graph->getOutNeighbours(i, ses);
        for (int j = 0; j < ses.size(); j++) {
            int jCID = vToCID[ses[j].first];
            // If in the same cluster, add edge to subgraph
            if ( iCID == jCID ) {
                // int jSID = vToSubGraphID[ ses[j].first ];
                //cout << "- ses[j].first =" << ses[j].first << ", jID=" << iID << ", jSID=" << jSID << endl;
                subGraphs[iCID]->addEdge( positionInC.at(i), positionInC.at(ses[j].first), labelSetToLabelID(ses[j].second) );
            } else {
				ROBI[i].push_back(ses[j]);
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

    int BNCount = 0;
    cout << "Boundary Nodes: ";
    for (int i = 0; i < N; i++) {
        if (isBoundaryNode[i]) {
            boundaryNodes.push_back(i);
        }
    }
    cout << boundaryNodes.size() << endl;


    cout << "Step 1 (clustering): " << print_digits( getCurrentTimeInMilliSec() - constStartTime, 6 ) << endl;

    // build an index for each cluster first
    for (int i = 0; i < subGraphs.size(); i++) {
        labeledBFSPerCluster(i, subGraphs[i]);
        getRBI(i, subGraphs[i], clusters);
        getRRBI(i, subGraphs[i], clusters);
    }

    cout << "Step 2 (RBI & RRBI indices built): " << print_digits( getCurrentTimeInMilliSec() - constStartTime, 6 ) << endl;
    // printRBI();
	printNewRBI();
    // printRRBI();

    initializeIndex();
    vector<bool> BNIndexed = vector<bool>(N, false);
    for (int i = 0; i < clusters.size(); i++) {
        labeledBFSAcrossClusters(i, clusters, BNIndexed);
        getRRCI(i);
    }

    // printRRCI();

    cout << "Step 3 (RRCI indices built): " << print_digits( getCurrentTimeInMilliSec() - constStartTime, 6 ) << endl;

    this->didComplete = true;
    constEndTime = getCurrentTimeInMilliSec();
    totalConstTime = constEndTime - constStartTime;

};

void NewIndex::labeledBFSPerCluster(int cID, Graph* sG) {

    int N = sG->getNumberOfVertices();
    // cout << "buildIndex sG->N=" << N << " ,cID=" << cID << endl;
    vector<bool> indexed = vector<bool>(N, false);

    for (int i = 0; i < N; i++) {
        tIn.at(i).clear();
    }

    for (int i = 0; i < N; i++) {
        // cout <<"At vertex: " << i << endl;
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
    while ( q.empty() == false ) {
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
                // cout << "using index of " << v1 << endl;
                for (int i = 0; i < tIn[v1].size(); i++) {
                    pair<VertexID, vector<LabelSet>> p = tIn[v1].at(i);
                    vector<LabelSet> lss = p.second;
                    for (int j = 0; j < lss.size(); j++) {
                        // cout << v << " " << p.first << " " << joinLabelSets(ls1, lss.at(j)) << endl;
                        tryInsert(v, p.first, joinLabelSets(ls1, lss.at(j)));
                    }
                }
                continue;
            }
        }

        SmallEdgeSet ses;
        sG->getOutNeighbours(v1, ses);

        for (int i = 0; i < ses.size(); i++) {
            VertexID v2 = ses[i].first;
            LabelSet ls2 = ses[i].second;
            // cout<<"Add to queue: "<<v2 <<" "<<ls2<<endl;
            LabelSet ls3 = joinLabelSets(ls1, ls2);

            if ( v2 == v1 ) {
                continue;
            }

            int dist = tr.dist;
            if ( ls3 != ls1 || ls3 != ls2 ) {
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

bool NewIndex::tryInsert(VertexID s, VertexID v, LabelSet ls) {
    if ( s == v ) {
        return true;
    }

    bool b2 = tryInsertLabelSetToIndex(ls, s, v);
    //cout << "tryInsert: w=" << w << ",v=" << v << ",ls=" << labelSetToString(ls) << ",b2=" << b2 << endl;
    return b2;
}

void NewIndex::printNewRBI() {
	cout << "New RBI" << endl;
	for (int i = 0; i < newRBI.size(); i++) {
		cout << "Vertex " << i << ": " << endl;
		for (auto& x : newRBI[i]) {
			cout << x.first << ": "; 
			unordered_set<VertexID> vs = x.second;
			for (auto& y : vs) {
				cout << y << " ";
			}
			cout << endl;
		}
	}
}

void NewIndex::printRBI() {
    cout << "RBI" << endl;
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
        int l2 = RRCI[i].size();
        if (l2 == 0) {
            break;
        }
        cout << "Cluster " << i << ": ";
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
    // printTIn(cID);
    for (int i = 0; i < N; i++) {
        vector<pair<VertexID, vector<LabelSet>>> RBPerVertex;
        vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(i);
        for (int j = 0; j < closure.size(); j++) {
            VertexID globalVID = clusters.at(cID).at(closure.at(j).first);
            if (isBoundaryNode[globalVID]) {
                RBPerVertex.push_back(make_pair(globalVID, closure.at(j).second)); 
            }
        }
        RBI.at(clusters.at(cID).at(i)) = RBPerVertex;
    }
	int L = sG->getNumberOfLabels();
	vector<unordered_set<LabelSet>> labelSetBuckets;
	// LabelSet with number of labels = X goes to Xth bucket
	// bucket[0] should be empty
	for (int j = 0; j <= L; j++) {
		unordered_set<LabelSet> s1;
		labelSetBuckets.push_back(s1);
	}
	
	cout << "L in cluster " << cID << ": " << L << endl;

	//vector<VertexID> cluster = clusters[cID];
	//for (int i = 0; i != N; i++) {
	//	vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(i);
	//	VertexID globalVID1 = cluster[i];
	//	unordered_map<LabelSet, unordered_set<VertexID>> lmap = newRBI[globalVID1];

	//	// organization
	//	for (int j = 0, sizeJ = closure.size(); j != sizeJ; j++) {
	//		VertexID globalVID2 = cluster[closure[j].first];
	//		LabelSets lss = closure[j].second;
	//		for (int k = 0, sizeK = lss.size(); k != sizeK; k++) {
	//			// propagation preparation
	//			labelSetBuckets[getNumberOfLabelsInLabelSet(lss[k])].insert(lss[k]);
	//			// end of preparation, continue organization
	//			if (lmap.count(lss[k])) {
	//				lmap[lss[k]].insert(globalVID2);
	//			}
	//			else {
	//				unordered_set<VertexID> s1({ globalVID2 });
	//				lmap.insert({ lss[k], s1 });
	//			}
	//		}
	//	}

	//	// propagation
	//	for (int j = 1; j <= L; j++) {
	//		unordered_set<VertexID> labelSetBucket = labelSetBuckets[j];
	//		for (const auto& ls : labelSetBucket) {
	//			vector<LabelID> labels;
	//			getLabelIDsFromLabelSet(ls, &labels);
	//			for (const auto& label : labels) {
	//				LabelSet reducedLs = ls - label;
	//				if (lmap.count(reducedLs)) {
	//					lmap[ls].insert(lmap[reducedLs]);
	//				}
	//			}
	//		}
	//	}

	vector<VertexID> cluster = clusters[cID];
	for (int i = 0; i != N; i++) {
		vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(i);
		VertexID globalVID1 = cluster[i];
		unordered_map<LabelSet, unordered_set<VertexID>> lmap = newRBI[globalVID1];

		// organization
		for (int j = 0, sizeJ = closure.size(); j != sizeJ; j++) {
			VertexID globalVID2 = cluster[closure[j].first];

			LabelSets lss = closure[j].second;
			for (int k = 0, sizeK = lss.size(); k != sizeK; k++) {
				// propagation preparation
				labelSetBuckets[getNumberOfLabelsInLabelSet(lss[k])].insert(lss[k]);
				// end of preparation, continue organization
				if (lmap.count(lss[k])) {
					lmap[lss[k]].insert(globalVID2);
				}
				else {
					unordered_set<VertexID> s1({ globalVID2 });
					lmap.insert({ lss[k], s1 });
				}
			}
		}

		cout << "LabelSetBucket: " << endl;
		for (int j = 0; j <= L; j++) {
			unordered_set<VertexID> labelSetBucket = labelSetBuckets[j];
			cout << j << ": ";
			for (auto& ls : labelSetBucket) {
				cout << ls << " ";
			}
			cout << endl;
		}

		// assume that the L in the cluster = the L in the graph
		// use a topdown approach to from full labelset, to get all possible labelsets
		// this can be extracted to be repeatedly used
		LabelSet lss = (1 << L) - 1;

		// insert the full-label into the lmap and bucket if it does not exist
		if (!lmap.count(lss)) {
			unordered_set<VertexID> s1;
			lmap.insert({ lss, s1 });
			labelSetBuckets[L].insert(lss);
		}
		
		// for each iteration, take the labelsets from the bucket 1 level above
		// reduce 1 label at a time to obtain new "less" labelset
		// a full labelset should have been built up
		for (int j = L - 1; j != 0; j--) {
			unordered_set<VertexID> labelSetBucket = labelSetBuckets[j + 1];
			for (const auto& ls : labelSetBucket) {
				vector<LabelID> labels;
				getLabelIDsFromLabelSet(ls, labels);
				for (const auto& label : labels) {
					LabelSet reducedLs = ls - (1 << label);
					if (!lmap.count(reducedLs)) {
						unordered_set<VertexID> s1;
						lmap.insert({ reducedLs, s1 });
						labelSetBuckets[j].insert(reducedLs);
					}
				}
			}
		}

		cout << "LabelSetBucket after population: " << endl;
		for (int j = 0; j <= L; j++) {
			unordered_set<VertexID> labelSetBucket = labelSetBuckets[i];
			cout << j << ": ";
			for (auto& ls : labelSetBucket) {
				cout << ls << " ";
			}
			cout << endl;
		}

		// propagation
		// for each labelset ls 
		// corresponding vertices = { vertices of ls - 1 label }
		for (int j = 2; j <= L; j++) {
			unordered_set<VertexID> labelSetBucket = labelSetBuckets[j];
			for (const auto& ls : labelSetBucket) {
				vector<LabelID> labels;
				getLabelIDsFromLabelSet(ls, labels);
				for (const auto& label : labels) {
					LabelSet reducedLs = ls - label;
					lmap[ls].insert(lmap[reducedLs].begin(), lmap[reducedLs].end());
				}
			}
		}
	}	
}

void NewIndex::getRRBI(int cID, Graph* sG, vector<vector<VertexID>> clusters) {
    int N = sG->getNumberOfVertices();
    vector<VertexID> cluster = clusters.at(cID);
    vector<VertexID> boundaryNodes = this->boundaryNodesPerCluster.at(cID);

    for (int i = 0; i < boundaryNodes.size(); i++) {
        VertexID globalBoundaryVID = cluster.at(boundaryNodes[i]);
        vector<pair<VertexID, vector<LabelSet>>> closure = tIn.at(boundaryNodes[i]);
        for (int j = 0; j < closure.size(); j++) {
            vector<LabelSet> RRBInstanceLabelSets = closure.at(j).second;
            VertexID globalVID = cluster.at(closure.at(j).first);
            RRBI.at(globalVID).push_back(make_pair(globalBoundaryVID, RRBInstanceLabelSets));
        }
    }
}

void NewIndex::labeledBFSAcrossClusters(int cID, vector<vector<VertexID>> clusters, vector<bool>& BNIndexed) {
    // cout << "Cluster ID: " << cID << endl;
    int N = graph->getNumberOfVertices();
    vector<VertexID> boundaryNodes = boundaryNodesPerCluster.at(cID);
    vector<VertexID> cluster = clusters.at(cID);
    for (int i = 0; i < boundaryNodes.size(); i++) {
        int globalVID = cluster.at(boundaryNodes.at(i));
        // cout << "Visiting boundary nodes " << globalVID << endl;
        priority_queue< NewIndexBitEntry, vector<NewIndexBitEntry>, NewIndexPQBitEntries > q;
        NewIndexBitEntry t;
        t.x = globalVID;
        t.ls = 0;
        t.dist = 0;

        q.push(t);

        while (!q.empty()) {
            NewIndexBitEntry tr = q.top();
            VertexID v1 = tr.x;
            LabelSet ls1 = tr.ls;
            q.pop();
            // cout << "Viewing queue element: " << v1 << ", " << ls1 << endl;
            if (v1 != globalVID) {
                if (!tryInsert(globalVID, v1, ls1)) {
                    continue;
                }
                if (BNIndexed[v1]) {
                    // cout << "Use indexed " << v1 << " with " << tIn[v1].size() << " index entries" << endl;
                    for (int k = 0; k < tIn[v1].size(); k++) {
                        pair<VertexID, vector<LabelSet>> p = tIn[v1].at(k);
                        vector<LabelSet> lss = p.second;
                        for (int j = 0; j < lss.size(); j++) {
                            // cout << "Trying to add " << p.first << ", " << lss.at(j) << " --> " << joinLabelSets(ls1, lss.at(j)) << endl;
                            tryInsert(globalVID, p.first, joinLabelSets(ls1, lss.at(j)));
                        }

                    }
                    continue;
                }
            }

            // cout << "BN from another cluster" << endl;

            SmallEdgeSet ses;
            graph->getOutNeighbours(v1, ses);
            for (int i = 0; i < ses.size(); i++) {
                VertexID v2 = ses[i].first;
                LabelSet ls2 = ses[i].second;
                LabelSet ls3 = joinLabelSets(ls1, ls2);

                if (vToCID.at(v2) != vToCID.at(v1)) {
                    if (isBoundaryNode[v2] == true) {
                        int dist = tr.dist;
                        if ( ls3 != ls1 || ls3 != ls2 ) {
                            dist += 1; // labels are added one by one
                        }

                        NewIndexBitEntry tr2;
                        tr2.x = v2;
                        tr2.ls = ls3;
                        tr2.dist = dist;
                        // cout << "Add to queue: " << v2 << " , " << ls3 << endl;
                        q.push( tr2 );
                    }
                }
            }

            // // No RBI paths should be added within source cluster
            // if (vToCID.at(v1) == cID) {
            //     continue;
            // }

            // cout << "BN within the same cluster" << endl;

            vector<pair<VertexID, vector<LabelSet>>> RBIv1 = RBI.at(v1);
            for (int i = 0; i < RBIv1.size(); i++) {
                VertexID v2 = RBIv1.at(i).first;
                if (isBoundaryNode[v2] == true) {
                    vector<LabelSet> lss = RBIv1.at(i).second;
                    for (int j = 0; j < lss.size(); j++) {
                        LabelSet ls2 = lss.at(j);
                        LabelSet ls3 = joinLabelSets(ls1, ls2);
                        int dist = tr.dist;
                        if ( ls3 != ls1 || ls3 != ls2 ) {
                            dist += 1; // labels are added one by one
                        }

                        NewIndexBitEntry tr2;
                        tr2.x = v2;
                        tr2.ls = ls3;
                        tr2.dist = dist;
                        // cout << "Add to queue: " << v2 << " , " << ls3 << endl;
                        q.push( tr2 );
                    }
                }
            }
        }
        BNIndexed[globalVID] = true;

        // cout << endl;
        // cout << "Index of " << globalVID << endl;
        // for (int k = 0; k < tIn[globalVID].size(); k++) {
        //     pair<VertexID, vector<LabelSet>> p = tIn[globalVID].at(k);
        //     vector<LabelSet> lss = p.second;
        //     for (int j = 0; j < lss.size(); j++) {
        //         cout << p.first << ", " << lss.at(j) << endl;
        //     }
        // }
        // cout << endl;
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
            if (cIDj != cID) {
                RRCI.at(cIDj).push_back(make_pair(cID, lss));
            }
        }
    }
}

bool NewIndex::query(VertexID source, VertexID target, LabelSet ls) {
    //cout << "NewIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();
    bool b = queryShell(source, target, ls);
    queryEndTime = getCurrentTimeInMilliSec();
    //cout << "NewIndex::query answer =" << b << endl;
    return b;
}

bool NewIndex::queryShell(VertexID source, VertexID target, LabelSet ls) {
    double queryStartTime = getCurrentTimeInMilliSec();
    if (source == target)
        return true;
    if (ls == 0)
        return false;

    int sourceCluster = vToCID[source];

    // same cluster
    if (sourceCluster == vToCID[target]) {
        // cout << "Same cluster " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 );
        unordered_set<int> possibleBoundaryNodes;
        vector<pair<VertexID, vector<LabelSet>>> RBIs = RBI[source];
        for (unsigned int i = 0, sizeI = RBIs.size(); i != sizeI; ++i) {
            vector<LabelSet> BNLabelSets = RBIs[i].second;
            for (unsigned int j = 0, sizeJ = BNLabelSets.size(); j != sizeJ; ++j) {
                if (isLabelSubset(BNLabelSets[j], ls)) {
                    possibleBoundaryNodes.insert((int)RBIs[i].first);
                    break;
                }
            }
        }

        // cout << "Possible boundary nodes: ";
        // for (unordered_set<int>::iterator i = possibleBoundaryNodes.begin(); i != possibleBoundaryNodes.end(); i++) {
        //     cout << *i << ", ";
        // }
        // cout << endl;

        vector<pair<VertexID, vector<LabelSet>>> RRBIt = RRBI[target];
        for (unsigned int i = 0, sizeI = RRBIt.size(); i != sizeI; ++i) {
            if (possibleBoundaryNodes.count(RRBIt[i].first) > 0) {
                vector<LabelSet> BNLabelSets = RRBIt[i].second;
                for (unsigned int j = 0, sizeJ = BNLabelSets.size(); j != sizeJ; ++j) {
                    if (isLabelSubset(BNLabelSets[j], ls)) {
                        return true;
                    }
                }
            }
        }

        queue< VertexID > q;
        q.push(source);

        // cout << "BFS within cluster " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 );
        // Use a BFS within cluster
        int N = graph->getNumberOfVertices();
        dynamic_bitset<> marked = dynamic_bitset<>(N);

        while ( q.empty() == false ) {
            VertexID x = q.front();
            q.pop();

            // cout << "Viewing queued element: " << x << endl;

            if ( x == target )
                return true;

            if ( marked[x] == 1 ) {
                continue;
            }
            marked[x] = 1;

            SmallEdgeSet ses;
            graph->getOutNeighbours(x, ses);
            for (unsigned int i = 0, sizeI = ses.size(); i != sizeI; ++i) {
				// 3 conditions: in the same cluster, not a boundary node (or a path can be confirmed in the front part), compatible labelset
				// same cluster condition can be merged into not a boundary node condition
                if (!isBoundaryNode[ses[i].first] && isLabelSubset(ses[i].second, ls)) {
                    // cout << "Add to queue: " << ses[i].first << endl;
                    q.push(ses[i].first);
                }
            }
        }
    }

    // Found out clusters that can reach target cluster
    // cout << "X " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 );
    unordered_set<VertexID> X;
    vector<pair<int, vector<LabelSet>>> RRCIt = RRCI.at(vToCID[target]);
    for (unsigned int i = 0, sizeI = RRCIt.size(); i != sizeI; ++i) {
        vector<LabelSet> lss = RRCIt[i].second;
        for (unsigned int j = 0, sizeJ = lss.size(); j != sizeJ; ++j) {
            if (isLabelSubset(lss.at(j), ls)) {
                X.insert(RRCIt.at(i).first);
                break;
            }
        }
    }

	X.insert(vToCID[target]);

    // If the source and target clusters can't be linked by another common cluster, return false
    if (X.count(vToCID[source]) == 0) {
        return false;
    }

    // Get boundary nodes that source can reach
	// BS_same stores nodes obtained by RBI in the same cluster
	// BS_diff stores neighbouring BNs from a different cluster
	// cout << "BS_same " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 );
    unordered_set<VertexID> BS_same;
	unordered_set<VertexID> BS_diff;

    /*vector<pair<VertexID, vector<LabelSet>>> RBIs = RBI.at(source);
    for (unsigned int i = 0, sizeI = RBIs.size(); i != sizeI; ++i) {
        vector<LabelSet> BNLabelSets = RBIs.at(i).second;
        for (unsigned int j = 0, sizeJ = BNLabelSets.size(); j != sizeJ; j++) {
            if (isLabelSubset(BNLabelSets.at(j), ls)) {
                BS_same.insert((int)RBIs.at(i).first);
                break;
            }
        }
    }*/

	BS_same = newRBI[source][ls];

    if (isBoundaryNode[source]) {
        BS_same.insert(source);
    }

    if (BS_same.empty()) {
        return false;
    }

    // cout << "BS_same: ";
    // for (unordered_set<int>::iterator i = BS_same.begin(); i != BS_same.end(); i++) {
    //     cout << *i << ", ";
    // }
    // cout << endl;

    // Get boundary nodes that can reach target
    // cout << "BT " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 );
    unordered_set<int> BT;
    vector<pair<VertexID, vector<LabelSet>>> RRBIt = RRBI.at(target);
    for (unsigned int i = 0, sizeI = RRBIt.size(); i != sizeI; ++i) {
        vector<LabelSet> BNLabelSets = RRBIt.at(i).second;
        for (unsigned int j = 0, sizeJ = BNLabelSets.size(); j != sizeJ; j++) {
            if (isLabelSubset(BNLabelSets.at(j), ls)) {
				if (BS_same.count(RRBIt.at(i).first) != 0) {
					return true;
				}
                BT.insert(RRBIt.at(i).first);
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

    // cout << "BT: ";
    // for (unordered_set<int>::iterator i = BT.begin(); i != BT.end(); i++) {
    //         cout << *i << ", ";
    // }
    // cout << endl;
    // cout << "BFS Across cluters " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 4 ) << endl;
    int N = graph->getNumberOfVertices();
    // map<VertexID, int> BNToID;
    // int count = 0;
    // for (int i = 0; i < N; i++) {
    //     if (isBoundaryNode[i]) {
    //         BNToID[i] = count;
    //         count ++;
    //     }
    // }
    dynamic_bitset<> visited = dynamic_bitset<>(N);

    int totalVisitedBN = 0;
	/*int totalNumofPasses = 0;*/

    while (!(BS_same.empty() && BS_diff.empty())) {
        unordered_set<VertexID> BS1_same;
		unordered_set<VertexID> BS1_diff;
		// inline the result check at the point of insertion
		//for (unordered_set<int>::iterator i = BS_same.begin(); i != BS_same.end(); ++i) {
		//	if (BT.count(*i) != 0) {
		//	/*	cout << "Total Visted BN: " << totalVisitedBN;
		//		cout << "    Total Number of passes: " << totalNumofPasses << " " << print_digits(getCurrentTimeInMilliSec() - queryStartTime, 4) << endl;*/
		//		return true;
		//	}
		//}
		//for (unordered_set<int>::iterator i = BS_diff.begin(); i != BS_diff.end(); ++i) {
		//	if (BT.count(*i) != 0) {
		//		/*cout << "Total Visted BN: " << totalVisitedBN;
		//		cout << "    Total Number of passes: " << totalNumofPasses << " " << print_digits(getCurrentTimeInMilliSec() - queryStartTime, 4) << endl;*/
		//		return true;
		//	}
		//}

		/*totalVisitedBN += BS_same.size() + BS_diff.size();
		totalNumofPasses += 1;*/

        for (unordered_set<VertexID>::iterator i = BS_same.begin(); i != BS_same.end(); ++i) {
			totalVisitedBN ++;
            cout << "Viewing queued element v in BS_same: " << *i << endl;

            SmallEdgeSet ses = ROBI[*i];
            // graph->getOutNeighbours(*i, ses);
            for (unsigned int j = 0, sizeJ = ses.size(); j != sizeJ; ++j) {
                VertexID v2 = ses[j].first;
                LabelSet ls2 = ses[j].second;

                // Unvisited boundary nodes from a different cluster that can reach target
                if (isLabelSubset(ls2, ls) && (visited[v2] == 0) && (X.find(vToCID[v2]) != X.end())) {
					// inline the result check for efficiency purpose
					if (BT.count(v2) != 0) {
						cout << totalVisitedBN << endl;
						return true;
					}
                    BS1_diff.insert(v2);
                    // cout << "Insert into BS1_diff: " << v2 << endl;
                    visited[v2] = 1;
                }
            }
        }

		for (unordered_set<VertexID>::iterator i = BS_diff.begin(); i != BS_diff.end(); ++i) {
			cout << "Viewing queued element v in BS_diff: " << *i << endl;
			totalVisitedBN++;
			//vector<pair<VertexID, vector<LabelSet>>> RBIi = RBI.at(*i);
			//// cout << "RBIi size: " << RBIi.size() << endl;
			//for (unsigned int j = 0, sizeJ = RBIi.size(); j != sizeJ; ++j) {
			//	unsigned int globalVID = RBIi.at(j).first;
			//	// cout << "reachable boundary node: " << globalVID << endl;
			//	if (visited[globalVID] == 0) {
			//		vector<LabelSet> lss = RBIi.at(j).second;
			//		for (unsigned int k = 0, sizeK = lss.size(); k != sizeK; k++) {
			//			if (isLabelSubset(lss.at(k), ls)) {
			//				// inline the result check for efficiency purpose
			//				if (BT.count(globalVID) != 0) {
			//					cout << totalVisitedBN << endl;
			//					return true;
			//				}
			//				BS1_same.insert(globalVID);
			//				// cout << "Insert into BS1_same: " << globalVID << endl;
			//				visited[globalVID] = 1;
			//				break;
			//			}
			//		}
			//	}
			//}

			unordered_set<VertexID> preliminaryBS = newRBI[*i][ls];
			for (const auto& v : preliminaryBS) {
				cout << "preliminary BS " << v << endl;
				if (visited[v] == 0) {
					if (BT.count(v) != 0) {
						cout << totalVisitedBN << endl;
						return true;
					}
					BS1_same.insert(v);
					visited[v] = 1;
				}
			}

			SmallEdgeSet ses = ROBI[*i];
			// graph->getOutNeighbours(*i, ses);
			for (unsigned int j = 0, sizeJ = ses.size(); j != sizeJ; ++j) {
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

				// Unvisited boundary nodes from a different cluster that can reach target
				if (isLabelSubset(ls2, ls) && (visited[v2] == 0) && (X.find(vToCID[v2]) != X.end())) {
					if (BT.count(v2) != 0) {
						cout << totalVisitedBN << endl;
						return true;
					}
					BS1_diff.insert(v2);
					// cout << "Insert into BS1_diff: " << v2 << endl;
					visited[v2] = 1;
				}
			}
		}

        BS_same = BS1_same;
		BS_diff = BS1_diff;

        // cout << "New BS_same: ";
        // for (unordered_set<int>::iterator i = BS_same.begin(); i != BS_same.end(); i++) {
        //     cout << *i << ", ";
        // }
        // cout << endl;
    }

    /*cout << "Total Visted BN: " << totalVisitedBN;
    cout << "    Total Number of passes: " << totalNumofPasses << " " << print_digits( getCurrentTimeInMilliSec() - queryStartTime, 7 ) << endl;*/

    return false;
};
