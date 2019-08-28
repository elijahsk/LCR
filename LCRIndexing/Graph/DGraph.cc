#include "DGraph.h"

#include <climits>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <queue>
#include <map>
#include <sys/time.h>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <unistd.h>
#include <unordered_set>

using namespace std;

void DGraph::construct(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges) {
    this->constStartTime = getCurrentTimeInMilliSec();

    this->allowMultipleEdges = allowMultipleEdges;
    if ( pN == -1 || pL == -1 ) {
        loadEdgeStats(edgeSet);
    } else {
        this->N = pN;
        this->L = pL;
        this->M = edgeSet->size();
    }

    buildGraph(edgeSet);
    this->constEndTime = getCurrentTimeInMilliSec();
    this->allowMultipleEdges = false;
    //cout << "DGraph |V|=" << N << ",|E|=" << M << ",|L|=" << L << endl;
    //cout << "DGraph size(byte)=" << getGraphSizeInBytes() << ", time(s)=" << getGraphConstructionTime() << endl;
}

DGraph::DGraph(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges) {
    construct(edgeSet, pN, pL, allowMultipleEdges);
};

DGraph::DGraph(EdgeSet* edgeSet, int pN, int pL) {
    construct( edgeSet, pN, pL, false );
};

DGraph::DGraph(EdgeSet* edgeSet) {
    construct( edgeSet, -1, -1, false );
};

DGraph::DGraph(string fileName) {
    construct( loadEdgeFile(fileName), -1, -1, false );
};

DGraph::DGraph(SmallEdgeSets inE, SmallEdgeSets outE, int pN, int pL, int pM, vector<int> weight) {
    this->constStartTime = getCurrentTimeInMilliSec();
    this->N = pN;
    this->L = pL;
    this->M = pM;
    this->inE = inE;
    this->outE = outE;
    this->weight = weight;
    this->allowMultipleEdges = false;
    this->constEndTime = getCurrentTimeInMilliSec();
}

DGraph::~DGraph() {

};

EdgeSet* DGraph::loadEdgeFile(string fileName) {
    cout << "DGraph fileName=" << fileName << endl;

    EdgeSet* edgeSet = new EdgeSet;

    string line;
    VertexID f, t;
    LabelID l2;

    ifstream edge_file (fileName);
    if (edge_file.is_open()) {
        while ( getline (edge_file, line) ) {
            istringstream iss(line);
            string from, label, to;
            iss >> from >> to >> label;
            istringstream (from) >> f;
            istringstream (to) >> t;
            istringstream (label) >> l2;

            l2 = labelIDToLabelSet(l2);

            //cout << "loadEdgeFile: (from,to,label) = " << to_string(f) << ","
            //    << to_string(t) << "," << to_string(l2) << endl;

            Edge edge = make_pair(f , make_pair(t, l2));
            edgeSet->push_back( edge );
        }
        edge_file.close();
    } else {
        cerr << "loadEdgeFile: Unable to open file " << fileName;
    }

    return edgeSet;
};

bool DGraph::findInsertablePosition(graphns::VertexID w, SmallEdgeSet& ses, int& pos) {
    if ( ses.size() == 0 ) {
        return false;
    }

    int low = 0;
    int high = ses.size() - 1;
    int mid = floor( (low + high) / 2 );

    while ( high >= low ) {
        mid = floor( (low + high) / 2 );
        //cout << "findTupleInTuples loop mid=" << mid << " low=" << low << " high=" << high << " tus[mid].first=" << tus[mid].first << endl;

        if ( ses[mid].first == w ) {
            pos = mid;
            return true;
        }

        if ( ses[mid].first > w ) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }

    }

    pos = mid;
    if ( ses[mid].first < w)
        pos++;

    //cout << "findTupleInTuples i=" << i << endl;
    return false;
};

void DGraph::insertEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel, SmallEdgeSet& ses) {
    // inserts an edge (v,w,newLabel) into ses keeping ses sorted
    int pos = 0;
    bool b1 = findInsertablePosition(w, ses, pos);

    if ( b1 == false ) {
        LabelSet ls = newLabel;
        ses.insert( ses.begin() + pos , make_pair(w, ls) );
    }
};

void DGraph::buildGraph(graphns::EdgeSet* edges) {
    // the edgeset does not have to be sorted but for performance this
    // is recommended

    // first initialize outE and inE
    outE.clear();
    inE.clear();

    for (int i = 0; i < N; i++) {
        SmallEdgeSet ses1 = SmallEdgeSet();
        SmallEdgeSet ses2 = SmallEdgeSet();
        outE.push_back(ses1);
        inE.push_back(ses2);
    }

    // add edges to outE and inE
    int i = 0;
    while ( i < edges->size() ) {
        Edge edge = edges->at(i);
        VertexID v = edge.first;
        VertexID w = edge.second.first;
        LabelID lID = edge.second.second;

        //cout << "buildGraph v=" << v << ",w=" << w << ",lID=" << lID << endl;

        insertEdge(v, w, lID, outE[v]);
        insertEdge(w, v, lID, inE[w]);

        i++;
    }

    //cout << "buildGraph completed" << endl;
};

void DGraph::loadEdgeStats(EdgeSet* edgeSet) {
    M = edgeSet->size();

    VertexID last = 0;
    VertexID maxvID = 0;
    set< LabelID > labels;

    int i = 0;
    countPerLabel = vector< long >();

    while ( i < edgeSet->size() ) {
        Edge e = edgeSet->at(i);
        VertexID v = e.first;
        VertexID w = e.second.first;
        LabelID label = e.second.second;

        if ( v > maxvID || w > maxvID ) {
            maxvID = std::max(v, w);
        }

        labels.insert(label);
        while ( countPerLabel.size() < label + 1 ) {
            countPerLabel.push_back( 0 );
        }
        countPerLabel[ label ] += 1;

        i++;
    }

    N = maxvID + 1;
    L = labels.size();

};

int DGraph::getGraphSizeInBytes() {
    int size = 0;
    for (int i = 0; i < N; i++) {
        //cout << "getGraphSizeInBytes outE[i].size()=" << outE[i].size() << ",i=" << i << endl;
        for (int j = 0; j < outE[i].size(); j++) {
            //cout << "getGraphSizeInBytes size=" << size << ",i=" << i << ",j=" << j << endl;
            size += sizeof(outE[i][j].first);
            size += sizeof(outE[i][j].second);
        }
    }

    size *= 2; // inE and outE have the same 'size'

    return size;
};

double DGraph::getGraphConstructionTime() {
    return max(constEndTime - constStartTime, 0.00000001);
};

void DGraph::getOutNeighbours(graphns::VertexID w, SmallEdgeSet& outNeighbours) {
    outNeighbours.clear();
    outNeighbours = outE[w];
};

void DGraph::getInNeighbours(graphns::VertexID w, SmallEdgeSet& inNeighbours) {
    inNeighbours.clear();
    inNeighbours = inE[w];
};

void DGraph::getAllNeighbours(graphns::VertexID w, SmallEdgeSet& allNeighbours) {
    allNeighbours.clear();
    allNeighbours = inE[w];
    allNeighbours.insert(  allNeighbours.end(), outE[w].begin(), outE[w].end() );
}

// prints stats of the graph
std::string DGraph::toString() {
    string output = "";
    output += "|V| = " + to_string(N) + "\n";
    output += "|E| = " + to_string(M) + "\n";
    output += "|L| = " + to_string(L) + "\n";
    output += "size in bytes = " + to_string(getGraphSizeInBytes()) + "\n";
    output += "--------------------------\n";

    //cout << "toString output=" << output << endl;

    for (int i = 0; i < N; i++) {
        SmallEdgeSet outN;
        getOutNeighbours(i, outN);

        output += "out(" + to_string(i) + ")= { \n";

        //cout << "i=" << i << ",outN.size()=" << outN.size() << endl;
        for (int j = 0; j < outN.size(); j++) {
            output += "(" + to_string(outN[j].first) + ",<" + labelSetToLetter(outN[j].second)
                      + "==" + labelSetToString(outN[j].second) + ">)\n";
        }

        output += "} \n";
    }

    for (int i = 0; i < N; i++) {
        SmallEdgeSet inN;
        getInNeighbours(i, inN);

        output += "in(" + to_string(i) + ")= { \n";

        //cout << "i=" << i << ",outN.size()=" << outN.size() << endl;
        for (int j = 0; j < inN.size(); j++) {
            output += "(" + to_string(inN[j].first) + ",<" + labelSetToLetter(inN[j].second)
                      + "==" + labelSetToString(inN[j].second) + ">)\n";
        }

        output += "} \n";
    }

    //cout << "toString2 output=" << output << endl;

    return output;
};

long DGraph::getCountPerLabel(LabelID l) {
    if ( l < 0 || l >= L ) {
        cerr << "getCountPerLabel l out of bounds" << endl;
        return -1;
    }

    return countPerLabel[l];
};

int DGraph::getNumberOfVertices() {
    return N;
};

int DGraph::getNumberOfLabels() {
    return L;
};

int DGraph::getNumberOfEdges() {
    return M;
};

void DGraph::getInE(SmallEdgeSets& tempInE) {
    tempInE = inE;
}

void DGraph::getOutE(SmallEdgeSets& tempOutE) {
    tempOutE = outE;
}

void DGraph::getInE(VertexID v, SmallEdgeSet& tempInE) {
    tempInE = inE[v];
}

void DGraph::getOutE(VertexID v, SmallEdgeSet& tempOutE) {
    tempOutE = outE[v];
}

void DGraph::addNode() {
    N += 1;
    outE.push_back( SmallEdgeSet() );
    inE.push_back( SmallEdgeSet() );
    weight.push_back(-1);
};

void DGraph::removeNode(graphns::VertexID w) {
    if ( w < 0 || w > N ) {
        cerr << " DGraph::removeNode out of bounds w" << w << endl;
        return;
    }

    // remove w and decrease id's of all successors
    for (int i = 0; i < N; i++) {
        if ( i == w )
            continue;

        for (int j = 0; j < outE[i].size(); j++) {
            VertexID w1 = outE[i][j].first;

            if ( w1 > w ) {
                w1--;
            } else if ( w1 == w ) {
                outE[i].erase( outE[i].begin() + j );
                j--;
            }
        }

        for (int j = 0; j < inE[i].size(); j++) {
            VertexID w1 = inE[i][j].first;

            if ( w1 > w ) {
                w1--;
            } else if ( w1 == w ) {
                inE[i].erase( inE[i].begin() + j );
                j--;
            }
        }
    }

    // delete edges of w
    M -= outE[w].size() + inE[w].size();

    outE.erase( outE.begin() + w );
    inE.erase( inE.begin() + w );
    N -= 1;
};

void DGraph::addEdge(VertexID v, VertexID w, LabelID newLabel) {
    if ( v < 0 || v > N || w < 0 || w > N ) {
        cerr << " DGraph::addEdge v or w out of bounds w=" << w  << ",v=" << v << endl;
        return;
    }

    int pos1 = 0;
    int pos2 = 0;

    bool b1 = findInsertablePosition(w, outE[v], pos1);
    bool b2 = findInsertablePosition(v, inE[w], pos2);

    // cout << "b1: " << b1 << endl;
    // cout << "pos1: " << pos1 << endl;
    // cout << "b2: " << b2 << endl;
    // cout << "pos2: " << pos2 << endl;

    // allowMultipleEdges = true;
    if ( (b1 == true || b2 == true) && allowMultipleEdges == false ) {
        // cerr << " DGraph::addEdge vw edge already exists w=" << w  << ",v=" << v << endl;
        return;
    }
    // allowMultipleEdges = false;

    LabelSet ls = labelIDToLabelSet( newLabel );

    outE[v].insert( outE[v].begin() + pos1, make_pair(w, ls));
    inE[w].insert( inE[w].begin() + pos2, make_pair(v, ls));

    M += 1;
};

void DGraph::addMultiEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelSet newLabelSet) {
    if ( v < 0 || v > N || w < 0 || w > N ) {
        //cerr << " DGraph::addMultiEdge v or w out of bounds w=" << w  << ",v=" << v << endl;
        return;
    }

    if ( hasMultiEdge(v, w, newLabelSet ) == true ) {
        //cerr << " DGraph::addMultiEdge newLabelSete already exists newLabelSet=" << newLabelSet << endl;
        return;
    }

    int pos1 = 0;
    int pos2 = 0;

    bool b1 = findInsertablePosition(w, outE[v], pos1);
    bool b2 = findInsertablePosition(v, inE[w], pos2);

    if ( (b1 == true || b2 == true) && allowMultipleEdges == false ) {
        //cerr << " DGraph::addMultiEdge vw edge already exists w=" << w  << ",v=" << v << endl;
        return;
    }

    outE[v].insert( outE[v].begin() + pos1, make_pair(w, newLabelSet));
    inE[w].insert( inE[w].begin() + pos2, make_pair(v, newLabelSet));
    M += 1;
};

void DGraph::removeEdge(graphns::VertexID v, graphns::VertexID w) {
    if ( v < 0 || v > N || w < 0 || w > N ) {
        cout << " DGraph::removeEdge v or w out of bounds w=" << w  << ",v=" << v << endl;
        return;
    }

    int pos1 = 0;
    int pos2 = 0;

    bool b1 = findInsertablePosition(w, outE[v], pos1);
    bool b2 = findInsertablePosition(v, inE[w], pos2);

    if ( b1 == false || b2 == false ) {
        cout << " DGraph::removeEdge vw edge does not exist w=" << w  << ",v=" << v << ",b1=" << b1 << ",b2=" << b2 << endl;
        cout << endl;
        return;
    }

    outE[v].erase( outE[v].begin() + pos1 );
    inE[w].erase( inE[w].begin() + pos2 );
    M -= 1;
    // cout << "removeEdge w=" << w  << ",v=" << v << ",b1=" << b1 << ",b2=" << b2 << ",pos1=" << pos1 << ",pos2=" << pos2 << endl;
};

void DGraph::changeLabel(graphns::VertexID v, graphns::VertexID w, LabelID newLabel) {
    if ( v < 0 || v > N || w < 0 || w > N ) {
        cerr << " DGraph::changeLabel v or w out of bounds w=" << w  << ",v=" << v << endl;
        return;
    }

    int pos1 = 0;
    int pos2 = 0;

    bool b1 = findInsertablePosition(w, outE[v], pos1);
    bool b2 = findInsertablePosition(v, inE[w], pos2);

    if ( b1 == false || b2 == false ) {
        cerr << " DGraph::changeLabel vw edge does not exist" << endl;
        return;
    }

    LabelSet ls = labelIDToLabelSet( newLabel );
    outE[v][pos1].second = ls;
    inE[w][pos2].second = ls;
};

LabelID DGraph::getLabelID(graphns::VertexID v , graphns::VertexID w) {
    if ( hasEdge(v, w) == false ) {
        cerr << " DGraph::getLabelID v or w out of bounds w=" << w  << ",v=" << v << endl;
        return 0;
    }

    int pos1 = 0;
    bool b1 = findInsertablePosition(w, outE[v], pos1);

    return labelSetToLabelID( outE[v][pos1].second );
}

bool DGraph::hasMultiEdge(graphns::VertexID v , graphns::VertexID w, graphns::LabelSet ls) {
    // verifies whether there already exists an edge (v,w,ls') s.t. ls' subset of ls
    if ( hasEdge(v, w) == false ) {
        return false;
    }

    int pos1 = 0;
    bool b1 = findInsertablePosition(w, outE[v], pos1);

    if ( pos1 > 0 ) {
        while (outE[v][pos1 - 1].first == w) {
            pos1 -= 1;
            if ( pos1 <= 0 ) {
                break;
            }
        }
    }

    while (outE[v][pos1].first == w) {
        LabelSet lsPrime = outE[v][pos1].second;

        if ( isLabelSubset(lsPrime, ls) == true ) {
            return true;
        }

        pos1 += 1;
        if ( pos1 >= outE[v].size() ) {
            break;
        }
    }

    return false;
}

void DGraph::setWeight(vector<int> w) {
    weight = w;
}

void DGraph::setWeight(VertexID v, int w) {
    weight[(int) v] = w;
}

void DGraph::getWeight(VertexID v, int& tempWeight) {
    if (v < N) {
        tempWeight = weight[(int) v];
    } else {
        tempWeight = -2;
    }

}

void DGraph::getWeights(vector<int>& weights) {
    weights = weight;
}

double DGraph::computeAverageDiameter() {
    // a simple and dumb way to compute the diameter
    int totalDiameter = 0;
    for (int i = 0; i < N; i++) {
        //cout << "computeDiameter: i=" << i << endl;

        dynamic_bitset<> marked = dynamic_bitset<>(N);
        queue< pair<VertexID, int> > q;
        q.push( make_pair(i, 0) );

        while ( q.empty() == false ) {
            auto p = q.front();
            q.pop();

            if ( marked[p.first] == 1 ) {
                continue;
            }
            marked[p.first] = 1;
            totalDiameter++;

            SmallEdgeSet ses;
            getOutNeighbours(p.first, ses);
            for (int j = 0; j < ses.size(); j++) {
                q.push( make_pair(ses[j].first, p.second + 1) );
            }

        }
    }

    return ( (double) totalDiameter / (double) N );
}

int DGraph::computerNumberOfTriangles() {
    int NoOfTriangles = 0;
    for (int i = 0; i < N; i++) {
        //cout << "computerNumberOfTriangles: i=" << i << endl;

        dynamic_bitset<> marked = dynamic_bitset<>(N);
        queue< pair<VertexID, int> > q;
        q.push( make_pair(i, 0) );

        while ( q.empty() == false ) {
            auto p = q.front();
            q.pop();

            if ( p.first == i && p.second == 3 ) {
                //cout << "NoOfTriangles=" << NoOfTriangles << endl;
                NoOfTriangles++;
            }

            if ( marked[p.first] == 1 || p.second >= 3 ) {
                continue;
            }
            marked[p.first] = 1;

            SmallEdgeSet ses;
            getOutNeighbours(p.first, ses);
            for (int j = 0; j < ses.size(); j++) {
                q.push( make_pair(ses[j].first, p.second + 1) );
            }
        }
    }
    return NoOfTriangles / 3;
}

bool DGraph::hasEdge(VertexID v, VertexID w) {
    SmallEdgeSet ses;
    getOutNeighbours(v, ses);

    for (int i = 0; i < ses.size(); i++) {
        if ( ses[i].first == w ) {
            return true;
        }
    }

    return false;
};

double DGraph::computeClusterCoefficient() {
    double clusterCoefficient = 0.0;
    for (int i = 0; i < N; i++) {
        // compute neighbourhood of i
        SmallEdgeSet outSes, inSes;
        getOutNeighbours(i, outSes);
        getInNeighbours(i, inSes);

        int Ki = outSes.size() + inSes.size();
        double localClusterCoefficient = 0.0;

        for (int j = 0; j < inSes.size(); j++) {
            VertexID v = inSes[j].first;
            for (int k = 0; k < outSes.size(); k++) {
                // we need a (w,v) edge for a triangle
                VertexID w = outSes[k].first;

                if ( hasEdge(w, v) == true ) {
                    localClusterCoefficient += 1.0;
                }
            }
        }

        localClusterCoefficient /= (Ki * (Ki - 1) * 1.0);
        if ( isnormal(localClusterCoefficient) == 0 ) {
            clusterCoefficient += localClusterCoefficient;
        }
    }

    return max(0.0, clusterCoefficient / N);
}

void DGraph::getNodesWithNoInDegree(DGraph* tempGraph, vector<VertexID>& nodesWithNoInDegree) {
    nodesWithNoInDegree.clear();
    // cout << N << endl;
    int size = N;
    SmallEdgeSets tempInE;
    // cout << "test1" << endl;
    tempGraph->getInE(tempInE);
    // cout << "test2" << endl;
    int weight;
    for (int i = 0; i < size; i++) {
        // tempInE[i].size() == 0 &&
        // cout << i << " " << endl;
        tempGraph->getWeight(i, weight);
        // cout << "test3" << endl;
        if (weight > 0) {
            // cout << "test4" << endl;
            nodesWithNoInDegree.push_back((VertexID) i);
        }
        // cout << "test5" << endl;
    }

    // for (int i = 0, sizeI = nodesWithNoInDegree.size(); i != sizeI; ++i) {
    //     cout << nodesWithNoInDegree[i] << " ";
    // }
    // cout << endl;
    // cout << "test6" << endl;
}


void DGraph::findLongestChain(DGraph* tempGraph, VertexID v, vector<VertexID>& chain) {
    // cout << "v: " << v << endl;
    vector<int> queue;
    vector<int> prev;
    vector<bool> visited;
    for (int i = 0, sizeI = tempGraph->getNumberOfVertices(); i != sizeI; ++i) {
        visited.push_back(false);
    }
    // cout << "test1: " << tempGraph->getNumberOfVertices() << endl;
    // ================================================
    SmallEdgeSets tempOutE;
    tempGraph->getOutE(tempOutE);
    SmallEdgeSets tempInE;
    tempGraph->getInE(tempInE);
    // Does chaning tempOutE affect the original one?
    // Ans: No
    // ================================================
    // cout << "test2" << endl;

    // use nodeID + 1 to avoid confusion caused by 0
    int N = tempGraph->getNumberOfVertices();
    queue.push_back(v + 1);
    prev.push_back(N);
    visited[v] = true;

    // Add head info in both directions
    SmallEdgeSet outEdges = tempOutE[v];
    SmallEdgeSet inEdges = tempInE[v];
    VertexID nextV;



    for (int i = 0, sizeI = outEdges.size(); i != sizeI; ++i) {
        nextV = outEdges[i].first;
        if (!visited[nextV]) {
            queue.push_back(nextV + 1);
            prev.push_back(0);
            visited[nextV] = true;
        }
    }

    for (int i = 0, sizeI = inEdges.size(); i != sizeI; ++i) {
        nextV = inEdges[i].first;
        if (!visited[nextV]) {
            queue.push_back(-nextV - 1);
            prev.push_back(0);
            visited[nextV] = true;
        }
    }


    int head = 1;
    while (head < queue.size()) {
        // cout << "head: " << head << endl;
        int currV = queue[head];
        // cout << "currV: " << currV << endl;
        if (currV > 0) {
            currV --;
            // cout << "currV > 0: " << currV << endl;
            outEdges = tempOutE[currV];
            // cout << "outEdges.size(): " << outEdges.size() << endl;
            for (int i = 0, j = outEdges.size(); i != j; ++i) {
                // cout << "i: " << i << endl;
                nextV = outEdges[i].first;
                if (!visited[nextV]) {
                    // cout << "next: " << nextV << endl;
                    // cout << "queue size: " << queue.size() << endl;

                    queue.push_back(nextV + 1);

                    // cout << "test3" << endl;
                    prev.push_back(head);
                    visited[nextV] = true;
                    // cout << "test4" << endl;
                }
                // cout << "done " << i << endl;
            }
        } else {
            currV = -currV - 1;
            // cout << "currV < 0: " << currV << endl;
            visited[currV] = true;
            inEdges = tempInE[currV];

            for (int i = 0, sizeI = inEdges.size(); i != sizeI; ++i) {
                nextV = inEdges[i].first;
                if (!visited[nextV]) {
                    queue.push_back(-nextV - 1);
                    prev.push_back(head);
                    visited[nextV] = true;
                }
            }
        }

        head += 1;
        // cout << "next head" << head << endl;
    }

    int head1 = head - 1;
    int head2 = head - 1;
    // find the last node that is in the path of in degree
    while (queue[head1] > 0) {
        head1 --;
    }
    // find the last node that is in the path of out degree
    while (queue[head2] < 0) {
        head2 --;
    }
    // cout << queue.size() << " " << head << endl;

    while (head2 >= 0 && prev[head2] != N) {
        chain.push_back(queue[head2] - 1);
        head2 = prev[head2];
        // cout << head << " ";
    }

    // cout << endl;

    chain.push_back(v);

    vector<VertexID> tempChain;
    while (head1 >= 0 && prev[head1] != N) {
        // cout << queue[head1] << endl;
        tempChain.push_back(-queue[head1] - 1);
        head1 = prev[head1];
    }

    reverse(tempChain.begin(), tempChain.end());

    chain.insert(chain.end(), tempChain.begin(), tempChain.end());

}

void DGraph::findLongestChain(DGraph* tempGraph, vector<VertexID> nodesWithNoInDegree, vector<VertexID>& chain) {
    chain.clear();
    // randomly choose X nodes to find the longest chain among them
    int populationSize = nodesWithNoInDegree.size();
    cout << "Population Size: " << populationSize << endl;
    // populationSize = N;
    int num = min((int) sqrt(populationSize), 30);
    int randomV = rand();
    cout << "Number of randomV: " << num << endl;
    for (int i = 0; i < num; i++) {
        // May give restrictions to the randomV, e.g. in & out degree
        srand(time(0) * randomV);
        randomV = rand() % populationSize;
        // cout << "randomV: " << nodesWithNoInDegree[randomV] << endl;
        // cout << "randomV: " << randomV << endl;
        vector<VertexID> tempChain;
        findLongestChain(tempGraph, nodesWithNoInDegree[randomV], tempChain);
        // vector<VertexID> tempChain = findLongestChain(tempGraph, randomV);
        // cout << "Chain length: " << tempChain.size() << endl;
        // for (int i = 0, sizeI = tempChain.size(); i != sizeI; ++i) {
        //     cout << tempChain[i] << " ";
        // }
        // cout << endl;
        if (tempChain.size() > chain.size()) {
            chain = tempChain;
        }
    }
    cout << "Chain length: " << chain.size() << endl;
}

// naive approach, length = m * r + remainder
// more complex approach length = (m + 1) segments of similar length
void DGraph::addSegments(int head, int tail, vector<VertexID> chain, vector<vector<VertexID>>& segments, int radius) {
    // cout << "Add segments" << endl;
    vector<VertexID> segment;
    int segmentsLength = tail - head;
    if (segmentsLength % radius == 0) {
        int count = 0;
        while (head < tail) {
            segment.push_back(chain[head]);
            count += 1;

            if (count == radius) {
                segments.push_back(segment);
                segment.clear();
                count = 0;
            }

            head += 1;
        }
    } else {
        int segmentCount = (int) segmentsLength / radius + 1;
        // cout << "segmentCount: " << segmentCount << endl;
        int approxSegmentLength = (int) segmentsLength / segmentCount;
        int remainder = segmentsLength % approxSegmentLength;
        int count = 0;
        // cout << "Approximate length: " << approxSegmentLength << endl;
        // cout << "Remainder: " << remainder << endl;
        while (head < tail) {
            segment.push_back(chain[head]);
            count += 1;

            if (remainder > 0) {
                if (count == approxSegmentLength + 1) {
                    segments.push_back(segment);
                    // cout << "Segment Size: " << segment.size() << endl;
                    segment.clear();
                    count = 0;
                    remainder -= 1;
                }
            } else if (count == approxSegmentLength) {
                segments.push_back(segment);
                // cout << "Segment Size: " << segment.size() << endl;
                segment.clear();
                count = 0;
            }

            head += 1;
        }
    }

    // cout << "current segments: " << endl;
    // for (int i = 0, sizeI = segments.size(); i != sizeI; ++i) {
    //     for (int j = 0, sizeJ = segments[i].size(); j != sizeJ; ++j) {
    //         cout << segments[i][j] << " ";
    //     }
    //     cout << endl;
    // }
}

// find out partitions in chain that do not contain any supernodes
// and further segment them in addSegments
void DGraph::segmentChain(vector<VertexID> chain, int radius, vector<int> weights, vector<vector<VertexID>>& segments) {
    int head = 0;
    int tail = 0;
    // cout << "segment chain" << endl;
    while (head < chain.size()) {
        // cout << "head & tail: " << head << " " << tail << endl;
        // cout << "tail id: " << chain[tail] << endl;
        // cout << "tail weight: " << weights[chain[tail]] << endl;
        // Find the first node that is a cluster, starting from head
        while (tail != chain.size() && weights[chain[tail]] == 1) {
            tail ++;
        }
        if (head < tail) {
            // cout << "head & tail: " << head << " " << tail << endl;
            addSegments(head, tail, chain, segments, radius);
        }
        head = tail;
        // Find the first node that is not a cluster, starting from tail (since )
        while (head != chain.size() && weights[chain[head]] != 1) {
            head ++;
        }
        tail = head;
    }
}

void DGraph::growSegment(DGraph* tempGraph, VertexID cID, vector<VertexID>& startVertices, int maxClusterSize, vector<vector<VertexID>>& clusters, vector<int>& vToCID) {

    int size = startVertices.size();
    int i = 0;
    int count = 0;
    SmallEdgeSets tempOutE;
    tempGraph->getOutE(tempOutE);
    SmallEdgeSets tempInE;
    tempGraph->getInE(tempInE);
    vector<VertexID> nextVertices;
    while (i < size && count <= maxClusterSize) {
        VertexID v = startVertices[i];
        SmallEdgeSet outEdges = tempOutE[(int) v];
        int edgeCount = outEdges.size();
        for (int j = 0; j != edgeCount; ++j) {
            VertexID currNode = outEdges[j].first;
            // non-supernode that has not been clustered
            // instead of tempGraph->getWeight(currNode) == 1, check if the node is newly added (&& currNode < N)
            // alternatively make cID of all supernodes another specific value?
            // Yep
            if (vToCID[currNode] == -1 ) {
                nextVertices.push_back(currNode);
                clusters[cID].push_back(currNode);
                vToCID[currNode] = cID;
                count += 1;
            }
        }

        SmallEdgeSet inEdges = tempInE[(int) v];
        edgeCount = inEdges.size();
        for (int j = 0; j != edgeCount; ++j) {

            if (count >= maxClusterSize) break;

            VertexID currNode = inEdges[j].first;
            // non-supernode that has not been clustered
            // instead of tempGraph->getWeight(currNode) == 1, check if the node is newly added (&& currNode < N)
            // alternatively make cID of all supernodes another specific value?
            // Yep
            if (vToCID[currNode] == -1 ) {
                nextVertices.push_back(currNode);
                clusters[cID].push_back(currNode);
                vToCID[currNode] = cID;
                count += 1;
            }
        }

        // Break only after all the neighbours are added, so that the current node will not become a boundary node
        if (count >= maxClusterSize) break;

        i += 1;
    }
    startVertices = nextVertices;
}

void DGraph::modifyGraph(DGraph* tempGraph, VertexID v, VertexID w, vector<vector<VertexID>>& clusters, vector<int>& vToCID) {
    // if both are nodes
    if (v < N && w < N) {
        tempGraph->addNode();
        int newN = tempGraph->getNumberOfVertices() - 1;

        vector<VertexID> temp;
        temp.push_back(v);
        temp.push_back(w);
        clusters.push_back(temp);

        int clusterID = clusters.size() - 1;

        vToCID.push_back(clusterID);
        vToCID[v] = clusterID;
        vToCID[w] = clusterID;
        tempGraph->setWeight(newN, 2);

        vector<VertexID> cluster = clusters[clusterID];
        SmallEdgeSets tempOutE;
        SmallEdgeSets tempInE;
        tempGraph->getOutE(tempOutE);
        tempGraph->getInE(tempInE);
        for (int j = 0; j < 2; j++) {
            // Reconnect from supernode to existing out-nodes
            // cout << "Add out edges:" << endl;
            SmallEdgeSet currOutE = tempOutE[cluster[j]];
            for (int k = 0; k < currOutE.size(); k++) {
                if (vToCID[currOutE[k].first] != clusterID) {
                    // cout << "Add edge from, to: " << newN << " " << currOutE[k].first << endl;
                    tempGraph->addEdge(newN, currOutE[k].first, labelSetToLabelID(currOutE[k].second));
                }
                // cout << "Remove edge from, to: " << cluster[j] << " " << currOutE[k].first;
                tempGraph->removeEdge(cluster[j], currOutE[k].first);
            }

            // Reconnect from existing in-nodes to supernode
            // cout << "Add in edges:" << endl;
            SmallEdgeSet currInE = tempInE[cluster[j]];
            for (int k = 0; k < currInE.size(); k++) {
                if (vToCID[currInE[k].first] != clusterID) {
                    // cout << "Add edge from, to: " << currInE[k].first << " " << newN << endl;
                    tempGraph->addEdge(currInE[k].first, newN, labelSetToLabelID(currInE[k].second));
                }
                // cout << "Remove edge from, to: " << currInE[k].first << " " << cluster[j];
                tempGraph->removeEdge(currInE[k].first, cluster[j]);
            }

            // Assign 0 weight to already clustered node
            tempGraph->setWeight(cluster[j], 0);
        }
    } else {
        // To deal with the situation where w is a cluster while v may / may not be a cluster
        // To ensure that we are always merging something into a cluster
        int weight;
        tempGraph->getWeight(w, weight);
        if (weight > 1) {
            VertexID temp = v;
            v = w;
            w = temp;
        }

        // Merge w to v, make w dummy node that directs to v
        // cout << "Merge " << w << " into " << v << endl;
        vToCID[w] = vToCID[v];
        int weightV;
        int weightW;
        tempGraph->getWeight(v, weightV);
        tempGraph->getWeight(w, weightW);
        tempGraph->setWeight(v, weightV + weightW);
        tempGraph->setWeight(w, 0);

        SmallEdgeSet currOutE;
        tempGraph->getOutE(w, currOutE);
        // cout << "Add out edges:" << endl;
        for (int k = 0, sizeK = currOutE.size(); k != sizeK; ++k) {
            // cout << k << endl;
            if (currOutE[k].first != v) {
                // cout << "Add edge from, to: " << v << " " << currOutE[k].first << endl;
                tempGraph->addEdge(v, currOutE[k].first, labelSetToLabelID(currOutE[k].second));
            }
            // cout << "Remove edge from, to: " << w << " " << currOutE[k].first;
            tempGraph->removeEdge(w, currOutE[k].first);
        }

        SmallEdgeSet currInE;
        tempGraph->getInE(w, currInE);
        // cout << "Add in edges:" << endl;
        for (int k = 0, sizeK = currInE.size(); k != sizeK; ++k) {
            // cout << k << endl;
            if (currInE[k].first != v) {
                // cout << "Add edge from, to: " << currInE[k].first << " " << v << endl;
                tempGraph->addEdge(currInE[k].first, v, labelSetToLabelID(currInE[k].second));
            }
            // cout << "Remove edge from, to: " << currInE[k].first << " " << w;
            tempGraph->removeEdge(currInE[k].first, w);
        }
    }
}


void DGraph::modifyGraph(DGraph* tempGraph, int segmentCount, vector<vector<VertexID>>& clusters, vector<int>& vToCID) {
    int size = clusters.size();
    int oldSize = size - segmentCount;

    for (int i = 0; i < segmentCount; i++) {
        int cID = oldSize + i;
        // cout << "Modify cluster " << cID << ": " << endl;

        // Add new node, set its weight to be cluster size
        tempGraph->addNode();
        VertexID newNodeID = (VertexID) tempGraph->getNumberOfVertices() - 1;
        tempGraph->setWeight(newNodeID, clusters[cID].size());
        // cout << "cluster size: " << clusters[cID].size() << endl;
        // cout << "new node id: " << newNodeID << endl;
        // cout << "weight of new node: " << tempGraph->getWeight(newNodeID) << endl;
        // Set CID to be the cluster ID of member nodes
        vToCID.push_back(cID);

        vector<VertexID> cluster = clusters[cID];
        int numNodes = cluster.size();
        for (int j = 0; j < numNodes; j++) {
            // cout << cluster[j] << " ";
            // cout << "node " << j << ": " << cluster[j] << endl;
            // Reconnect from supernode to existing out-nodes
            SmallEdgeSet currOutE;
            tempGraph->getOutE(cluster[j], currOutE);

            // cout << "Add out edges: " << endl;
            for (int k = 0, sizeK = currOutE.size(); k != sizeK; k++) {
                // not in the same cluster
                // cout << "k: " << k << endl;
                // cout << currOutE[k].first << endl;
                if (vToCID[currOutE[k].first] != cID) {
                    // cout << "Add edge: from, to, label  " <<  newNodeID << ", " << currOutE[k].first << ", " << currOutE[k].second << endl;
                    tempGraph->addEdge(newNodeID, currOutE[k].first, labelSetToLabelID(currOutE[k].second));
                    // SmallEdgeSet tempInE;
                    // tempGraph->getInE(currOutE[k].first, tempInE);
                    // for (int m = 0, sizeM = tempInE.size(); m != sizeM; ++m) {
                    //     cout << tempInE[m].first << " " << tempInE[m].second << endl;
                    // }

                }
                tempGraph->removeEdge(cluster[j], currOutE[k].first);

            }

            // Reconnect from existing in-nodes to supernode
            // cout << "Add in edges: " << endl;
            SmallEdgeSet currInE;
            tempGraph->getInE(cluster[j], currInE);
            for (int k = 0, sizeK = currInE.size(); k != sizeK; k++) {
                // not in the same cluster
                if (vToCID[currInE[k].first] != cID) {
                    // cout << "Add edge: from, to, label  " <<  currInE[k].first << ", " << newNodeID << ", " << currInE[k].second << endl;
                    tempGraph->addEdge(currInE[k].first, newNodeID, labelSetToLabelID(currInE[k].second));
                    // cout << "checking 797.." << endl;
                    // curroutE = tempGraph->getOutE(797);
                }
                tempGraph->removeEdge(currInE[k].first, cluster[j]);
                // cout << "checking 798.." << endl;
                // curroutE = tempGraph->getOutE(798);
            }

            // Assign 0 weight to already clustered node
            tempGraph->setWeight(cluster[j], 0);
        }
        // cout << "Done with cluster " << cID << endl;

        // tempGraph->getStatus();
        // break;
    }

    // tempGraph->getStatus();
}

void DGraph::getStatus() {
    for (int i = 0, sizeI = N; i != N; i++) {
        cout << "vertex: " << i << " weight: " << weight[i] << endl;
        cout << "outE size: " << outE[i].size() << endl;
        for (int j = 0, sizeJ = outE[i].size(); j != sizeJ; ++j) {
            cout << outE[i][j].first << " " << outE[i][j].second << endl;
        }
        cout << " inE size: " << inE[i].size() << endl;
        for (int j = 0, sizeJ = inE[i].size(); j != sizeJ; j++) {
            cout << inE[i][j].first << " " << inE[i][j].second << endl;
        }
    }
}

void DGraph::newClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID, int radius, int maxClusterSize, int minClusterSize) {
    // *
    // Where should initialization be performed?
    // Create another graph for all the merging and substitution
    // Keep the original copy
    // *

    DGraph* tempGraph = new DGraph(inE, outE, N, L, M, vector<int>(N, 1));
    vector<VertexID> nodesWithNoInDegree;
    vector<VertexID> chain;


    // What is the termination condition?
    // chain length? segmentCount?
    while (true) {
        // cout << "Start" << endl;

        getNodesWithNoInDegree(tempGraph, nodesWithNoInDegree);

        if (N - nodesWithNoInDegree.size() < clusters.size() * minClusterSize) {
            break;
        }

        // cout << "done get nodes" << endl;

        // for (int i = 0, sizeI = nodesWithNoInDegree.size(); i != sizeI; ++i) {
        //     cout << nodesWithNoInDegree[i] << " ";
        // }
        // cout << endl;

        findLongestChain(tempGraph, nodesWithNoInDegree, chain);

        for (int i = 0, sizeI = chain.size(); i != sizeI; ++i) {
            cout << chain[i] << " ";
        }
        cout << endl;

        vector<vector<VertexID>> segments;
        vector<int> weights;
        tempGraph->getWeights(weights);
        segmentChain(chain, radius, weights, segments);

        // cout << "Final segments: " << endl;
        // for (int i = 0, sizeI = segments.size(); i != sizeI; ++i) {
        //     for (int j = 0, sizeJ = segments[i].size(); j != sizeJ; ++j) {
        //         cout << segments[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        int segmentCount = segments.size();
        int oldClusterSize = clusters.size();

        // termination condition
        if (segmentCount == 0) {
            break;
        }

        // // include all elements in each segments into respective clusters
        // for (int i = 0; i < segments.size(); i++) {
        //     vector<VertexID> segment = segments[i];
        //     clusters.push_back(segment);
        //     int cID = clusters.size();
        //     for (int j = 0; j < segment.size(); j++) {
        //         vToCID[(int)segment[j]] = cID;
        //     }
        // }

        // include only the centre nodes into cluster
        // Adjust max cluster size based on radius and segment length
        // to make sure that there is always a place for other nodes in the segment
        vector<vector<VertexID>> startVertices;
        for (int i = 0; i != segmentCount; i++) {
            vector<VertexID> segment = segments[i];
            // VertexID middleID = segment[segment.size() / 2];
            // vector<VertexID> temp{middleID};
            // startVertices.push_back(temp);
            // clusters.push_back(temp);
            startVertices.push_back(segment);
            clusters.push_back(segment);
            for (int j = 0, sizeJ = segment.size(); j != sizeJ; ++j) {
                vToCID[segment[j]] = oldClusterSize + i;
            }
        }

        // cout << "segment count & cluster size: " << segmentCount << " " << clusters.size() << endl;
        // cout << "Middle nodes: ";
        // for (int i = oldClusterSize, sizeI = clusters.size(); i != sizeI; ++i) {
        //     cout << clusters[i][0] << " ";
        // }
        // cout << endl;


        for (int i = 0; i < radius; i++) {
            for (int j = 0; j < segmentCount; j++) {
                // cout << "radius: " << i << endl;
                vector<VertexID> segment = segments[j];
                // VertexID cID = vToCID[segment[segment.size() / 2]];
                VertexID cID = vToCID[segment[0]];
                // total - existing cluster - nodes on segment that have not been included
                // int maxCurrentClusterSize = maxClusterSize - (int) clusters[(int) cID].size() - max((int) segment.size() - (i * 2 + 1), 0);
                int maxCurrentClusterSize = maxClusterSize - clusters[(int) cID].size();

                // cout << "Members of cluster " << cID << ": ";
                // for (int k = 0, sizeK = clusters[(int) cID].size(); k != sizeK; ++k) {
                //     cout << clusters[(int) cID][k] << " ";
                // }
                // cout << endl;

                // cout << "segment " << j << " current max cluster size: " << maxCurrentClusterSize << endl;
                // grow segment when there are vacancy in the cluster and there are nodes to develop
                if (maxClusterSize > 0 && startVertices[j].size() > 0) {
                    growSegment(tempGraph, cID, startVertices[j], maxCurrentClusterSize, clusters, vToCID);
                }

                // if (i + 1 + segment.size() / 2 < segment.size()) {
                //     VertexID v1 = segment[segment.size() / 2 + i + 1];
                //     if (vToCID[v1] == -1) {
                //         vToCID[v1] == cID;
                //         clusters[cID].push_back(v1);
                //         startVertices[j].push_back(v1);
                //     }
                // }

                // if (segment.size() / 2 >= i + 1) {
                //     VertexID v2 = segment[segment.size() / 2 - i - 1];
                //     if (vToCID[v2] == -1) {
                //         vToCID[v2] = cID;
                //         clusters[cID].push_back(v2);
                //         startVertices[j].push_back(v2);
                //     }
                // }
            }
        }

        // Add rest of the segments to the cluster, if they have not been added in growSegment
        // for (int i = 0; i < segmentCount; i++) {
        //     vector<VertexID> segment = segments[i];
        //     if (clusters[oldClusterSize + i].size() < maxClusterSize) {
        //         // for segments with even no. of nodes 0, 1, 2, ... 2n - 1,
        //         // middle node is n, n nodes in front, n - 1 nodes at the back
        //         // check segment[0] and the rest should be in pairs
        //         // for segments with odd no. of nodes, the nodes should already be paired
        //         int head = 0;
        //         int tail = segment.size() - 1;

        //         if (segments.size() % 2 == 0) {
        //             if (vToCID[segment[0]] == -1) {
        //                 vToCID[segment[0]] = oldClusterSize + i;
        //                 clusters[oldClusterSize + i].push_back(segment[0]);
        //                 cout << "Add " << segment[0] << " to cluster " << oldClusterSize + i << endl;
        //                 head = 1;
        //             } else {
        //                 continue;
        //             }
        //         }

        //         // if any of the pair has been clustered, all the more centered nodes have been clustered
        //         bool notClustered = true;
        //         while (notClustered) {
        //             if (vToCID[segment[head]] == -1) {
        //                 vToCID[segment[head]] = oldClusterSize + i;
        //                 clusters[oldClusterSize + i].push_back(segment[head]);
        //                 cout << "Add " << segment[head] << " to cluster " << oldClusterSize + i << endl;
        //             } else {
        //                 notClustered = false;
        //             }
        //             if (vToCID[segment[tail]] == -1) {
        //                 vToCID[segment[tail]] = oldClusterSize + i;
        //                 clusters[oldClusterSize + i].push_back(segment[tail]);
        //                 cout << "Add " << segment[tail] << " to cluster " << oldClusterSize + i << endl;
        //             } else {
        //                 notClustered = false;
        //             }
        //         }
        //     }
        // }


        // vector<VertexID> cluster = clusters[clusters.size() - 1];
        // cout << "Cluster " << clusters.size() - 1 << ": ";
        // for (int j = 0, sizeJ = cluster.size(); j != sizeJ; ++j) {
        //     cout << cluster[j] << " ";
        // }
        // cout << endl;


        cout << "Modify graph" << endl;

        // the new node shouldn't be involved in growing or segmenting, but only finding the longest chain
        modifyGraph(tempGraph, segmentCount, clusters, vToCID);

        cout << "Done modifying graph" << endl;

        // for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
        //     vector<VertexID> cluster = clusters[i];
        //     cout << "Cluster " << i << ": ";
        //     for (int j = 0, sizeJ = cluster.size(); j != sizeJ; ++j) {
        //         cout << cluster[j] << " ";
        //     }
        //     cout << endl;
        // }

        // cout << "vToCID: " << endl;
        // for (int i = 0, sizeI = tempGraph->getNumberOfVertices(); i != sizeI; ++i) {
        //     cout << i << " " << vToCID[i] << endl;
        // }
    }

    cout << "clean up" << endl;

    // merge small clusters / nodes
    // nodes that have not been included into any cluster / cluster that does not reach the minimum cluster size
    // 1. connect to neighbouring unclustered node
    //  - assign a new cluster id
    //  - modify edges to form a supernode
    // 2. if there are no unclustered node, connect to the smallest nearby cluster


    // Bring unclustered nodes to clusters
    int weight;
    SmallEdgeSet tempOutE;
    SmallEdgeSet tempInE;

    for (int i = 0; i < N; i++) {
        if (vToCID[i] == -1) {
            // cout << "Unclustered node: " << i << endl;
            tempGraph->getOutE(i, tempOutE);
            VertexID minVertexID = -1;
            int minWeight = N + 1;
            for (size_t j = 0, size = tempOutE.size(); j != size; ++j) {
                VertexID v = tempOutE[j].first;
                // all nodes that can be connected should have weight > 0
                tempGraph->getWeight(v, weight);
                if (weight < minWeight) {
                    minWeight = weight;
                    minVertexID = v;
                }
            }

            tempGraph->getInE(i, tempInE);
            for (size_t j = 0, sizeJ = tempInE.size(); j != sizeJ ; ++j) {
                VertexID v = tempInE[j].first;
                // all nodes that can be connected should have weight > 0
                tempGraph->getWeight(v, weight);
                if (weight < minWeight) {
                    minWeight = weight;
                    minVertexID = v;
                }
            }

            if (minVertexID != -1) {
                // cout << " merge with " << minVertexID << endl;
                modifyGraph(tempGraph, minVertexID, (VertexID) i, clusters, vToCID);
            }
            // cout << endl;
        }
    }

    int newN = tempGraph->getNumberOfVertices();
    vector<VertexID> smallClusters;
    for (int i = N; i < newN; i++) {
        tempGraph->getWeight(i, weight);
        tempGraph->getOutE(i, tempOutE);
        tempGraph->getInE(i, tempInE);
        if (weight != 0 && weight < minClusterSize && !(tempOutE.empty() && tempInE.empty())) {
            smallClusters.push_back(i);
        }
    }

    // cout << "Small Clusters: ";
    // for (int i = 0, sizeI = smallClusters.size(); i != sizeI; ++i) {
    //     cout << smallClusters[i] << " ";
    // }
    // cout << endl;

    // cout << "Total & Small: " << clusters.size() << " " << smallClusters.size() << endl;

    if (!smallClusters.empty()) {
        for (size_t i = 0, sizeI = smallClusters.size(); i != sizeI; ++i) {
            // the node will only be merged once regardless
            // 1. node become empty if merging is required
            // 2. the node can't be merged as it has no neighbours
            // 3. the node, after some previous merging, has become big enough
            // 4. the node, after all the merging, still does not reach the minimum size,
            //    will be merged with a node that exceeds minimum size.
            tempGraph->getWeight(smallClusters[i], weight);
            if (weight < minClusterSize) {
                // cout << "small cluster: " << smallClusters[i];
                tempGraph->getOutE(smallClusters[i], tempOutE);
                VertexID minVertexID = -1;
                int minWeight = N + 1;
                for (size_t j = 0, sizeJ = tempOutE.size(); j != sizeJ ; ++j) {
                    VertexID v = tempOutE[j].first;
                    // all nodes that can be connected should have weight > 0
                    tempGraph->getWeight(v, weight);
                    if (weight < minWeight) {
                        minWeight = weight;
                        minVertexID = v;
                    }
                }

                tempGraph->getInE(smallClusters[i], tempInE);
                for (size_t j = 0, sizeJ = tempInE.size(); j != sizeJ ; ++j) {
                    VertexID v = tempInE[j].first;
                    // all nodes that can be connected should have weight > 0
                    tempGraph->getWeight(v, weight);
                    if (weight < minWeight) {
                        minWeight = weight;
                        minVertexID = v;
                    }
                }

                if (minVertexID != -1) {
                    // cout << " merge with " << minVertexID << " with weight " << minWeight << endl;
                    modifyGraph(tempGraph, minVertexID, smallClusters[i], clusters, vToCID);
                }
                // cout << endl;
            }
        }

    }

    // For original node 0 ... N - 1
    // 1. vToCID contains the original cluster it is in
    // 2. vToCID == -1, indicating that the node is not connected to any other nodes

    // For cluster node N ...
    // 1. vToCID contains the true cluster id
    // 2. vToCID contains the cluster id of its ancester

    // Given a cluster id I, the position of the cluster node is N + I

    // vToCID[i] -> cID -> vToCID[cID + N] == cID ? cID : clusterID = vToCID[cID + N]

    // Update VToCID based on cluster redirection
    for (int i = 0; i != N; ++i) {
        int clusterID = vToCID[i];
        if (clusterID == -1) {
            vToCID[i] = clusters.size();
            clusters.push_back(vector<VertexID>(i));
        } else {
            while (vToCID[clusterID + N] != clusterID) {
                clusterID = vToCID[clusterID + N];
            }

            // path compression
            vToCID[vToCID[i] + N] = clusterID;

            vToCID[i] = clusterID;
        }
    }

    cout << "done union find" << endl;

    // the cluster ids may not be consecutive due to merging
    // Need to serialize the vToCID and consolidate cluster info
    int oldSize = clusters.size();
    int newSize = 0;
    vector<int> cIDToConsecutiveCID(oldSize, -1);
    for (int i = 0; i != N; ++i) {
        if (cIDToConsecutiveCID[vToCID[i]] == -1) {
            cIDToConsecutiveCID[vToCID[i]] = newSize;
            newSize ++;
        }
    }

    cout << "done serialization" << endl;
    cout << "old size: " << oldSize << endl;
    cout << "new size: " << newSize << endl;

    clusters.clear();
    vector<vector<VertexID>> temp(newSize, vector<VertexID>());
    clusters = temp;
    for (int i = 0; i != N; ++i) {
        vToCID[i] = cIDToConsecutiveCID[vToCID[i]];
        clusters[vToCID[i]].push_back(i);
    }

    cout << "done clustering" << endl;


    // for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
    //     vector<VertexID> cluster = clusters[i];
    //     cout << "Cluster " << i << ": ";
    //     for (int j = 0, sizeJ = cluster.size(); j != sizeJ; ++j) {
    //         cout << cluster[j] << " ";
    //     }
    //     cout << endl;
    // }

    // cout << "vToCID: " << endl;
    // for (int i = 0, sizeI = N; i != sizeI; ++i) {
    //     cout << i << " " << vToCID[i] << endl;
    // }

}

void DGraph::getNodesRequiredToUnBN(int index, pair<VertexID, vector<VertexID>>& entry) {
    unordered_set<VertexID>::iterator itr;
    unordered_set<VertexID> requiredNodeSet;
    vector<VertexID> requiredNodeVector;

    for (int j = 0, sizeJ = outE[index].size(); j != sizeJ; ++j) {
        // cout << i << " " << j << " " << outE[i][j].first << endl;
        requiredNodeSet.insert(outE[index][j].first);
    }
    for (int j = 0, sizeJ = inE[index].size(); j != sizeJ; ++j) {
        // cout << i << " " << j << " " << inE[i][j].first << endl;
        requiredNodeSet.insert(inE[index][j].first);
    }

    for (itr = requiredNodeSet.begin(); itr != requiredNodeSet.end(); itr++) {
        requiredNodeVector.push_back(*itr);
    }

    sort(requiredNodeVector.begin(), requiredNodeVector.end());

    entry = make_pair((VertexID) index, requiredNodeVector);

}

void getParentCluster(vector<int> parentCluster, int& currCID) {
    while (parentCluster[currCID] != -1) {
        cout << currCID << ", ";
        currCID = parentCluster[currCID];
    }
}

int findPos(int v, vector<VertexID> vs) {
    int head = 0;
    int tail = vs.size() - 1;
    int mid = (head + tail) / 2;

    while (head <= tail) {
        if ((int) vs[mid] < v) {
            head = mid + 1;
        } else if ((int) vs[mid] == v) {
            return mid;
        } else {
            tail = mid - 1;
        }

        mid = (head + tail) / 2;
    }

    if (head > tail) {
        return -1;
    }
}


void DGraph::minBoundaryNodesClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID, int maxClusterSize) {

    struct queueComparator {

        // if number of nodes in a is larger than that of b, a is after b
        bool operator()(pair<VertexID, vector<VertexID>> const& a, pair<VertexID, vector<VertexID>> const& b) const {
            return (a.second).size() < (b.second).size();
        }
    };



    priority_queue<pair<VertexID, vector<VertexID>>, vector<pair<VertexID, vector<VertexID>>>, queueComparator> requiredNodesQueue;
    vector<vector<VertexID>> requiredNodesVector;
    pair<VertexID, vector<VertexID>> entry;


    for (int i = 0; i != N; ++i) {
        getNodesRequiredToUnBN(i, entry);
        requiredNodesQueue.push(entry);
        requiredNodesVector.push_back(entry.second);
    }


    dynamic_bitset<> visited = dynamic_bitset<>(N);
    // stores the parent cluster of clusters[i]
    // vector<int> parentCluster;
    vector<int> clusterWeight;
    vector<VertexID> unmergeableNodes;
    int count = 0;

    while (count < N) {
        entry = requiredNodesQueue.top();
        requiredNodesQueue.pop();
		// the node
        int v1 = (int) entry.first;
		// a vector of required nodes
        vector<VertexID> v2s = entry.second;
        // cout << "entry: " << entry.first << endl;

        if (visited[v1] == 1) {
            continue;
        }
        visited[v1] = 1;
        count += 1;

        // check total weight
        int existingWeight = 0;
		// use set to remove duplicates
        unordered_set<int> cIDs;

        // v1 
		// v1 does not belong to any cluster
        if (vToCID[v1] == -1) {
            existingWeight += 1;
        } else {
            int currCID = vToCID[v1];
            // getParentCluster(parentCluster, currCID);
            // cout << "parent cluster: " << currCID << endl;
            // path compression
            // if (parentCluster[vToCID[v1]] != -1) {
            //     parentCluster[vToCID[v1]] = currCID;
            // }
            cIDs.insert(currCID);
        }

        // v2s
        for (int i = 0, sizeI = v2s.size(); i != sizeI; ++i) {
            int v2 = (int) v2s[i];
            // cout << "v2 & cID: " << v2 << " " << vToCID[v2] << endl;
            // unclusterd node has weight of 1
            if (vToCID[v2] == -1) {
                existingWeight ++;
                // cout << "updated weight: "  << existingWeight << endl;
            } else {
                // find parent cluster
                int currCID = vToCID[v2];
                // cout << "get CID: ";
                // getParentCluster(parentCluster, currCID);
                // cout << currCID << endl;
                // path compression
                // if (parentCluster[vToCID[v2]] != -1) {
                //     parentCluster[vToCID[v2]] = currCID;
                // }
                cIDs.insert(currCID);
                // get weight
            }
        }
        unordered_set<int>::iterator itr;
        for (itr = cIDs.begin(); itr != cIDs.end(); itr ++) {
            existingWeight += clusterWeight[*itr];
            // cout << "member cluster and updated weight: "  << *itr << " " << existingWeight << endl;
        }

        // cout << "total weight of entry " << entry.first << ": " << existingWeight << endl;

        // proceed to merging if total doesn't exceeds
        if (existingWeight < maxClusterSize) {
            int clusterID;

            // create a new cluster
            if (vToCID[v1] == -1) {
                clusterID = clusters.size();
                clusters.push_back(vector<VertexID> {(VertexID) v1});
                // parentCluster.push_back(-1);
                clusterWeight.push_back(existingWeight);
                vToCID[v1] = clusterID;
            }
            // merge to existing CID
            else {
                // path compression should have been done in weight checking
                // if (parentCluster[vToCID[v1]] != -1) {
                //     clusterID = parentCluster[vToCID[v1]];
                // } else {
                //     clusterID = vToCID[v1];
                // }
                clusterID = vToCID[v1];
                clusterWeight[clusterID] = existingWeight;
            }

            // cout << "the growing cluster: " << clusterID << endl;
            // cout << "cluster size: " << clusters.size() << endl;

            vector<VertexID> cluster = clusters[clusterID];

			// merge the required nodes into the cluster
            for (int i = 0, sizeI = v2s.size(); i != sizeI; ++i) {
                int v2 = (int) v2s[i];
                int oldCID = vToCID[v2];
                // unclusterd node
                if (oldCID == -1) {
                    vToCID[v2] = clusterID;
                    clusters[clusterID].push_back(v2);
                } else {
                    if (oldCID != clusterID) {
                        // parentCluster[oldCID] = clusterID;
                        clusters[clusterID].insert(clusters[clusterID].end(), clusters[oldCID].begin(), clusters[oldCID].end());
                        // sort(clusters[clusterID].begin(), clusters[clusterID].end()); // should not have overlapping nodes
                        clusterWeight[oldCID] = 0;
                        for (int j = 0, sizeJ = clusters[oldCID].size(); j != sizeJ; ++j) {
                            vToCID[clusters[oldCID][j]] = clusterID;
                        }
                    }
                }
            }

			// update required nodes
            for (int i = 0, sizeI = clusters[clusterID].size(); i != sizeI; ++i) {
                // cout << "i: " << i << endl;
                int v = (int) clusters[clusterID][i];
                if (visited[v] == 1) {
                    continue;
                }

                // enumerate all nodes in the same cluster
				// remove link if there was a link between nodes in the same cluster
                bool clusterNodesUpdated = false;
                for (int j = 0, sizeJ = clusters[clusterID].size(); j != sizeJ; ++j) {
                    if (i != j) {
                        int pos = findPos(clusters[clusterID][j], requiredNodesVector[v]);
                        if (pos != -1) {
                            requiredNodesVector[v].erase(requiredNodesVector[v].begin() + pos);
                            clusterNodesUpdated = true;
                        }
                    }
                }

                // for nodes of neighbours of unvisited nodes of the cluster, 
                // remove redundant edges if they have neighbours in the same cluster
                unordered_set<int> clusterIDs;
                for (int j = 0, sizeJ = requiredNodesVector[v].size(); j != sizeJ; ++j) {
                    // cout << "j: " << j << "-";
                    VertexID v3 = requiredNodesVector[v][j];
                    bool neighbourUpdated = false;
                    int k = 0;
                    while (k < requiredNodesVector[v3].size()) {
                        // cout << k << " ";
                        VertexID v4 = requiredNodesVector[v3][k];
                        int cID = vToCID[v4];
                        if (cID != -1) {
                            if (clusterIDs.find(cID) != clusterIDs.end()) {
                                requiredNodesVector[v3].erase(requiredNodesVector[v3].begin() + k);
                                neighbourUpdated = true;
                                k --;
                            } else {
                                clusterIDs.insert(cID);
                            }
                        }
                        k++;
                    }
                    // cout << endl;
					// push a new entry into the queue if the node is updated
                    if (neighbourUpdated) {
                        requiredNodesQueue.push(make_pair(v3, requiredNodesVector[v3]));
                    }
                }
				// push a new entry into the queue if the node is updated
                if (clusterNodesUpdated) {
                    requiredNodesQueue.push(make_pair(v, requiredNodesVector[v]));
                }
            }
            

            // for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
            //     cout << "Cluster " << i << ": ";
            //     for (int j = 0, sizeJ = clusters[i].size(); j != sizeJ; ++j) {
            //         cout << clusters[i][j] << " ";
            //     }
            //     cout << endl;
            //     cout << "cluster weight: " << clusterWeight[i] << endl;
            // }

            // cout << endl;
        } else {
            // mark the node as not mergeable
            if (vToCID[v1] == -1) {
                unmergeableNodes.push_back(v1);
            }
        }
    }

    // for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
    //     cout << "Cluster " << i << ": ";
    //     for (int j = 0, sizeJ = clusters[i].size(); j != sizeJ; ++j) {
    //         cout << clusters[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    int oldSize = clusters.size();
    int newSize = 0;
    vector<int> cIDToConsecutiveCID(oldSize, -1);
    for (int i = 0; i != oldSize; ++i) {
        // int parentCID = i;
        // getParentCluster(parentCluster, parentCID);
        if (clusterWeight[i] != 0 && cIDToConsecutiveCID[i] == -1) {
            cIDToConsecutiveCID[i] = newSize;
            newSize ++;
        }
    }

    clusters.clear();
    for (int i = 0; i != newSize; ++i) {
        clusters.push_back(vector<VertexID>());
    }

    for (int i = 0; i < N; i++) {
        int cID = cIDToConsecutiveCID[vToCID[i]];
        clusters[cID].push_back(i);
        vToCID[i] = cID;
    }

	// merge single nodes to their smallest neighbouring nodes / clusters
	for (int i = 0, sizeI = unmergeableNodes.size(); i != sizeI; ++i) {
		// if the node is unclustered
		if (vToCID[(int) unmergeableNodes[i]] == -1) {
			int v1 = unmergeableNodes[i];
			int minNeighbourSize = N;
			int minV2 = -1;
			// for all its out neighbours
			for (int j = 0, sizeJ = outE[v1].size(); j != sizeJ; ++j) {
				int v2 = outE[v1][j].first;
				// if the out neighbour is a single node, 
				// immediately return this node as the target to be merged
				if (vToCID[v2] == -1) {
					minV2 = v2;
					minNeighbourSize = 1;
					break;
				}
				// otherwise update the smallest neighbour by its cluster size
				if (clusters[vToCID[v2]].size() < minNeighbourSize) {
					minNeighbourSize = clusters[vToCID[v2]].size();
					minV2 = v2;
				}
			}
			// if the node does not have any neighbours
			// make this node a new single node cluster
			if (minV2 == -1) {
				vToCID[v1] = clusters.size();
				clusters.push_back(vector<VertexID> {(VertexID) v1});
			} 
			// if the neighbour is a single node
			// make the two nodes a new two-node cluster
			else if (minNeighbourSize == 1) {
				vToCID[v1] = clusters.size();
				vToCID[minV2] = vToCID[v1];
				clusters.push_back(vector<VertexID> {(VertexID) v1, (VertexID) minV2});
			}
			// the neighbour is a cluster
			// add the node into the cluster
			else if (minNeighbourSize < maxClusterSize){
				vToCID[v1] = vToCID[minV2];
				clusters[vToCID[minV2]].push_back((VertexID) v1);
			}
		}
	}

	// Add unmergeable nodes into the clusters as single-node clusters
	for (int i = 0, sizeI = unmergeableNodes.size(); i != sizeI; ++i) {
		if (vToCID[(int)unmergeableNodes[i]] == -1) {
			vToCID[(int)unmergeableNodes[i]] = clusters.size();
			clusters.push_back(vector<VertexID> {unmergeableNodes[i]});
		}
	}

    // for (int i = 0, sizeI = clusters.size(); i != sizeI; ++i) {
    //     cout << "Cluster " << i << ": ";
    //     for (int j = 0, sizeJ = clusters[i].size(); j != sizeJ; ++j) {
    //         cout << clusters[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    //cout << "Done!" << endl;

}


void DGraph::initializeUnionFind(vector<VertexID>& parent) {
    for (int i = 0; i < N; i++) {
        parent.push_back((VertexID)i);
    }
}

int DGraph::find(vector<VertexID> parent, VertexID v) {
    if (parent[v] != v) {
        // path compression
        parent[v] = find(parent, parent[v]);
    }

    return parent[v];
}

void DGraph::connect(vector<VertexID>& parent, VertexID v, VertexID w) {
    int vSet = find(parent, v);
    int wSet = find(parent, w);
    parent[wSet] = vSet;
}


void DGraph::randomClustering(vector<vector<VertexID>>& clusters, vector<int>& vToCID) {
    // initialize all vertices to be in different clusters, where the Vertex ID is the Cluster ID
    cout << "No of Vertices: " << N << endl;
    cout << "No. of Edges: " << M << endl;
    vector<VertexID> parent;
    initializeUnionFind(parent);

    // randomly choose k edges
    int randomK = rand() % M;
    // currently we merge half of the edges
    // randomK = 20;
    cout << "Random number of edges to be merged: " << randomK << endl;

    for (int i = 0; i < randomK; i++) {
        int randomV = rand() % N;
        int numOfOutE = outE.at(randomV).size();
        while (numOfOutE == 0) {
            randomV = rand() % N;
            numOfOutE = outE.at(randomV).size();
        }
        int randomW = outE.at(randomV).at(rand() % numOfOutE).first;

        // cout << "Random edge with vertex IDs: " << randomV << ", " << randomW << endl;

        // link up parent of w to v
        connect(parent, randomV, randomW);
        if (i % 200 == 0) {
            cout << i << " ";
        }
    }
    cout << endl;

    int clusterCount = 0;
    map<int, int> ancestorToCID;

    for (int i = 0; i < N; i++) {
        int ancestor = find(parent, i);
        // Create a new cluster if the given ancestor is not in the map yet
        if (ancestorToCID.find(ancestor) == ancestorToCID.end()) {
            ancestorToCID[ancestor] = clusterCount;
            vToCID[i] = clusterCount;
            clusterCount ++;
            vector<VertexID> temp;
            temp.push_back((VertexID)i);
            clusters.push_back(temp);
        } else {
            clusters.at(ancestorToCID[ancestor]).push_back((VertexID)i);
            vToCID[i] = ancestorToCID[ancestor];
        }
    }

    // Print cluster info for verification
    // cout << "Clustering: " << endl;
    // for (int i = 0; i < clusterCount; i++) {
    //     cout << "Cluster " << i << ": ";
    //     vector<VertexID> currCluster = clusters.at(i);
    //     int cSize = currCluster.size();
    //     for (int j = 0; j < cSize; j++) {
    //         cout << currCluster.at(j) << ", ";
    //     }
    //     cout << endl;
    // }
}

void DGraph::tarjan(vector< vector<VertexID> >& SCCs) {
    int index = 0;
    stack<VertexID> q;
    vector< int > indexPerNode = vector< int >(N, -1);
    vector< int > lowlinkPerNode = vector< int >(N, -1);
    vector< bool > onStack = vector< bool >(N, false);

    for (int i = 0; i < N; i++) {
        if ( indexPerNode[i] == -1 ) {
            tarjanStrongConnect(i, index, q, indexPerNode, lowlinkPerNode, onStack, SCCs);
        }
    }
}

void DGraph::tarjanStrongConnect(int v, int& index, stack<VertexID>& q, vector< int >& indexPerNode,
                                 vector< int >& lowlinkPerNode, vector< bool >& onStack, vector< vector<VertexID> >& SCCs) {
    //cout << "v=" << v << ",index=" << index << endl;

    indexPerNode[v] = index;
    lowlinkPerNode[v] = index;
    index += 1;

    q.push( v );
    onStack[v] = true;

    SmallEdgeSet ses;
    getOutNeighbours(v, ses);
    for (int i = 0; i < ses.size(); i++) {
        VertexID w = ses[i].first;

        //cout << "v=" << v << ",w=" << w << ",indexPerNode[w]=" << indexPerNode[w] << ",onStack[w]=" << onStack[w] << ",q.size()=" << q.size() << endl;

        if ( indexPerNode[w] == -1 ) {
            tarjanStrongConnect(w, index, q, indexPerNode, lowlinkPerNode, onStack, SCCs);
            lowlinkPerNode[v] = min(lowlinkPerNode[v], lowlinkPerNode[w]);
        } else {
            if ( onStack[w] == true ) {
                lowlinkPerNode[v] = min(lowlinkPerNode[v], indexPerNode[w]);
            }
        }
    }

    if ( lowlinkPerNode[v] == indexPerNode[v] ) {
        vector< VertexID > SCC;
        VertexID w = 0;

        do {
            w = q.top();
            q.pop();
            //cout << "SCC w = " << w << endl;

            onStack[w] = false;
            SCC.push_back(w);
        } while ( w != v && q.empty() == false );

        SCCs.push_back(SCC);
    }
};

std::string DGraph::getStats() {
    string s = "";
    s += "N, " + to_string(N) + "\n";
    s += "M, " + to_string(M) + "\n";
    s += "L, " + to_string(L) + "\n";

    int noOfTriangles = computerNumberOfTriangles();
    double ratio = computeClusterCoefficient();
    s += "Number of triangles, " + to_string(noOfTriangles) + "\n";
    s += "Average cluster coefficient, " + to_string(ratio) + "\n";

    cout << "getStats: noOfTriangles=" << noOfTriangles << ",ratio=" << ratio << endl;

    vector< vector<VertexID> > SCCs;
    tarjan(SCCs);
    int maxSCC = 0;
    for (int i = 0; i < SCCs.size(); i++) {
        if ( SCCs[i].size() > maxSCC ) {
            maxSCC = SCCs[i].size();
        }
    }

    ratio = ((double) maxSCC) / ((double) N);
    cout << "getStats: SCCs.size()=" << SCCs.size() << ",ratio=" << ratio << endl;

    int diameter = computeDiameter();

    s += "Largest SCC #nodes/N, " + to_string(ratio) + "\n";
    s += "Number of SCCs, " + to_string(SCCs.size()) + "\n";
    s += "Diameter, " + to_string(diameter) + "\n";

    cout << "getStats: diameter=" << diameter << endl;

    s += "\n";
    s += "vertex id, out degree, in degree \n";

    vector< int > outFreqs;
    vector< int > inFreqs;
    vector< int > labelFreqs = vector<int>(L, 0);

    // in- and outdegree per vertex
    // and frequency of in- and outdegree
    // and label frequencies
    for (int i = 0; i < N; i++) {
        SmallEdgeSet outSes;
        getOutNeighbours(i, outSes);
        SmallEdgeSet inSes;
        getInNeighbours(i, inSes);
        int outDegree = outSes.size();
        int inDegree = inSes.size();

        //cout << to_string(i) + "," + to_string(outDegree) + "," + to_string(inDegree) + "\n";
        while ( outFreqs.size() < outDegree + 1 ) {
            outFreqs.push_back(0);
        }
        outFreqs[ outDegree ]++;

        while ( inFreqs.size() < inDegree + 1 ) {
            inFreqs.push_back(0);
        }
        inFreqs[ inDegree ]++;

        for (int j = 0; j < outSes.size(); j++) {
            LabelID lID = labelSetToLabelID(outSes[j].second);
            labelFreqs[ lID ]++;
        }

        s += to_string(i) + "," + to_string(outDegree) + "," + to_string(inDegree) + "\n";
    }

    // add frequencies to table
    s += "\n\noutdegree, frequency \n";
    for (int i = 0; i < outFreqs.size(); i++) {
        s += to_string(i) + "," + to_string(outFreqs[i]) + "\n";
    }

    s += "\n\nindegree, frequency \n";
    for (int i = 0; i < inFreqs.size(); i++) {
        s += to_string(i) + "," + to_string(inFreqs[i]) + "\n";
    }

    s += "\n";
    for (int i = 0; i < L; i++) {
        s += to_string(i);
        if ( i < L - 1 ) {
            s += ",";
        }
    }
    s += "\n";

    for (int i = 0; i < L; i++) {
        s += to_string( labelFreqs[i] );
        if ( i < L - 1 ) {
            s += ",";
        }
    }
    s += "\n";

    cout << "getStats: stats generated" << endl;

    return s;
};

int DGraph::computeDiameter() {
    // Gives an approximate diameter by taking N/SAMPLE_SIZE nodes
    int SAMPLE_SIZE = min(200, N / 10);
    vector<int> distances = vector<int>(N);
    int diameter = 0;

    for (int i = 0; i < N; i++) {
        if ( i % SAMPLE_SIZE != 0 && N >= 100 )
            continue;

        int newDiameter = findLongestShortestPath(i);
        if ( newDiameter > diameter ) {
            //cout << "i=" << i << " ,diameter=" << diameter << endl;
            diameter = newDiameter;
        }


    }

    return diameter + 1;
};
