// graph.cpp
#include "graph.h"
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
using namespace std;

Graph::Graph() = default;

Graph::Graph(int numVertices) { resize(numVertices); }

void Graph::resize(int numVertices) { adj_.assign(numVertices, {}); }

void Graph::addEdge(int u, int v, int baseCost) {
    // TODO: add edge (u, v) between GCell u and v with baseCosts to the graph
}

std::vector<int> dijkstra(const Graph &g, int source, vector<int> *outPrev) {
    const int n = g.numVertices();
    vector<int> dist(n, INF);
    vector<int> prev(n, -1);

    // TODO: implement Dijkstra's algorithm to compute shortest paths from
    // source

    return dist;
}
