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
    adj_[u].push_back({v, baseCost});
    adj_[v].push_back({u, baseCost});
}
