// graph.h
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <limits>

// Simple directed edge for adjacency list
struct Edge {
    int to;          // destination vertex id
    int baseCost;    // physical length (W_j, H_i, or via)
    // students can ignore or extend this structure as needed
};

class Graph {
public:
    Graph();
    explicit Graph(int numVertices);

    void resize(int numVertices);
    int numVertices() const { return static_cast<int>(adj_.size()); }

    void addEdge(int u, int v, int baseCost);

    const std::vector<Edge>& adj(int u) const { return adj_[u]; }

private:
    std::vector<std::vector<Edge>> adj_;
};

/// A minimal Dijkstra interface that students can call or modify.
/// They can also write their own version if they prefer.
std::vector<int> dijkstra(
    const Graph &g,
    int source,
    std::vector<int> *outPrev = nullptr   // optional predecessor tree
);

const int INF = std::numeric_limits<int>::max() / 4;

#endif // GRAPH_H
