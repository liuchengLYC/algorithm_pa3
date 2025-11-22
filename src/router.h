// router.h
#ifndef ROUTER_H
#define ROUTER_H

#include "graph.h"
#include "grid.h"
#include "types.h"
#include <vector>

struct Segment {
    Coord3D from;
    Coord3D to;
};

struct RoutedNet {
    std::string name;
    std::vector<Segment> segments; // continuous path in order
};

struct RoutingResult {
    std::vector<RoutedNet> nets;
};

class Router {
  public:
    Router();
    explicit Router(const Grid &grid);
    /// Build a graph from the current grid, using preferred directions.
    /// Each GCell becomes one vertex; edges represent allowed moves.
    Graph buildGraphFromGrid(const Grid &grid);

    /// Compute per-vertex congestion cost based on grid.demand().
    /// Students can modify this function to experiment with other cost models.
    void computeVertexCost(const Grid &grid, std::vector<int> &costs);

    /// Run routing for all nets.  Students will mainly implement this.
    RoutingResult runRouting(Grid &grid, const std::vector<Net> &nets);

    bool writeRouteFile(const std::string &filename,
                        const RoutingResult &result);

    /// A minimal Dijkstra interface that students can call or modify.
    /// They can also write their own version if they prefer.
    void dijkstra(const Graph &g, int source, int target);

  private:
    std::vector<int> dist;
    std::vector<int> prev;
    std::vector<int> stamp;
    std::vector<int> costs;
};
#endif // ROUTER_H
