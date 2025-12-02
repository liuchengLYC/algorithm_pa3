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
    void computeVertexCost(const Grid &grid);

    /// Run routing for all nets.  Students will mainly implement this.
    RoutingResult runRouting(Grid &grid, const std::vector<Net> &nets);

    bool writeRouteFile(const std::string &filename,
                        const RoutingResult &result);

    /// A minimal Dijkstra interface that students can call or modify.
    /// They can also write their own version if they prefer.
    void dijkstra(const Graph &g, int source, int target);

    std::pair<int, int> compute_gcell_overflow(const Grid &grid);

    /// compute all net scores to determine who to be ripped-up
    void compute_all_net_scores(const Grid &grid, const size_t netsize);

    /// @param topk topk != 0: select base on top k nets -> return vector size =
    /// k
    /// @param threshold threshold != 0: select nets with overflow greater than
    /// threshold
    std::vector<std::pair<int, int>> select_ripup_nets(int topk, int threshold);

    void rip_up_and_reroute(Grid &grid, int netId, const Net &net,
                            RoutedNet &routed, Graph &graph);

  private:
    // stamp is quite useless now
    std::vector<int> dist;
    std::vector<int> prev;
    std::vector<int> stamp;
    std::vector<int> costs;
    std::vector<int> cell_overflow;
    std::vector<std::pair<int, int>> cell_score;
    std::vector<std::vector<Coord3D>> net_path;
};
#endif // ROUTER_H
