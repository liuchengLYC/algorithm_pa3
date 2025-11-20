// router.h
#ifndef ROUTER_H
#define ROUTER_H

#include "grid.h"
#include "graph.h"
#include "types.h"
#include <vector>

struct Segment {
    Coord3D from;
    Coord3D to;
};

struct RoutedNet {
    std::string name;
    std::vector<Segment> segments;  // continuous path in order
};

struct RoutingResult {
    std::vector<RoutedNet> nets;
};

/// Build a graph from the current grid, using preferred directions.
/// Each GCell becomes one vertex; edges represent allowed moves.
Graph buildGraphFromGrid(const Grid &grid);

/// Compute per-vertex congestion cost based on grid.demand().
/// Students can modify this function to experiment with other cost models.
std::vector<int> computeVertexCost(const Grid &grid);

/// Run routing for all nets.  Students will mainly implement this.
RoutingResult runRouting(
    Grid &grid,
    const std::vector<Net> &nets
);

bool writeRouteFile(const std::string &filename, const RoutingResult &result);

#endif // ROUTER_H
