// router.cpp
#include "router.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>

namespace {

constexpr int OVERFLOW_WEIGHT = 100000;

int totalVertices(const Grid &grid) {
    return grid.numLayers() * grid.xSize() * grid.ySize();
}

char normalizeDir(char d) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(d)));
}

std::vector<Coord3D> reconstructPath(
    const Grid &grid,
    int sourceIdx,
    int targetIdx,
    const std::vector<int> &prev
) {
    // TODO: Starting from targetIdx, follow the `prev` parent array (produced by
    // dijkstra) until you reach sourceIdx or hit -1. You should store each vertex
    // index you visit so you can reverse the order later. Use grid.fromIndex(idx)
    // to convert each flat index back to a Coord3D {layer, col, row}. The resulting
    // vector should go from the source coordinate to the target coordinate.
    // Return an empty vector if the target is unreachable (i.e., prev chain does
    // not lead back to the source).
    return {};
}

std::vector<Coord3D> buildFallbackPath(const Grid &grid, const Net &net) {
    // TODO: Produce a deterministic Manhattan-style backup route when Dijkstra fails. 
    // Start from net.pin1 (source) and peek at grid.layerInfo(layer).direction to decide whether you need to swap layers
    // so that the preferred direction (H or V) matches the axis you are about to walk. 
    // Move one gcell at a time by incrementing/decrementing Coord3D::col for horizontal moves and Coord3D::row for vertical moves until you reach net.pin2.
    // Whenever you change layer, append the intermediate Coord3D so downstream code can emit via segments.
    // Return the ordered list of coordinates from source to target.
    return {};
}

void updateDemandAlongPath(
    Grid &grid,
    int netId,
    const std::vector<Coord3D> &coords
) {
    for (const Coord3D &c : coords) {
        grid.addDemandForNetGCell(netId, c.layer, c.col, c.row);
    }
}

} // namespace

Graph buildGraphFromGrid(const Grid &grid) {
    Graph g(totalVertices(grid));

    // TODO: Build the routing graph by iterating over layers/rows/cols and calling
    // grid.gcellIndex(layer, col, row) to convert 3D coordinates into vertex ids.
    // Use grid.layerInfo(layer).direction to determine whether the preferred track direction is horizontal or vertical, 
    // and add edges with g.addEdge(u, v, cost) using the wire length returned by grid.horizontalDist(col) / grid.verticalDist(row).
    // Remember to add via edges between layers when grid.numLayers() > 1, using grid.wlViaCost() as the edge weight.

    return g;
}

std::vector<int> computeVertexCost(const Grid &grid) {
    const int total = totalVertices(grid);
    std::vector<int> costs(total, 0);

    // TODO: Translate the 1D vertex index back to (layer, col, row) with grid.fromIndex(idx) 
    // and penalize vertices whose demand would overflow their capacity. 
    // Use grid.demand() / grid.capacity() and OVERFLOW_WEIGHT as the per-unit overflow penalty.

    return costs;
}

RoutingResult runRouting(
    Grid &grid,
    const std::vector<Net> &nets
) {
    RoutingResult result;
    result.nets.reserve(nets.size());

    grid.resetDemand();
    Graph graph = buildGraphFromGrid(grid);

    for (size_t netIdx = 0; netIdx < nets.size(); ++netIdx) {
        const Net &net = nets[netIdx];
        RoutedNet routed;
        routed.name = net.name;

        int src = grid.gcellIndex(net.pin1.layer, net.pin1.col, net.pin1.row);
        int dst = grid.gcellIndex(net.pin2.layer, net.pin2.col, net.pin2.row);

        std::vector<int> predecessors;
        std::vector<int> costs = computeVertexCost(grid);
        // TODO: integrate the per-vertex penalty into your search.

        std::vector<Coord3D> path;
        // TODO: Run Dijkstra on `graph` from src -> dst (see graph.h for Graph usage).
        // Use reconstructPath(grid, src, dst, predecessors) once you have the parent
        // array produced by dijkstra(...) and only fall back if no legal path exists.

        if (path.empty()) {
            std::cerr << "Warning: using fallback routing for " << net.name << "\n";
            path = buildFallbackPath(grid, net);
        }

        for (size_t i = 1; i < path.size(); ++i) {
            routed.segments.push_back(Segment{path[i - 1], path[i]});
        }
        result.nets.push_back(routed);

        // TODO: After a successful route, call updateDemandAlongPath() so that subsequent nets see the updated usage when computeVertexCost(...) is rerun.
    }

    return result;
}

bool writeRouteFile(const std::string &filename, const RoutingResult &result) {
    std::ofstream fout(filename);
    if (!fout) return false;

    for (const RoutedNet &net : result.nets) {
        fout << net.name << "\n";
        fout << "(\n";
        for (const Segment &seg : net.segments) {
            fout << seg.from.layer << " " << seg.from.col << " " << seg.from.row << " "
                 << seg.to.layer << " " << seg.to.col << " " << seg.to.row << "\n";
        }
        fout << ")\n";
    }

    return true;
}
