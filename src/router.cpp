// router.cpp
#include "router.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <utility>
using namespace std;
const bool if_reroute = true;

#define debug false

namespace {

constexpr int OVERFLOW_WEIGHT = 300000;
constexpr double ALMOST_FULL_WEIGHT = 3000.0;
constexpr int HISTORY_INCREMENT = 4000;
constexpr double COST_BASE_POW = 2.0;

using Prefix2D = vector<vector<long long>>;

int totalVertices(const Grid &grid) { return grid.gridSize(); }

// char normalizeDir(char d) {
//     return static_cast<char>(toupper(static_cast<unsigned char>(d)));
// }

vector<Coord3D> reconstructPath(const Grid &grid, int src, int dst,
                                const vector<int> &prev) {
    vector<Coord3D> tmp;
    int cur = dst;
    while (cur != -1 && cur != src) {
        tmp.push_back(grid.fromIndex(cur));
        cur = prev[cur];
    }
    if (cur != src)
        return {};
    else {
        tmp.push_back(grid.fromIndex(src));
        reverse(tmp.begin(), tmp.end());
        return tmp;
    }
    // TODO: Starting from targetIdx, follow the `prev` parent array
    // (produced by dijkstra) until you reach sourceIdx or hit -1. You
    // should store each vertex index you visit so you can reverse the order
    // later. Use grid.fromIndex(idx) to convert each flat index back to a
    // Coord3D {layer, col, row}. The resulting vector should go from the
    // source coordinate to the target coordinate. Return an empty vector if
    // the target is unreachable (i.e., prev chain does not lead back to the
    // source).
}

vector<Coord3D> buildFallbackPath(const Grid &grid, const Net &net) {
    vector<Coord3D> ans;
    char dir = grid.layerInfo(net.pin1.layer).direction;
    int dx = (net.pin1.col < net.pin2.col) ? 1 : -1;
    int dy = (net.pin1.row < net.pin2.row) ? 1 : -1;

    Coord3D cur = net.pin1;
    ans.push_back(cur);

    if (dir == 'H') {
        while (cur.col != net.pin2.col) {
            cur.col += dx;
            ans.push_back(cur);
        }
        if (cur.row != net.pin2.row) {
            if (grid.layerInfo(cur.layer).direction != 'V') {
                cur.layer ^= 1;
                ans.push_back(cur);
            }
            while (cur.row != net.pin2.row) {
                cur.row += dy;
                ans.push_back(cur);
            }
        }
        if (cur.layer != net.pin2.layer) {
            cur.layer ^= 1;
            ans.push_back(cur);
        }

    } else {
        while (cur.row != net.pin2.row) {
            cur.row += dy;
            ans.push_back(cur);
        }

        if (cur.col != net.pin2.col) {
            if (grid.layerInfo(cur.layer).direction != 'H') {
                cur.layer ^= 1;
                ans.push_back(cur);
            }
            while (cur.col != net.pin2.col) {
                cur.col += dx;
                ans.push_back(cur);
            }
        }
        if (cur.layer != net.pin2.layer) {
            cur.layer ^= 1;
            ans.push_back(cur);
        }
    }
    return ans;
}

void updateDemandAlongPath(Grid &grid, int netId,
                           const vector<Coord3D> &coords) {
    grid.beginAddDemandForNet(netId);
    for (const Coord3D &c : coords) {
        grid.addDemandForNetGCell(netId, c.layer, c.col, c.row);
    }
}

void removeDemandAlongPath(Grid &grid, int netId,
                           const vector<Coord3D> &coords) {
    grid.beginRemoveDemandForNet(netId);
    for (const Coord3D &c : coords) {
        grid.removeDemandForNetGCell(netId, c.layer, c.col, c.row);
    }
}

void path_to_seg(vector<Coord3D> &tmp) {
    if (tmp.size() <= 2)
        return;
    vector<Coord3D> ans = {tmp[0]};
    for (size_t i = 1; i < tmp.size() - 1; i++) {
        if ((tmp[i - 1].col == tmp[i + 1].col &&
             tmp[i - 1].row == tmp[i + 1].row) ||
            (tmp[i - 1].row == tmp[i + 1].row &&
             tmp[i - 1].layer == tmp[i + 1].layer) ||
            (tmp[i - 1].layer == tmp[i + 1].layer &&
             tmp[i - 1].col == tmp[i + 1].col))
            continue;
        else
            ans.push_back(tmp[i]);
    }
    ans.push_back(tmp.back());
    tmp = ans;
}

Prefix2D build_capacity_prefix(const Grid &grid) {
    const int X = grid.xSize();
    const int Y = grid.ySize();
    const int L = grid.numLayers();

    Prefix2D ps(Y + 1, vector<long long>(X + 1, 0));

    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            int capSum = 0;
            for (int l = 0; l < L; ++l) {
                capSum += grid.capacity(l, j, i);
            }
            ps[i + 1][j + 1] = ps[i + 1][j] + ps[i][j + 1] - ps[i][j] + capSum;
        }
    }
    return ps;
}

double bbox_avg_capacity(const Prefix2D &ps, int xmin, int xmax, int ymin,
                         int ymax) {
    if (xmax < xmin || ymax < ymin)
        return 0.0;

    long long sum = ps[ymax + 1][xmax + 1] - ps[ymin][xmax + 1] -
                    ps[ymax + 1][xmin] + ps[ymin][xmin];

    long long area = 1LL * (xmax - xmin + 1) * (ymax - ymin + 1);
    if (area <= 0)
        return 0.0;

    return static_cast<double>(sum) / static_cast<double>(area);
}

pair<double, double> netSortKey(const Prefix2D &ps, const Net &n) {
    const double ratio = 1000.0;
    int dx = abs(n.pin1.col - n.pin2.col);
    int dy = abs(n.pin1.row - n.pin2.row);
    int dz = abs(n.pin1.layer - n.pin2.layer);
    double man_dist = dx + dy + dz;

    int xmin = min(n.pin1.col, n.pin2.col);
    int xmax = max(n.pin1.col, n.pin2.col);
    int ymin = min(n.pin1.row, n.pin2.row);
    int ymax = max(n.pin1.row, n.pin2.row);

    double avgCap = bbox_avg_capacity(ps, xmin, xmax, ymin, ymax);
    if (avgCap <= 0.0)
        avgCap = 1.0;

    double difficulty = 1.0 / avgCap;

    return make_pair(man_dist, man_dist + ratio * difficulty);
}

} // namespace

Router::Router() = default;

Router::Router(const Grid &grid) {
    const int total = totalVertices(grid);
    dist.assign(total, INF);
    prev.assign(total, -1);
    costs.assign(total, INF);
    cell_overflow.assign(total, 0);
    history_costs.assign(total, 0);
}

Graph Router::buildGraphFromGrid(const Grid &grid) {
    Graph g(totalVertices(grid));
    auto gi = [&](int x, int y, int z) { return grid.gcellIndex(x, y, z); };
    for (int lay = 0; lay < grid.numLayers(); lay++) {
        for (int j = 0; j < grid.xSize(); j++) {
            for (int i = 0; i < grid.ySize(); i++) {
                if (lay > 0)
                    g.addEdge(gi(lay - 1, j, i), gi(lay, j, i),
                              grid.wlViaCost());
                if (grid.layerInfo(lay).direction == 'H' &&
                    i < grid.ySize() - 1 && j < grid.xSize() - 1)
                    g.addEdge(gi(lay, j, i), gi(lay, j + 1, i),
                              grid.horizontalDist(j));
                else if (i < grid.ySize() - 1 && j < grid.xSize() - 1)
                    g.addEdge(gi(lay, j, i), gi(lay, j, i + 1),
                              grid.verticalDist(i));
            }
        }
    }
    return g;
}

void Router::computeVertexCost(const Grid &grid) {
    const int total = totalVertices(grid);
    auto cfunc = [&](int x) {
        int dm = grid.demandByIndex(x);
        int cap = grid.capacity(grid.fromIndex(x));
        int hist = history_costs[x];
        int base_cost = 1;

        if (dm >= cap) {
            int overflow = dm - cap;
            double congestion_penalty =
                OVERFLOW_WEIGHT * std::pow(COST_BASE_POW, overflow);
            return static_cast<int>(base_cost + hist + congestion_penalty);
        } else {
            double ratio_cost =
                (static_cast<double>(dm) / static_cast<double>(cap)) *
                ALMOST_FULL_WEIGHT;
            return static_cast<int>(base_cost + hist + ratio_cost);
        }
    };

    for (int i = 0; i < total; i++) {
        costs[i] = cfunc(i);
    }
}

void Router::dijkstra(const Graph &g, int source, int target) {
    // dist -> distance(or cost) in dijkstra, cost -> cost to step on gcell
    fill(dist.begin(), dist.end(), INF);
    fill(prev.begin(), prev.end(), -1);
    using pii = pair<int, int>;
    auto cmp = [&](pii p, pii q) { return p.second > q.second; };
    priority_queue<pii, vector<pii>, decltype(cmp)> pq(cmp);
    pq.push({source, costs[source]});
    dist[source] = costs[source];
    while (!pq.empty()) {
        auto [v, c] = pq.top();
        pq.pop();
        if (c != dist[v])
            continue;
        else if (v == target)
            return;
        for (auto &e : g.adj(v)) {
            int nc = c + costs[e.to] + e.baseCost;
            if (nc < dist[e.to]) {
                pq.push({e.to, nc});
                prev[e.to] = v;
                dist[e.to] = nc;
            }
        }
    }
}

pair<int, int> Router::compute_gcell_overflow(const Grid &grid) {
    const int total = totalVertices(grid);
    int total_of = 0, max_of = -1;
    for (int i = 0; i < total; i++) {
        int dm = grid.demandByIndex(i), cap = grid.capacity(grid.fromIndex(i));
        int of = max(dm - cap, 0);
        cell_overflow[i] = of;
        total_of += of;
        max_of = max(max_of, of);
    }
    return make_pair(total_of, max_of);
}

void Router::compute_all_net_scores(const Grid &grid, const size_t netsize) {
    // score function
    auto scorefunc = [&](int of) { return of; };
    if (!cell_score.size())
        cell_score.resize(netsize);
    for (size_t i = 0; i < netsize; i++) {
        int score = 0;
        for (auto &cell : net_path[i]) {
            score += scorefunc(
                cell_overflow[grid.gcellIndex(cell.layer, cell.col, cell.row)]);
        }
        cell_score[i] = make_pair(i, score);
    }
}

vector<pair<int, int>> Router::select_ripup_nets(int topk, int threshold) {
    if (topk != -1) {
        vector<pair<int, int>> tmp = cell_score;
        if (static_cast<size_t>(topk) > cell_score.size()) {
            cerr << "Topk is too large" << '\n';
            return tmp;
        }
        partial_sort(tmp.begin(), tmp.begin() + topk, tmp.end(),
                     [&](auto x, auto y) { return x.second > y.second; });
        return vector<pair<int, int>>(tmp.begin(), tmp.begin() + topk);
    } else if (threshold != -1) {
        vector<pair<int, int>> ret;
        copy_if(cell_score.begin(), cell_score.end(), back_inserter(ret),
                [&](auto x) { return x.second >= threshold; });
        // sort(ret.begin(), ret.end(),
        //      [&](auto x, auto y) { return x.second > y.second; });
        return ret;
    } else
        return {};
}

void Router::rip_up_and_reroute(Grid &grid, int netId, const Net &net,
                                RoutedNet &routed, Graph &graph) {
    auto path = net_path[netId];
    removeDemandAlongPath(grid, netId, path);
    path.clear();
    routed.segments.clear();

    int src = grid.gcellIndex(net.pin1.layer, net.pin1.col, net.pin1.row);
    int dst = grid.gcellIndex(net.pin2.layer, net.pin2.col, net.pin2.row);

    computeVertexCost(grid);
    dijkstra(graph, src, dst);
    path = reconstructPath(grid, src, dst, prev);

    if (path.empty()) {
        cerr << "Warning: using fallback routing for " << net.name << "\n";
        path = buildFallbackPath(grid, net);
    }

    net_path[netId] = path;
    updateDemandAlongPath(grid, netId, path);
    path_to_seg(path);

    for (size_t i = 1; i < path.size(); ++i) {
        routed.segments.push_back(Segment{path[i - 1], path[i]});
    }
}

RoutingResult Router::runRouting(Grid &grid, const vector<Net> &nets) {
    RoutingResult result;
    result.nets.reserve(nets.size());

    grid.resetDemand();
    fill(history_costs.begin(), history_costs.end(), 0);
    Graph graph = buildGraphFromGrid(grid);
    vector<Net> nnets = nets;
    Prefix2D ps = build_capacity_prefix(grid);
    net_path.resize(nets.size());

    sort(nnets.begin(), nnets.end(), [&](const Net &a, const Net &b) {
        auto ka = netSortKey(ps, a);
        auto kb = netSortKey(ps, b);
        // swap(ka.first, ka.second), swap(kb.first, kb.second);
        if (ka.first != kb.first)
            return ka.first < kb.first;
        else
            return ka.second > kb.second;
    });

    // all process is based on sorted nets(i.e. nnets) from now on
    for (size_t netIdx = 0; netIdx < nnets.size(); ++netIdx) {
        const Net &net = nnets[netIdx];
        RoutedNet routed;
        routed.name = net.name;

        int src = grid.gcellIndex(net.pin1.layer, net.pin1.col, net.pin1.row);
        int dst = grid.gcellIndex(net.pin2.layer, net.pin2.col, net.pin2.row);

        computeVertexCost(grid);
        dijkstra(graph, src, dst);

        vector<Coord3D> path = reconstructPath(grid, src, dst, prev);

        if (path.empty()) {
            cerr << "Warning: using fallback routing for " << net.name << "\n";
            path = buildFallbackPath(grid, net);
        }

        net_path[netIdx] = path;
        updateDemandAlongPath(grid, netIdx, path);
        path_to_seg(path);

        for (size_t i = 1; i < path.size(); ++i) {
            routed.segments.push_back(Segment{path[i - 1], path[i]});
        }
        result.nets.push_back(routed);
    }

    pair<int, int> stat = compute_gcell_overflow(grid);

    if (if_reroute) {
        for (int i = 0; i < 50; i++) {
            stat = compute_gcell_overflow(grid);
            if (stat.first == 0)
                break;

            for (size_t idx = 0; idx < cell_overflow.size(); ++idx) {
                if (cell_overflow[idx] > 0) {
                    history_costs[idx] += HISTORY_INCREMENT;
                }
            }

            compute_all_net_scores(grid, nnets.size());
            // vector<pair<int, int>> reroute_list =
            //     select_ripup_nets(nnets.size() / 10, -1);
            vector<pair<int, int>> reroute_list = select_ripup_nets(-1, 1);

#if debug
            cerr << i << " " << stat.first << " " << stat.second;
            cerr << " " << reroute_list.size() << '\n';
#endif
            for (auto &[netId, _] : reroute_list) {
                rip_up_and_reroute(grid, netId, nnets[netId],
                                   result.nets[netId], graph);
            }
        }
    }
#if debug
    stat = compute_gcell_overflow(grid);
    cerr << "After routing, total overflow is: " << stat.first << '\n';
#endif
    return result;
}

bool Router::writeRouteFile(const string &filename,
                            const RoutingResult &result) {
    ofstream fout(filename);
    if (!fout)
        return false;
    for (const RoutedNet &net : result.nets) {
        fout << net.name << "\n";
        fout << "(\n";
        for (const Segment &seg : net.segments) {
            fout << seg.from.layer << " " << seg.from.col << " " << seg.from.row
                 << " " << seg.to.layer << " " << seg.to.col << " "
                 << seg.to.row << "\n";
        }
        fout << ")\n";
    }

    return true;
}