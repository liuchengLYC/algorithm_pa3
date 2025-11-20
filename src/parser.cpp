// parser.cpp
#include "parser.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <vector>

namespace {

std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
        ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
        --e;
    return s.substr(b, e - b);
}

bool parseCoordLine(const std::string &line, Coord3D &coord) {
    std::istringstream iss(line);
    char ch;
    if (!(iss >> ch) || ch != '(')
        return false;
    if (!(iss >> coord.layer))
        return false;
    if (!(iss >> ch) || ch != ',')
        return false;
    if (!(iss >> coord.col))
        return false;
    if (!(iss >> ch) || ch != ',')
        return false;
    if (!(iss >> coord.row))
        return false;
    if (!(iss >> ch) || ch != ')')
        return false;
    return true;
}

bool readNonEmptyLine(std::ifstream &fin, std::string &out) {
    while (std::getline(fin, out)) {
        std::string t = trim(out);
        if (!t.empty()) {
            out = t;
            return true;
        }
    }
    return false;
}

} // namespace

bool parseCapFile(const std::string &filename, Grid &grid) {
    std::ifstream fin(filename);
    if (!fin)
        return false;

    int numLayers = 0;
    int xSize = 0;
    int ySize = 0;
    if (!(fin >> numLayers >> xSize >> ySize))
        return false;
    if (numLayers != 2)
        return false;

    grid.resize(xSize, ySize);

    int viaCost = 0;
    if (!(fin >> viaCost))
        return false;
    grid.setViaCost(viaCost);

    std::vector<int> horizontal(std::max(0, xSize - 1));
    for (int j = 0; j < xSize - 1; ++j) {
        if (!(fin >> horizontal[j]))
            return false;
    }
    grid.setHorizontalDistances(horizontal);

    std::vector<int> vertical(std::max(0, ySize - 1));
    for (int i = 0; i < ySize - 1; ++i) {
        if (!(fin >> vertical[i]))
            return false;
    }
    grid.setVerticalDistances(vertical);

    for (int l = 0; l < numLayers; ++l) {
        std::string layerName;
        std::string dir;
        if (!(fin >> layerName >> dir))
            return false;
        LayerInfo info{layerName, dir.empty() ? 'H' : dir[0]};
        grid.setLayerInfo(l, info);

        for (int row = 0; row < ySize; ++row) {
            for (int col = 0; col < xSize; ++col) {
                int cap = 0;
                if (!(fin >> cap))
                    return false;
                grid.setCapacity(l, col, row, cap);
            }
        }
    }

    grid.resetDemand();
    return true;
}

bool parseNetFile(const std::string &filename, std::vector<Net> &nets) {
    std::ifstream fin(filename);
    if (!fin)
        return false;

    std::string line;
    while (readNonEmptyLine(fin, line)) {
        if (!line.empty() && line[0] == '#')
            continue;

        Net net;
        net.name = line;

        if (!readNonEmptyLine(fin, line) || line.find('(') == std::string::npos)
            return false;

        Coord3D pins[2];
        for (int idx = 0; idx < 2; ++idx) {
            if (!readNonEmptyLine(fin, line))
                return false;
            if (!parseCoordLine(line, pins[idx]))
                return false;
        }

        if (!readNonEmptyLine(fin, line) || line.find(')') == std::string::npos)
            return false;

        net.pin1 = pins[0];
        net.pin2 = pins[1];
        nets.push_back(net);
    }

    return true;
}

bool parseInputFiles(const std::string &capFilename,
                     const std::string &netFilename, ParsedInput &out) {
    Grid grid;
    std::vector<Net> nets;

    if (!parseCapFile(capFilename, grid))
        return false;
    if (!parseNetFile(netFilename, nets))
        return false;

    out.grid = std::move(grid);
    out.nets = std::move(nets);
    return true;
}
