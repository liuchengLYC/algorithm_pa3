// grid.cpp
#include "grid.h"
#include <algorithm>
#include <stdexcept>

Grid::Grid() = default;

void Grid::resize(int xSize, int ySize) {
    xSize_ = xSize;
    ySize_ = ySize;
    const int total = numLayers() * xSize_ * ySize_;
    gridlSize_ = total;
    capacity_.assign(total, 0);
    demand_.assign(total, 0);
    W_.assign(xSize_ > 0 ? xSize_ - 1 : 0, 0);
    H_.assign(ySize_ > 0 ? ySize_ - 1 : 0, 0);
}

int Grid::gcellIndex(int l, int j, int i) const {
    if (l < 0 || l >= numLayers() || j < 0 || j >= xSize_ || i < 0 ||
        i >= ySize_)
        throw std::out_of_range("Invalid gcell coordinate");
    return l * (xSize_ * ySize_) + i * xSize_ + j;
}

Coord3D Grid::fromIndex(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(capacity_.size()))
        throw std::out_of_range("Invalid gcell index");
    int perLayer = xSize_ * ySize_;
    Coord3D c;
    c.layer = idx / perLayer;
    int rem = idx % perLayer;
    c.row = rem / xSize_;
    c.col = rem % xSize_;
    return c;
}

int Grid::capacity(int l, int j, int i) const {
    return capacity_[gcellIndex(l, j, i)];
}

void Grid::setCapacity(int l, int j, int i, int cap) {
    capacity_[gcellIndex(l, j, i)] = cap;
}

int Grid::demand(int l, int j, int i) const {
    return demand_[gcellIndex(l, j, i)];
}

void Grid::resetDemand() { std::fill(demand_.begin(), demand_.end(), 0); }

void Grid::addDemandForNetGCell(int /*netId*/, int l, int j, int i) {
    ++demand_[gcellIndex(l, j, i)];
}

int Grid::demandByIndex(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(demand_.size()))
        throw std::out_of_range("Invalid gcell index");
    return demand_[idx];
}

int Grid::horizontalDist(int j) const {
    if (j < 0 || j >= static_cast<int>(W_.size()))
        throw std::out_of_range("Invalid horizontal distance index");
    return W_[j];
}

int Grid::verticalDist(int i) const {
    if (i < 0 || i >= static_cast<int>(H_.size()))
        throw std::out_of_range("Invalid vertical distance index");
    return H_[i];
}

void Grid::setHorizontalDistances(const std::vector<int> &distances) {
    if (static_cast<int>(distances.size()) != (xSize_ > 0 ? xSize_ - 1 : 0))
        throw std::runtime_error("Horizontal distance vector size mismatch");
    W_ = distances;
}

void Grid::setVerticalDistances(const std::vector<int> &distances) {
    if (static_cast<int>(distances.size()) != (ySize_ > 0 ? ySize_ - 1 : 0))
        throw std::runtime_error("Vertical distance vector size mismatch");
    H_ = distances;
}

void Grid::setLayerInfo(int l, const LayerInfo &info) {
    if (l < 0 || l >= numLayers())
        throw std::out_of_range("Invalid layer index");
    layers_[l] = info;
}
