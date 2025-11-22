// grid.h
#ifndef GRID_H
#define GRID_H

#include "types.h"
#include <string>
#include <vector>

struct LayerInfo {
    std::string name; // e.g., "Metal1"
    char direction;   // 'H' or 'V'
};

class Grid {
  public:
    Grid();

    // dimensions
    int numLayers() const { return 2; }
    int xSize() const { return xSize_; }
    int ySize() const { return ySize_; }
    int gridSize() const { return gridlSize_; }

    // index mapping between (l, j, i) and a flat GCell index
    int gcellIndex(int l, int j, int i) const;
    Coord3D fromIndex(int idx) const;

    // capacity & demand on GCells
    int capacity(int l, int j, int i) const;
    int capacity(Coord3D c) const;
    void setCapacity(int l, int j, int i, int cap);

    int demand(int l, int j, int i) const;
    int demand(Coord3D c) const;
    void resetDemand();
    void addDemandForNetGCell(int netId, int l, int j, int i);
    int demandByIndex(int idx) const;

    // distances / via cost
    int wlViaCost() const { return wlViaCost_; }
    void setViaCost(int cost) { wlViaCost_ = cost; }
    int horizontalDist(int j) const; // W_j
    int verticalDist(int i) const;   // H_i
    void setHorizontalDistances(const std::vector<int> &distances);
    void setVerticalDistances(const std::vector<int> &distances);

    const LayerInfo &layerInfo(int l) const { return layers_[l]; }
    void setLayerInfo(int l, const LayerInfo &info);

    void resize(int xSize, int ySize);

  private:
    int xSize_ = 0;
    int ySize_ = 0;
    int gridlSize_ = 0;

    // size: [2 * xSize * ySize]
    std::vector<int> capacity_;
    std::vector<int> demand_;

    std::vector<int> W_; // size xSize-1
    std::vector<int> H_; // size ySize-1

    int wlViaCost_ = 0;
    LayerInfo layers_[2];
};

#endif // GRID_H
