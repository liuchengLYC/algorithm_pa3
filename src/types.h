// types.h
#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <vector>

struct Coord3D {
    int layer; // 0 or 1
    int col;   // j
    int row;   // i
};

struct Net {
    std::string name;
    Coord3D pin1;
    Coord3D pin2;
};

#endif // TYPES_H
