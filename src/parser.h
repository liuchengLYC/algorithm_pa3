// parser.h
#ifndef PARSER_H
#define PARSER_H

#include "grid.h"
#include "types.h"
#include <string>
#include <vector>

struct ParsedInput {
    Grid grid;
    std::vector<Net> nets;
};

bool parseCapFile(const std::string &filename, Grid &grid);
bool parseNetFile(const std::string &filename, std::vector<Net> &nets);

/// Convenience wrapper: parse both files.
bool parseInputFiles(
    const std::string &capFilename,
    const std::string &netFilename,
    ParsedInput &out
);

#endif // PARSER_H
