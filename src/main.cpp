// main.cpp
#include "parser.h"
#include "router.h"
#include <iostream>

int main(int argc, char **argv) {
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0]
                  << " --cap case.cap --net case.net --out case.route\n";
        return 1;
    }

    std::string capFile, netFile, outFile;
    // simple argument parsing (TA implemented)
    for (int i = 1; i < argc; i += 2) {
        std::string opt = argv[i];
        if (opt == "--cap")
            capFile = argv[i + 1];
        else if (opt == "--net")
            netFile = argv[i + 1];
        else if (opt == "--out")
            outFile = argv[i + 1];
    }

    ParsedInput input;
    if (!parseInputFiles(capFile, netFile, input)) {
        std::cerr << "Error: failed to parse input files.\n";
        return 1;
    }

    Router router;
    RoutingResult result = router.runRouting(input.grid, input.nets);

    if (!router.writeRouteFile(outFile, result)) {
        std::cerr << "Error: failed to write route file.\n";
        return 1;
    }

    return 0;
}
