#include <fstream>
#include <iostream>
#include "Grid.h"
#include "RockProperties.h"
#include "FluidProperties.h"
#include "TransmissibilityCalculator.h"
#include "Simulator.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main(int argc, char** argv){
    //look for config.json in the build folder
    const std::string cfgFile = (argc > 1) ? argv[1] : "config.json";

    std::ifstream f(cfgFile);
    if (!f) {
        std::cerr << "Cannot open config file: " << cfgFile << '\n';
        return 1;
    }

    json cfg;
    try {
        cfg = json::parse(f);
    } catch (json::parse_error& e) {
        std::cerr << "JSON parse error in " << cfgFile << ": " << e.what() << '\n';
        return 1;
    }

    /*//Error catching for json components Grid check
    if (cfg.contains("grid")) {
        const auto& g = cfg["grid"];
        std::cout << "Grid nx = " << g["nx"] << ", ny = " << g["ny"] << "\n";
        // *-----Grid can be build here-----*       !! Add check for value type !!
    } else {
        std::cerr << "Config missing 'grid' section!\n";
        return 1;
    }*/

    std::cout << "Configuration loaded successfully!\n";

    //Build grid
    const auto& g = cfg["grid"];
    Grid grid(g["nx"], g["ny"], g["nz"], g["dx"], g["dy"], g["dz"], g["depth"]);

    //Build rock properties
    const auto& rp = cfg["rock"];
    RockProperties rock(grid.numCells(),
                        rp["porosity"],
                        {rp["permeability"][0], rp["permeability"][1], rp["permeability"][2]});

    //Fluid properties
    const auto& fp = cfg["fluid"];
    double API = fp["API"];
    double gamma_o = 141.5 / (131.5 + API);
    double rhoOSC = fp["rho_o"];

    FluidProperties fluid(fp["Bo"], fp["Bg"],
                          fp["mu_o"], fp["mu_g"], fp["mu_w"],
                          fp["rho_o"], fp["rho_g"], fp["rho_w"],
                          fp["Pbub"], API, fp["gasSG"], fp["temperature"],
                          fp["Pc"], fp["Tc"], fp["MgO"], fp["MgG"],
                          fp["Psep"], fp["Tsep"], gamma_o, rhoOSC);

    TransmissibilityCalculator transCalc(grid, rock);
    Simulator simulator(grid, rock, fluid, transCalc);

    const auto& sm = cfg["sim"];
    simulator.run(sm["totalTime"], sm["dt"]);

    return 0;
}