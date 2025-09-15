#include "TransmissibilityCalculator.h"
#include <cmath>
#include <array>


//Member initializer list for private members - Constructor
TransmissibilityCalculator::TransmissibilityCalculator(const Grid& grid, const RockProperties& rock)
    :grid_(grid), rock_(rock){
        computeTransmissibilities();
};

//Method that uses getters from Grid
void TransmissibilityCalculator::computeTransmissibilities(){
    int nx = grid_.nx();
    int ny = grid_.ny();
    int nz = grid_.nz();
    double dx = grid_.dx();
    double dy = grid_.dy();
    double dz = grid_.dz();

    //Lambda function for 3D to 1D indexing - same as in to1D from Grid
    auto index = [nx, ny](int i, int j, int k){
        return i + j * nx + k * nx * ny;
    };

    for(int k = 0; k<nz; ++k){
        for(int j=0; j<ny; ++j){
            for(int i=0; i<nx; ++i){
                int cell = index(i, j, k);
                std::array<double, 3> perm = rock_.getPermeability(cell);

                //X-Direction neighbor
                if(i < nx-1){
                    int neighbor = index(i + 1, j, k);
                    std::array<double, 3> permN = rock_.getPermeability(neighbor);
                    //Calculate Harmonic Avarage Permeability
                        //Avoiding division by zero
                        double k1 = perm[0];
                        double k2 = permN[0];

                        //Setting Harmonic Average Permeability calculation by ternary operation - Plain harmonic average for equal-distant cells (no Δx weighting)
                        double k_harmonic = (k1 + k2 > 0.0) ? (2.0 * k1 * k2 / (k1+k2)) : 0.0;

                    //Calclate Geometric Transmissibility
                    double T = k_harmonic * dy * dz / dx;
                    //Save transmisibility
                    transmissibilities_.push_back({cell, neighbor, T, 0.0});
                }

                //Y-Direction neighbor
                if (j < ny - 1) {
                    int neighbor = index(i, j + 1, k);
                    std::array<double, 3> permN = rock_.getPermeability(neighbor);
                    //Calculate Harmonic Avarage Permeability
                        //Avoiding division by zero
                        double k1 = perm[1];
                        double k2 = permN[1];

                        //Setting Harmonic Average Permeability calculation by ternary operation - Plain harmonic average for equal-distant cells (no Δx weighting)
                        double k_harmonic = (k1 + k2 > 0.0) ? (2.0 * k1 * k2 / (k1+k2)) : 0.0;

                    //Calclate Geometric Transmissibility
                    double T = k_harmonic * dx * dz / dy;
                    //Save transmisibility
                    transmissibilities_.push_back({cell, neighbor, T, 0.0});
                }

                //Z-Direction neighbor
                if (k < nz - 1) {
                    int neighbor = index(i, j, k + 1);
                    std::array<double, 3> permN = rock_.getPermeability(neighbor);
                    //Calculate Harmonic Avarage Permeability
                        //Avoiding division by zero
                        double k1 = perm[2];
                        double k2 = permN[2];

                        //Setting Harmonic Average Permeability calculation by ternary operation - Plain harmonic average for equal-distant cells (no Δx weighting)
                        double k_harmonic = (k1 + k2 > 0.0) ? (2.0 * k1 * k2 / (k1+k2)) : 0.0;

                    //Calclate Geometric Transmissibility
                    double T = k_harmonic * dx * dy / dz;

                    //Cell depth difference for gravity effect
                    double z1 = grid_.cellDepth(cell);
                    double z2 = grid_.cellDepth(neighbor);
                    double dz_gravity = z2 - z1;
                    //Save transmisibility
                    transmissibilities_.push_back({cell, neighbor, T, dz_gravity});
                }
            }
        }

    //Computing transmissibilities for positive side ONLY as they are symmetric.
    }
}
const std::vector<TransmissibilityCalculator::Transmissibility>& 
TransmissibilityCalculator::getTransmissibilities() const {
    return transmissibilities_;
}