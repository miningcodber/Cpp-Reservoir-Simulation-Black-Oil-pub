#pragma once
#include "Grid.h"
#include "RockProperties.h"
#include <vector>
#include <tuple>

class TransmissibilityCalculator{
    public:
        TransmissibilityCalculator(const Grid& grid, const RockProperties& rock);

        struct Transmissibility{
            int cell1;
            int cell2;
            double value;
            double dz; //vertical difference between cell centers (cell2 - cell1)
        };

        const std::vector<Transmissibility>& getTransmissibilities() const;

    private:
        void computeTransmissibilities();
        const Grid& grid_;
        const RockProperties& rock_;
        std::vector<Transmissibility> transmissibilities_;
};