#pragma once
#include <vector>
#include <array>

class RockProperties{
    public:
        RockProperties(int numCells, double defaultPorosity, std::array<double, 3> defaultPermeability);
        void setPorosity(int cellId, double value);
        void setPermeability(int cellId, std::array<double, 3> value);

        double getPorosity(int cellId) const;
        std::array<double, 3> getPermeability(int cellId) const;

    private:
        std::vector<double> porosity_;                      //size = numCells
        std::vector<std::array<double, 3>> permeability_;    //kx, ky, kz per cell
};