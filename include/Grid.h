#pragma once
#include <vector>
#include <array>

class Grid{

    public:
        Grid(int nx, int ny, int nz, double dx, double dy, double dz, double depthref);

        int numCells() const;
        std::array<int, 3> dimensions() const;

        double getVolume(int cellId) const;
        std::array<double, 3> getCenter(int cellId) const;
        std::vector<int> getNeighbors(int cellId) const;

        //Coordinates and Size Grid getters
        int nx() const;
        int ny() const;
        int nz() const;
        double dx() const;
        double dy() const;
        double dz() const;

        //Cell centers for gravity
        double cellDepth(int cellId) const;

        //Getter for depth reference
        double getdepthref() const;

        //Helper Get Vertical Neighbour Above
        int getVerticalNeighborAbove(int cellId) const;

    private:
        int nx_, ny_, nz_;
        double dx_, dy_, dz_;
        double depthref_; //in meters

        int to1D(int i, int j, int k) const;
        std::array<int, 3> to3D(int cellId) const;
};