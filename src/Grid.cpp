#include "Grid.h"

//Member initializer list to directly initialize private members
Grid::Grid(int nx, int ny, int nz, double dx, double dy, double dz, double depthref)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), depthref_(depthref) {}

    //Nx - Ny - Nz are the number of cells in the x - y - z direction.
    int Grid::numCells() const{
        return nx_ *ny_ *nz_;
    }

    std::array<int, 3> Grid::dimensions() const{
        return{nx_, ny_, nz_};
    }

    //Dx - Dy - Dz is the size of each cell respectively
        //For now getVolume doesnt utilize cellId, as all cells are considered to be equal size
    double Grid::getVolume(int cellId) const{
        return dx_ *dy_ *dz_;
    }

    //Returns coordinates of cell centers in an array of integers
    std::array<double, 3> Grid::getCenter(int cellId) const {
        auto [i, j, k] = to3D(cellId);
        return{(i+0.5)*dx_, (j+0.5)*dy_, (k+0.5)*dz_};
    }
    
    //Returns vector called neighbors of cellId - integers
    std::vector<int> Grid::getNeighbors(int cellId) const{
        std::vector<int> neighbors;
        auto [i, j, k] = to3D(cellId);
        
        //Arrays of possible neighbors to cell
        const std::array<std::array<int, 3>, 6> directions = {{
            {-1, 0, 0}, {1, 0, 0},
            {0, -1, 0}, {0, 1, 0},
            {0, 0, -1}, {0, 0, 1}
        }};

        //Adds the i - j - k coordinates toe each array to check the existance of neighbouring cells
        for(const auto& d: directions){
            int ni = i+d[0];
            int nj = j+d[1];
            int nk = k+d[2];
            
            //Neighbor condition check
            if (ni >= 0 && ni < nx_ && nj >= 0 && nj < ny_ && nk >= 0 && nk < nz_) {
            neighbors.push_back(to1D(ni, nj, nk));
            }
        }
        return neighbors;
    }

    //Converts cartesian coordinates to cellId
    int Grid::to1D(int i, int j, int k) const{
        return i+nx_*(j+ny_*k);
    }

    //Converts cellId to i - j - k cartesian coordinates
    std::array<int, 3> Grid::to3D(int cellId) const{
        int k = cellId / (nx_* ny_);
        int j = (cellId %(nx_*ny_)) / nx_;
        int i = cellId % nx_;
        return{i, j, k};
    }

    //Coordinates and Size Grid Getters
    int Grid::nx() const { return nx_; }
    int Grid::ny() const { return ny_; }
    int Grid::nz() const { return nz_; }
    double Grid::dx() const { return dx_; }
    double Grid::dy() const { return dy_; }
    double Grid::dz() const { return dz_; }


    //Cell depth for gravity
    double Grid::cellDepth(int cellId) const{
        return getCenter(cellId)[2]; // Z coordinate of the cell center
    }

    //Getter for depth reference
    double Grid::getdepthref() const{
        return depthref_; //in meters
    }
    
    //Helper to get Vertical Neighbors
    int Grid::getVerticalNeighborAbove(int cellId) const{
        auto [i, j, k] = to3D(cellId);
        if(k+1 < nz_){
            return to1D(i, j, k+1);
        } else {
            return -1; //or throw
        }
    }