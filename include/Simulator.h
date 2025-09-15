#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Grid.h"
#include "RockProperties.h"
#include "FluidProperties.h"
#include "TransmissibilityCalculator.h"
#include <vector>

class Simulator{
    public:
        Simulator(const Grid& grid, 
                  const RockProperties& rock, 
                  const FluidProperties& fluid, 
                  const TransmissibilityCalculator& transmissibilityCalculator);
        
        void run(double totalTime, double timeStep);

        //Accessors Pressure - Saturation
        double pressure(int cell) const;
        double waterSaturation(int cell) const;
        double gasSaturation(int cell) const;

    private:
            void initializeState();
            void solveTimestep(double dt);
            void computeResidual(std::vector<double>& residual, double dt);
            void applyBoundaryConditions(std::vector<double>& residual);

            void assembleJacobian(Eigen::SparseMatrix<double>& J, std::vector<double>& residual, double dt);

            double getgravity() const;

            const Grid& grid_;
            const RockProperties& rock_;
            const FluidProperties& fluid_;
            const TransmissibilityCalculator& transmissibilityCalculator_;

            //Primary variables
            std::vector<double> pressure_;
            std::vector<double> Sw_;
            std::vector<double> Sg_;
            std::vector<double> Rs_;

            //Previous State variables
            std::vector<double> pressure_prev_;
            std::vector<double> Sw_prev_;
            std::vector<double> Sg_prev_;
            std::vector<double> Rs_prev_;

            //Hydrostatic Pressure variables
            int nx_;
            int ny_;
            int nz_;
            double dx_;
            double dy_;
            double dz;

            //Physics units
            double gravity_ = 9.81; //m/s2

            //Flags saturation system
            std::vector<int> flags_; //0: Oil, 1:Oil + Rs, 2:Free gas present

            //Savign transmisibilities:
    std::vector<TransmissibilityCalculator::Transmissibility> transmissibilities_;

};