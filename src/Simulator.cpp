#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include "Simulator.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Logger.h"

Simulator::Simulator(const Grid& grid, 
                     const RockProperties& rock, 
                     const FluidProperties& fluid, 
                     const TransmissibilityCalculator& transmissibilityCalculator):
                     grid_(grid), rock_(rock), fluid_(fluid), transmissibilityCalculator_(transmissibilityCalculator){

    //Transmisdibilities variables
    transmissibilities_ = transmissibilityCalculator.getTransmissibilities();

initializeState();
}

void Simulator::initializeState(){
    int n = grid_.numCells();
    pressure_.resize(n, 200.0e5); //Pressure in Pascal
    Sw_.resize(n, 0.2);
    Sg_.resize(n, 0.0);
    flags_.resize(grid_.numCells(), 0);

    //Calculate Hydrostatic pressure

        //Dimension Getters
        int nx_ = grid_.nx();
        int ny_ = grid_.ny();
        int nz_ = grid_.nz();

        double dx_ = grid_.dx();
        double dy_ = grid_.dy();
        double dz_ = grid_.dz();

        //Calculating Hydrostatic Pressure
        int nCells = grid_.numCells();
        std::vector<double> hspressure(nCells); //hydrostatic pressure

        for(int cellId = 0; cellId < nCells; ++cellId){
            double depth = grid_.cellDepth(cellId);                                                                                     //the depth is in meters
            double So_ = 1 - ( Sw_[cellId] + Sg_[cellId]);                                                                              //Oil Saturation i calculated as 1 - (Sw + Sg)
            double rho_avg = (Sw_[cellId] * fluid_.getRhoWdef()) + (Sg_[cellId] * fluid_.getRhoGdef()) + (So_ * fluid_.getRhoOdef());   //calculates average density (Sw*rho water) + (Sg*rho gas) + (So*rho oil)
            double g = getgravity();
            double depthref = grid_.getdepthref();
            double Pref = fluid_.getRhoWdef() * g * depthref;                                                                           //Pressure reference is currently calculated
            hspressure[cellId] = Pref + rho_avg * g * (depth - depthref);
        }
    
    //Reset Pressure to Hydrostatic Pressure
    pressure_ = hspressure;                                                                                                             //Pressure doubles vector is set to Hydrostatic pressure vector

    //Save Rs to implement Flag logic
    Rs_.resize(n);
    Rs_prev_.resize(n);
    for(int cell = 0; cell< n; ++cell){
            //Converting to Si from field units(scf/STB -> m3/m3)
            double convrs_si = 0.1781;
        Rs_[cell] = fluid_.getRs(pressure_[cell]) * convrs_si;
        Rs_prev_[cell] = Rs_[cell];
    }
}

void Simulator::run(double totalTime, double timeStep){
    Logger logger("simulation_log.jsonl");                                                                                               //Open the logger once
    
    //Compute cell centers
    std::vector<std::array<double, 3>> positions(grid_.numCells());
    for(int cell = 0; cell<grid_.numCells(); ++cell){
        positions[cell] = grid_.getCenter(cell);
    }

    int steps = static_cast<int>(totalTime / timeStep);
    for(int step = 0; step < steps; ++step){
        std::cout << "Step "<< step + 1 <<" / "<< steps << std::endl;
        solveTimestep(timeStep);                                                                                                         //Solve timestep

        //Compute residual norm for logging
        std::vector<double> residual(3 * grid_.numCells());
        computeResidual(residual, timeStep);

        double norm = 0.0;
        for(double r : residual){
            norm += r*r;
        }
        norm = std::sqrt(norm);

        //Preparing Scalars for logging
        std::map<std::string, double> scalars;
        scalars["residual norm"] = norm;

        //Log step
        logger.logStep(step, step * timeStep, scalars, pressure_, Sw_, Sg_, positions);
    }
}

void Simulator::solveTimestep(double dt){
    const int maxIters = 10;
    const double tol = 1e-5;
    int n = grid_.numCells();
    int size = 3 * n;

    //Save current state before computing residual
    pressure_prev_ = pressure_;
    Sw_prev_= Sw_;
    Sg_prev_ = Sg_;

    std::vector<double> residual(3 * n);    //water, oil, gas residuals per cell

    Rs_prev_ = Rs_;                         //Saving previous Rs

    for(int iter = 0; iter < maxIters; ++iter){
        computeResidual(residual, dt);

        //Compute residual norm;
        double norm = 0.0;
        for(double r : residual){
            norm += r*r;
        }

        norm = std::sqrt(norm);
        
        std::cout<< " Newton iter " << iter <<  ", residual norm = "<< norm << std::endl;

        if(norm < tol){
            std::cout<< " Converged.\n";
            return;
        }
        
        //Assemble Jacobian Matrix
        Eigen::SparseMatrix<double> J(size, size);
        assembleJacobian(J, residual, dt);

        //Build RHS
        Eigen::VectorXd rhs(size);
        for(int i = 0; i<size; ++i){
            rhs[i] = -residual[i];
        }

        //Solve J * dx = -residual

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(J);
        if(solver.info() != Eigen::Success){
            std::cerr << "Failed to factorize Jacobian. \n";

            break;
        }

        Eigen::VectorXd dx = solver.solve(rhs);
        if(solver.info() != Eigen::Success){
            std::cerr << "Failed to solve linear system.\n";

            break;
        }


        //Updating primary variables Using Newton Damping -- Damping dx
        
            //Newton Damping / Line search
            double alpha = 1.0;
            const double alpha_min = 1e-3;

            //Compute current residual before update
            double base_norm = norm;
            bool accepted = false;

            while (alpha >=alpha_min && !accepted){
                //Save previous state for backroll
                auto p_old = pressure_;
                auto Sw_old = Sw_;
                auto Sg_old = Sg_;

                //Appyly temporary update with damping
                for(int cell = 0; cell < n; ++cell){
                    pressure_[cell] = p_old[cell] + alpha * dx[3* cell + 0];
                    Sw_[cell] = std::clamp(Sw_old[cell] + alpha * dx[3 * cell + 1], 0.0, 1.0);
                    Sg_[cell] = std::clamp(Sg_old[cell] + alpha * dx[3 * cell + 2], 0.0, 1.0);
                }

                //Recompute residual with temporary state
                std::vector<double> r(size);
                computeResidual(r, dt);
                double new_norm = 0.0;
                for(double v : r){
                    new_norm += v*v;}
                new_norm = std::sqrt(new_norm);

                if(new_norm < base_norm){
                    //accept step
                    accepted = true;
                    std::cout << "Accepted step wit alpha=" << alpha
                              << " , new residual="<< new_norm << std::endl;
                } else {
                    // rollback to previous state and try smaller alpha
                    pressure_ = p_old;
                    Sw_ = Sw_old;
                    Sg_ = Sg_old;
                    alpha *= 0.5;
                }
            }
        if (!accepted){
            std::cerr<<" Line search failed, step not accepted, \n";
        } 
    }

    std::cerr<< " Newton convergance failed.\n";
}

void Simulator::computeResidual(std::vector<double>& residual, double dt){
    int n = grid_.numCells();
    residual.assign(3 * n, 0.0);                //reseting residual
    flags_.resize(grid_.numCells(), 0);

    //Computing Accumulation, Looping over cells
    for(int cell = 0; cell < n; ++cell){
        double phi = rock_.getPorosity(cell);
        double V = grid_.getVolume(cell);
        double p = pressure_[cell];
        double Sw = Sw_[cell];
        double Sg = Sg_[cell];
        double So = 1.0 - Sw -Sg;
        double Pb = fluid_.getBubblePointPressure();

        if(Sg > 1e-6){
            flags_[cell] = 2;
        } else if(p < Pb){
            flags_[cell] = 1;
        }else{
            flags_[cell] = 0;
        }

        //Previous state
        double p0 = pressure_prev_[cell];
        double Sw0 = Sw_prev_[cell];
        double Sg0 = Sg_prev_[cell];
        double So0 = 1.0 - Sw0 - Sg0;

        //Properties
            //Cnverting Bo and Bg from field units to Si (rb/scf -> m3/m3 & rb/STB -> m3/m3)
            double convbg_si = 5.615;
            double convbo_si = 1.2;
        double Bo = fluid_.getBo(p) * convbo_si;
        double Bg = fluid_.getBg(p) * convbg_si;
        /*double Rs = fluid_.getRs(p);*/

        //Implementing Flag logic on Rs
            double Rs;
            if(flags_[cell] == 2){
                //Converting Rs from scf/STB -> m3/m3 (field units to SI)
                double convrs_si = 0.1781;
                Rs = convrs_si * Rs_prev_[cell];            //Lock Rs if free gas present
            } else {
                double convrs_si = 0.1781;
                Rs = convrs_si * fluid_.getRs(p);           //otherwise update from pressure
                Rs_[cell] = Rs;                           //update current Rs
            }

        double phi_Bo = fluid_.getBo(p0) * convbo_si;
        double phi_Bg = fluid_.getBg(p0) * convbg_si;

        //Implementing Flag logic on Previous Rs
            double Rs0 = Rs_prev_[cell];

        //Accumulation
        double accW = phi * V * (Sw - Sw0);
       
        double accO = phi * V * ((So / Bo) - (So0 / phi_Bo)); //Oil component
        
        double accG = phi * V * ((Rs * So / Bo + Sg / Bg) - (Rs0 * So0 / phi_Bo + Sg0 / phi_Bg)); //Gas component

        residual[3 * cell + 0] += accW;
        residual[3 * cell + 1] += accO;
        residual[3 * cell + 2] += accG;
    }

    //Loop over transmissibilities & Flux computing
    for(const auto& T : transmissibilities_) {
        int i = T.cell1, j = T.cell2;
        double Tij = T.value;
            //Converting Harmonic Average Permeability from mD to m^2
            double Tij_SI = Tij * 9.869233e-16;

        //Geometry
        auto xi = grid_.getCenter(i);
        auto xj = grid_.getCenter(j);
        double dz = xj[2] - xi[2]; //for gravity

        //States
        double pi = pressure_[i], pj = pressure_[j];
        double Swi = Sw_[i], Swj = Sw_[j];
        double Sgi = Sg_[i], Sgj = Sg_[j];
        double Soi = 1 - Swi - Sgi, Soj = 1 - Swj - Sgj;

        //Mobilities
        double krwi = fluid_.getKrw(Swi), krwj = fluid_.getKrw(Swj);
        double kroi = fluid_.getKro(Soi), kroj = fluid_.getKro(Soj);
        double krgi = fluid_.getKrg(Sgi), krgj = fluid_.getKrg(Sgj);

        //Viscosities
        double muw = fluid_.getMuW();
        double muo = fluid_.getMuO((pi + pj) / 2.0);
        double mug = fluid_.getMuG((pi + pj) / 2.0);
            //Converting from cp to Pa*sec
            double muw_si = muw / 1000;
            double muo_si = muo;                //Already in SI - Pa*sec
            double mug_si = mug / 1000;

        double lamb_wi = krwi / muw_si, lamb_wj = krwj / muw_si;
        double lamb_oi = kroi / muo_si, lamb_oj = kroj / muo_si;
        double lamb_gi = krgi / mug_si, lamb_gj = krgj / mug_si;

        //Upwind pressure difference -- Central Averaging
        double dp = pj - pi;
        
        //Gravity (assume unit vector in -z and density constants for now)
        double g = getgravity();

        //Converting Densities from g/cm3 to kg/m3 (APART from rhoW, for now it returns default rhoW)

        double rho_w = fluid_.getRhoW(), rho_o = fluid_.getRhoO((pi + pj) / 2) * 1000, rho_g = fluid_.getRhoG((pi + pj) / 2) * 1000;
        double dpw = dp + rho_w * g * dz;
        double dpo = dp + rho_o * g * dz;
        double dpg = dp + rho_g * g * dz;

        double lamb_w_up = (dpw >= 0.0) ? lamb_wi : lamb_wj;
        double lamb_o_up = (dpo >= 0.0) ? lamb_oi : lamb_oj;
        double lamb_g_up = (dpg >= 0.0) ? lamb_gi : lamb_gj;

        double fluxW = Tij_SI * lamb_w_up * dpw;
        double fluxO = Tij_SI * lamb_o_up * dpo;
        double fluxG = Tij_SI * lamb_g_up * dpg;

        //Add Flux to residual (sign convention: positive outflow)
        residual[3 * i + 0] -= dt * fluxW;
        residual[3 * i + 1] -= dt * fluxO;
        residual[3 * i + 2] -= dt * fluxG;

        residual[3 * j + 0] += dt * fluxW;
        residual[3 * j + 1] += dt * fluxO;
        residual[3 * j + 2] += dt * fluxG;
    };

}

//Relative Epsilon adapted
void Simulator::assembleJacobian(Eigen::SparseMatrix<double>& J, std::vector<double>& residual, double dt){
    const double eps_rel = 1e-6; //relative pertubation
    const double eps_abs = 1e-8; //minimum pertubation

    int n = grid_.numCells();
    int size = 3 * n;

    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<double> baseResidual = residual;

    //Loop over each uknown
    for(int cell = 0; cell <n; ++cell){
        for(int var = 0; var < 3; ++var){
            //Backup
            double* target = nullptr;
            if(var == 0){
                target = &pressure_[cell];
            }
            else if(var == 1){
                target = &Sw_[cell];
            }
            else{
                target = &Sg_[cell];
            }
            double backup = *target;

            //Adaptive pertubation decision
            double eps = std::max(eps_abs, eps_rel * std::fabs(backup));
            *target += eps;

            //Compute pertubed residual
            std::vector<double> pertubedResidual(size);
            computeResidual(pertubedResidual, dt);

            //Finite difference derivative
            for (int eq = 0; eq < size; ++eq){
                double deriv = (pertubedResidual[eq] - baseResidual[eq]) / eps;
                triplets.emplace_back(eq, 3 * cell + var, deriv);
            }
            
            //Restore
            *target = backup;
        }
    }

    J.resize(size, size);
    J.setFromTriplets(triplets.begin(), triplets.end());
};

//Accessors for Pressure - Saturation
double Simulator::pressure(int cell) const{
    return pressure_[cell];
}

double Simulator::waterSaturation(int cell) const{
    return Sw_[cell];
}

double Simulator::gasSaturation(int cell) const{
    return Sg_[cell];
}

//Physics Getters
double Simulator::getgravity() const{
    return gravity_;
}