#include "RockProperties.h"

//Member initializer list for private members -- porosity_ and permeability_ vectors are constructed to default value
RockProperties::RockProperties(int numCells, double defaultPorosity, std::array<double, 3> defaultPermeability)
    : porosity_(numCells, defaultPorosity), permeability_(numCells, defaultPermeability){}

//For now homogeneous porosity
void RockProperties::setPorosity(int cellId, double value){
    porosity_[cellId] = value;
}

//For now homogeneous permeability
void RockProperties::setPermeability(int cellId, std::array<double, 3> value){
    permeability_[cellId] = value;
}

//Method to change porosity on individual cell level
double RockProperties::getPorosity(int cellId) const{
    return porosity_[cellId];
}

//Method to change permeability on individual cell level
std::array<double, 3> RockProperties::getPermeability(int cellId) const{
    return permeability_[cellId];
}