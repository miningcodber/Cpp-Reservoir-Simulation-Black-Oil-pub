#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <map>

class Logger{

public:
    Logger(const std::string& filename);
    ~Logger();

    void logStep(int step, double time, 
                 const std::map<std::string, double>& scalars, 
                 const std::vector<double>& pressures, 
                 const std::vector<double>& Sw, 
                 const std::vector<double>& Sg,
                 const std::vector<std::array<double, 3>>& positions);

private:
    std::ofstream file_;
};