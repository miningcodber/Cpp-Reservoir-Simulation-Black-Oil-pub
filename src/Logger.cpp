#include"Logger.h"
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Logger::Logger(const std::string& filename){
    file_.open(filename, std::ios::out);
}

Logger::~Logger(){
    if(file_.is_open()){
        file_.close();
    }
}

void Logger::logStep(int step, double time, const std::map<std::string, double>& scalars, const std::vector<double>& pressures, const std::vector<double>& Sw, const std::vector<double>& Sg, const std::vector<std::array<double, 3>>& positions){
    json entry;
    entry["step"] = step;
    entry["time"] = time;

    for(auto& kv : scalars){
        entry["scalars"][kv.first] = kv.second;
    }

    entry["fields"]["pressure"] = pressures;
    entry["fields"]["Sw"] = Sw;
    entry["fields"]["Sg"] = Sg;

    // Compute So
    std::vector<double> So(pressures.size());
    for(size_t i=0; i<Sw.size(); ++i)
        So[i] = 1.0 - Sw[i] - Sg[i];
    entry["fields"]["So"] = So;

    entry["positions"] = positions;

    file_ << entry.dump() << std::endl;
    file_.flush();
}