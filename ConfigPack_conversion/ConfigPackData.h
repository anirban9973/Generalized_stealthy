#ifndef CONFIGPACK_DATA_H
#define CONFIGPACK_DATA_H

#include <string>
#include <vector>
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>

class ConfigPackData {
private:
    ConfigurationPack cp;
    int total_configs;

public:
    // Constructor
    ConfigPackData(const std::string& configpack_name) 
        : cp(configpack_name) 
    {
        total_configs = cp.NumConfig();
    }
    
    // Remove 'const' from all these functions
    int NumConfigs() {
        return total_configs;
    }
    
    int GetDimension(int config_idx) {
        Configuration c = cp.GetConfig(config_idx);
        return c.GetDimension();
    }
    
    size_t GetNumParticles(int config_idx) {
        Configuration c = cp.GetConfig(config_idx);
        return c.NumParticle();
    }
    
    double GetVolume(int config_idx) {
        Configuration c = cp.GetConfig(config_idx);
        return c.PeriodicVolume();
    }
    
    std::vector<std::vector<double>> GetBasisVectors(int config_idx) {
        Configuration c = cp.GetConfig(config_idx);
        int dim = c.GetDimension();
        
        std::vector<std::vector<double>> basis(dim, std::vector<double>(dim));
        for (int i = 0; i < dim; ++i) {
            GeometryVector bv = c.GetBasisVector(i);
            for (int j = 0; j < dim; ++j) {
                basis[i][j] = bv.x[j];
            }
        }
        return basis;
    }
    
    std::vector<std::vector<double>> GetPositions(int config_idx) {
        Configuration c = cp.GetConfig(config_idx);
        int dim = c.GetDimension();
        size_t num_particles = c.NumParticle();
        
        std::vector<std::vector<double>> positions(num_particles, std::vector<double>(dim));
        for (size_t i = 0; i < num_particles; ++i) {
            GeometryVector pos = c.GetCartesianCoordinates(i);
            for (int j = 0; j < dim; ++j) {
                positions[i][j] = pos.x[j];
            }
        }
        return positions;
    }
};

#endif
