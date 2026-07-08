/*
 * Parse ConfigPack files into a single HDF5 file.
 * Output: configs_N<total>_chi<chi>.h5   (e.g. configs_N2500_chi0.55.h5)
 *
 * File structure
 *   Global attributes : dim, N, n_configs, basis (dim×dim float64)
 *   Datasets          : config_0, config_1, …  each (N, dim) float64
 *
 * Positions are Cartesian coordinates at full double precision.
 *
 * Usage:
 *   configpack_to_h5           quenched:  gather *_Success.ConfigPack
 *                              -> configs_N<total>_chi<chi>.h5
 *   configpack_to_h5 thermal   thermal:   gather crystalmd_run*.ConfigPack (raw MD snapshots)
 *                              -> configs_thermal_N<total>_chi<chi>.h5
 * All matching packs in the current directory are accumulated across array tasks.
 */
#include "ConfigPackData.h"
#include <highfive/H5File.hpp>
#include <highfive/H5DataSpace.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <omp.h>

namespace fs = std::filesystem;

int main(int argc, char ** argv)
{
    // Mode: default = quenched (*_Success.ConfigPack); "thermal" = raw MD snapshots
    // (crystalmd_run*.ConfigPack). The mode also tags the output filename so the two
    // never collide in the same folder.
    bool thermal = (argc > 1 && std::string(argv[1]) == "thermal");
    const std::string tag = thermal ? "thermal_" : "";
    std::cout << "Mode: " << (thermal ? "thermal (crystalmd_run*.ConfigPack)"
                                       : "quenched (*_Success.ConfigPack)") << "\n\n";

    auto ends_with = [](const std::string& s, const std::string& suf) {
        return s.size() >= suf.size() &&
               s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
    };

    // ----------------------------------------------------------------
    // Collect matching ConfigPack names
    // ----------------------------------------------------------------
    std::vector<std::string> configpacks;
    for (const auto& entry : fs::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        std::string fn = entry.path().filename().string();

        bool match;
        if (thermal) {
            // raw MD snapshots: crystalmd_run<i>.ConfigPack; skip quenched/init packs
            match = fn.rfind("crystalmd_run", 0) == 0 &&
                    ends_with(fn, ".ConfigPack") &&
                    fn.find("_Success")    == std::string::npos &&
                    fn.find("_InitConfig") == std::string::npos;
        } else {
            match = ends_with(fn, "_Success.ConfigPack");
        }

        if (match)
            configpacks.push_back(
                fn.substr(0, fn.size() - std::string(".ConfigPack").size()));
    }
    std::sort(configpacks.begin(), configpacks.end());

    std::cout << "Found " << configpacks.size() << " matching ConfigPack(s):\n";
    for (const auto& cp : configpacks) std::cout << "  " << cp << "\n";
    std::cout << "\n";

    if (configpacks.empty()) {
        std::cerr << "No matching ConfigPack files found. Exiting.\n";
        return 1;
    }

    // ----------------------------------------------------------------
    // Pass 1 (sequential): count total configs and compute offsets
    // ----------------------------------------------------------------
    int n_packs = (int)configpacks.size();
    std::vector<int> offsets(n_packs + 1, 0);
    for (int k = 0; k < n_packs; ++k) {
        ConfigPackData tmp(configpacks[k]);
        offsets[k + 1] = offsets[k] + tmp.NumConfigs();
    }
    int total_configs = offsets[n_packs];
    std::cout << "Total configs to export: " << total_configs << "\n\n";

    // ----------------------------------------------------------------
    // Read global metadata from the first config of the first pack
    // ----------------------------------------------------------------
    int dim_g, N_g;
    std::vector<std::vector<double>> basis_g;
    {
        ConfigPackData first(configpacks[0]);
        dim_g   = first.GetDimension(0);
        N_g     = (int)first.GetNumParticles(0);
        basis_g = first.GetBasisVectors(0);
    }

    // ----------------------------------------------------------------
    // Read chi and N-per-side from params.dat (line 1: N per side, line 2: chi)
    // so each chi's configs.h5 is self-identifying.
    // ----------------------------------------------------------------
    double chi_val   = -1.0;
    int    N_perside = -1;
    {
        std::ifstream pf("params.dat");
        if (pf) {
            pf >> N_perside >> chi_val;
            std::cout << "params.dat: N_per_side = " << N_perside
                      << ", chi = " << chi_val << "\n\n";
        } else {
            std::cout << "WARNING: params.dat not found; "
                         "chi/N_per_side attributes set to -1.\n\n";
        }
    }

    // ----------------------------------------------------------------
    // Create HDF5 file. Filename encodes the total particle number (N_g = N*N)
    // and chi so files stay distinguishable by name, e.g. configs_N2500_chi0.55.h5
    // ----------------------------------------------------------------
    std::ostringstream h5name;
    h5name << "configs_" << tag << "N" << N_g << "_chi" << chi_val << ".h5";

    HighFive::File h5file(h5name.str(), HighFive::File::Truncate);

    h5file.createAttribute("dim",        dim_g);
    h5file.createAttribute("N",          N_g);          // total particles (from config)
    h5file.createAttribute("N_per_side", N_perside);    // from params.dat line 1
    h5file.createAttribute("chi",        chi_val);      // from params.dat line 2
    h5file.createAttribute("n_configs",  total_configs);
    h5file.createAttribute("basis",      basis_g);      // written as (dim, dim) float64

    // Pre-create all datasets (serial, before parallel section)
    {
        HighFive::DataSpace pos_space({(size_t)N_g, (size_t)dim_g});
        for (int i = 0; i < total_configs; ++i)
            h5file.createDataSet<double>("config_" + std::to_string(i), pos_space);
    }

    // ----------------------------------------------------------------
    // Pass 2 (parallel reads, serialized HDF5 writes)
    // ----------------------------------------------------------------
    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < n_packs; ++k) {

        ConfigPackData data(configpacks[k]);
        int tot = data.NumConfigs();

        #pragma omp critical
        std::cout << "Thread " << omp_get_thread_num()
                  << " processing: " << configpacks[k]
                  << " (" << tot << " configs, offset=" << offsets[k] << ")\n";

        for (int i = 0; i < tot; ++i) {
            auto pos = data.GetPositions(i);   // vector<vector<double>> (N, dim)
            std::string dsname = "config_" + std::to_string(offsets[k] + i);
            #pragma omp critical
            {
                h5file.getDataSet(dsname).write(pos);
            }
        }

        #pragma omp critical
        std::cout << "  Done: " << configpacks[k] << "\n";
    }

    std::cout << "\nTotal configs exported: " << total_configs << "\n";
    std::cout << "Output: " << h5name.str() << "\n";
    return 0;
}
