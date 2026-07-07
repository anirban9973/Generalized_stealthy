/*
 * Parse ConfigPack files into a single HDF5 file.
 * Output: configs.h5
 *
 * File structure
 *   Global attributes : dim, N, n_configs, basis (dim×dim float64)
 *   Datasets          : config_0, config_1, …  each (N, dim) float64
 *
 * Positions are Cartesian coordinates at full double precision.
 * ConfigPack files matched: *_Success.ConfigPack in the current directory
 * (e.g. crystal_run1_Success.ConfigPack, ..., accumulated across all array tasks).
 */
#include "ConfigPackData.h"
#include <highfive/H5File.hpp>
#include <highfive/H5DataSpace.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <omp.h>

namespace fs = std::filesystem;

int main()
{
    const std::string suffix = "_Success.ConfigPack";

    // ----------------------------------------------------------------
    // Collect matching ConfigPack names
    // ----------------------------------------------------------------
    std::vector<std::string> configpacks;
    for (const auto& entry : fs::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        std::string fn = entry.path().filename().string();
        if (fn.size() >= suffix.size() &&
            fn.compare(fn.size() - suffix.size(), suffix.size(), suffix) == 0)
        {
            configpacks.push_back(
                fn.substr(0, fn.size() - std::string(".ConfigPack").size()));
        }
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
    // Create HDF5 file
    // ----------------------------------------------------------------
    HighFive::File h5file("configs.h5", HighFive::File::Truncate);

    h5file.createAttribute("dim",       dim_g);
    h5file.createAttribute("N",         N_g);
    h5file.createAttribute("n_configs", total_configs);
    h5file.createAttribute("basis",     basis_g);   // written as (dim, dim) float64

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
    std::cout << "Output: configs.h5\n";
    return 0;
}
