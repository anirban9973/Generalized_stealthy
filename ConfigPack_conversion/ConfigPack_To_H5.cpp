/*
 * Parse ConfigPack files into a single HDF5 file.
 * Output: ../h5/configs.h5
 *
 * File structure
 *   Global attributes : dim, N, n_configs, basis (dim×dim float64)
 *   Datasets          : config_0, config_1, …  each (N, dim) float64
 *
 * Positions are Cartesian coordinates at full double precision.
 * ConfigPack files matched: *GS_Success.ConfigPack in the current directory.
 */
#include "ConfigPackData.h"
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <H5Cpp.h>
#include <omp.h>

namespace fs = std::filesystem;

int main()
{
    const std::string suffix = "GS_Success.ConfigPack";

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
    // Create HDF5 file in the current working directory
    // ----------------------------------------------------------------
    fs::path h5path = "configs.h5";
    H5::H5File h5file(h5path.string(), H5F_ACC_TRUNC);

    // Write global attributes
    {
        H5::DataSpace scalar(H5S_SCALAR);

        H5::Attribute a_dim = h5file.createAttribute(
            "dim",       H5::PredType::NATIVE_INT, scalar);
        H5::Attribute a_N   = h5file.createAttribute(
            "N",         H5::PredType::NATIVE_INT, scalar);
        H5::Attribute a_nc  = h5file.createAttribute(
            "n_configs", H5::PredType::NATIVE_INT, scalar);

        a_dim.write(H5::PredType::NATIVE_INT, &dim_g);
        a_N  .write(H5::PredType::NATIVE_INT, &N_g);
        a_nc .write(H5::PredType::NATIVE_INT, &total_configs);

        // basis as (dim, dim) float64 attribute
        std::vector<double> basis_flat;
        for (auto& row : basis_g)
            for (double v : row) basis_flat.push_back(v);
        hsize_t bdims[2] = { (hsize_t)dim_g, (hsize_t)dim_g };
        H5::DataSpace bspace(2, bdims);
        H5::Attribute a_basis = h5file.createAttribute(
            "basis", H5::PredType::NATIVE_DOUBLE, bspace);
        a_basis.write(H5::PredType::NATIVE_DOUBLE, basis_flat.data());
    }

    // Pre-create all datasets (serial, before parallel section)
    {
        hsize_t pos_dims[2] = { (hsize_t)N_g, (hsize_t)dim_g };
        H5::DataSpace pos_space(2, pos_dims);
        for (int i = 0; i < total_configs; ++i) {
            h5file.createDataSet(
                "config_" + std::to_string(i),
                H5::PredType::NATIVE_DOUBLE,
                pos_space);
        }
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

            // Flatten row-major: pos[p][d]
            std::vector<double> flat;
            flat.reserve(N_g * dim_g);
            for (auto& p : pos)
                for (double v : p) flat.push_back(v);

            std::string dsname = "config_" + std::to_string(offsets[k] + i);
            #pragma omp critical
            {
                H5::DataSet ds = h5file.openDataSet(dsname);
                ds.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
            }
        }

        #pragma omp critical
        std::cout << "  Done: " << configpacks[k] << "\n";
    }

    h5file.close();
    std::cout << "\nTotal configs exported: " << total_configs << "\n";
    std::cout << "Output: " << h5path << "\n";
    return 0;
}
