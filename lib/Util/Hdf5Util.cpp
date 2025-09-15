#include <stdexcept>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "Hdf5Util.h"

Hdf5Util::Hdf5Util(const std::string& fileName)
{
    m_fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (m_fileId < 0) {
        throw std::runtime_error("Failed to open HDF5 file: " + fileName);
    }

    if (H5Lexists(m_fileId, m_RESULT_DIR.c_str(), H5P_DEFAULT) <= 0) {
        hid_t group = H5Gcreate2(m_fileId, m_RESULT_DIR.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group);
    }
}

Hdf5Util::~Hdf5Util()
{
    H5Fclose(m_fileId);
}

void Hdf5Util::readIntConfig(const std::string& configName, int& value)
{
    const auto fullPath = m_CONFIG_DIR + configName;
    H5LTread_dataset_int(m_fileId, fullPath.c_str(), &value);
}

void Hdf5Util::readDoubleConfig(const std::string& configName, double& value)
{
    const auto fullPath = m_CONFIG_DIR + configName;
    H5LTread_dataset_double(m_fileId, fullPath.c_str(), &value);
}

void Hdf5Util::saveResult(const Field2d& data, const std::string& resultName, int timeIndex)
{
    const size_t ny = data.size();
    const size_t nx = data[0].size();

    std::vector<double> flat(ny * nx);
    for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i)
            flat[j * nx + i] = data[j][i];

    const std::string path = m_RESULT_DIR + resultName;

    hid_t dataset;
    if (H5Lexists(m_fileId, path.c_str(), H5P_DEFAULT) <= 0) {
        hsize_t dims[3] = {0, ny, nx};
        hsize_t maxdims[3] = {H5S_UNLIMITED, ny, nx};
        hsize_t chunkdims[3] = {1, ny, nx};
        hid_t dataspace = H5Screate_simple(3, dims, maxdims);
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, 3, chunkdims);

        dataset = H5Dcreate2(m_fileId, path.c_str(), H5T_NATIVE_DOUBLE,
                             dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

        H5Pclose(plist);
        H5Sclose(dataspace);
    } else {
        dataset = H5Dopen2(m_fileId, path.c_str(), H5P_DEFAULT);
    }

    hsize_t newdims[3] = {timeIndex + 1, ny, nx};
    H5Dset_extent(dataset, newdims);

    hsize_t start[3] = { (hsize_t)timeIndex, 0, 0 };
    hsize_t count[3] = { 1, ny, nx };
    hid_t dataspace = H5Dget_space(dataset);
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

    hid_t memspace = H5Screate_simple(3, count, nullptr);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, flat.data());

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}