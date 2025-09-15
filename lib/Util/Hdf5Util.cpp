#include <stdexcept>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "Hdf5Util.h"

Hdf5Util::Hdf5Util(const std::string& fileName)
{
    m_fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (m_fileId < 0) {
        throw std::runtime_error("Failed to open HDF5 file: " + fileName);
    }
}

Hdf5Util::~Hdf5Util()
{
    H5Fclose(m_fileId);
}

void Hdf5Util::readIntConfig(const std::string& configName, int& value)
{
    const auto fullPath = m_configDir + configName;
    H5LTread_dataset_int(m_fileId, fullPath.c_str(), &value);
}

void Hdf5Util::readDoubleConfig(const std::string& configName, double& value)
{
    const auto fullPath = m_configDir + configName;
    H5LTread_dataset_double(m_fileId, fullPath.c_str(), &value);
}