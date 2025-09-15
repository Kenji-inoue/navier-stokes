#include <string>
#include <hdf5.h>
#include "typedef.h"

class Hdf5Util
{
public:
    Hdf5Util(const std::string& fileName);
    ~Hdf5Util();

    void readIntConfig(const std::string& configName, int& value);
    void readDoubleConfig(const std::string& configName, double& value);
    void saveResult(const Field2d& data, const std::string& resultName, int timeIndex);
private:
    hid_t m_fileId;
    const std::string m_CONFIG_DIR = "/config/";
    const std::string m_RESULT_DIR = "/result/";
};