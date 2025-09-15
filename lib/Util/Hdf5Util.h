#include <string>
#include <hdf5.h>

class Hdf5Util
{
public:
    Hdf5Util(const std::string& fileName);
    ~Hdf5Util();

    void readIntConfig(const std::string& configName, int& value);
    void readDoubleConfig(const std::string& configName, double& value);
private:
    hid_t m_fileId;
    const std::string m_configDir = "/config/";
};