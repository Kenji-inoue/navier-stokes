#include <iostream>
#include "FileUtil.h"

FileUtil::FileUtil(const char* fileName) {
    fopen_s(&m_fprot, fileName, "w");
}

FileUtil::~FileUtil() {
    if (m_fprot) {
        fclose(m_fprot);
    }
}

void FileUtil::saveField(const Field2d& f, const char* title, int time, int interval) {
    if (time % interval != 0) {
        return;
    }

    fprintf(m_fprot, "iteration time = %4d : %s \n", time, title);
    for (const auto& row : f) {
        for (const auto& value : row) {
            fprintf(m_fprot, "%10.2e, ", value);
        }
        fprintf(m_fprot, "\n");
    }
}
