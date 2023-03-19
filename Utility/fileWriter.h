#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <cstdlib>
#include <system_error>

class fileWriter
{
public:
    fileWriter(){};
    fileWriter(const std::string &file) : _filepath(file){};
    ~fileWriter(){};

    bool write(const std::string &text);
    void addHeader(const std::string &header);

private:
    std::string _filepath;
    bool isDirectoryDeleted = false;
};

#endif