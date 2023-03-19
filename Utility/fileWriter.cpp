#include "fileWriter.h"

bool fileWriter::write(const std::string &text)
{
    std::ofstream file(_filepath, std::ios::app);

    if (!file.is_open())
    {
        std::cout << "Could not open file " << _filepath << std::endl;
        return false;
    }

    file << text;
    file.close();
    return true;
}


void fileWriter::addHeader(const std::string &text)
{
    //Create new file with _filename
    std::ofstream file(_filepath);

    //open file
    if (!file.is_open())
    {
        std::cout << "Could not open file " << _filepath << std::endl;
        return;
    }

    //write header
    file << text;
    file.close();
}


