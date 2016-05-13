//
// Created by Peter Rodenkirch on 13.05.16.
//

#include "ParameterParser.h"

#include <sstream>
#include <fstream>


ParameterParser::ParameterParser()
{

}


std::map<std::string, double> ParameterParser::parseFile(std::string path)
{
    std::map<std::string, double> parameterMap;

    // Loads the file content in a string
    std::ifstream t(path);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    std::istringstream fileString(str);

    std::string line;
    while(std::getline(fileString, line))
    {
        std::istringstream is_line(line);
        std::string key;

        if(std::getline(is_line, key, '='))
        {
            std::string value;
            if(std::getline(is_line, value))
            {
                parameterMap[key] = ::atof(value.c_str());
            }
        }
    }

    return parameterMap;
}