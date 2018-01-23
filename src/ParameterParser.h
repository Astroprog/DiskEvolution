//
// Created by Peter Rodenkirch on 13.05.16.
//

#ifndef DISKEVOLUTION_PARAMETERPARSER_H
#define DISKEVOLUTION_PARAMETERPARSER_H

#include <string>
#include <map>

class ParameterParser {

public:
    ParameterParser();
    static std::map<std::string, double> parseFile(std::string path);

};


#endif //DISKEVOLUTION_PARAMETERPARSER_H
