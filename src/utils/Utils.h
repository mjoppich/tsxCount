//
// Created by mjoppich on 11/1/17.
//

#ifndef TSXCOUNT_SPLITUTILS_H_H
#define TSXCOUNT_SPLITUTILS_H_H

#include <string>
#include <vector>
#include <sstream>

#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <stdlib.h>
#include <functional>

class Utils
{
public:

    static std::vector<std::string> split(const std::string &sString, char cDelim) {

        std::vector<std::string> vret;

        std::stringstream sStringStream( sString );
        std::string sItem;
        while (std::getline(sStringStream, sItem, cDelim)) {
            vret.push_back(sItem);
        }
        return vret;
    }


    static bool file_exists(const std::string* name)
    {

        if (name == NULL)
            return false;

        struct stat buffer;
        return (stat (name->c_str(), &buffer) == 0);

    }
};

#endif //TSXCOUNT_SPLITUTILS_H_H
