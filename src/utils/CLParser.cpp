/*
 * Splitter Application Package - some toolkit to work with gtf/gff files
 * Copyright (C) 2015  Markus Joppich
 *
 * The Splitter Application Package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Splitter Application Package is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cstring>

#include "CLParser.h"

CLParser::CLParser(std::string sArgs)
{

    const std::string sInput = std::string(sArgs);
    std::vector<std::string> vArgs = Utils::split(sInput, ' ');
    vArgs.insert(vArgs.begin(), std::string("."));

    char** ppArgs = new char*[vArgs.size()];

    for (int i = 0; i < vArgs.size(); ++i)
    {
        ppArgs[i] = new char[vArgs.at(i).length()];
        memcpy(ppArgs[i], vArgs.at(i).c_str(), vArgs.at(i).length());
    }

    this->initialize(vArgs.size(), ppArgs);

    for (int i = 0; i < vArgs.size(); ++i)
    {
        delete ppArgs[i];
    }

    delete ppArgs;
}

CLParser::CLParser(int argc, char** argv) {

    this->initialize(argc, argv);

    //this->printAllArguments();

}

CLParser::CLParser(const CLParser& orig) {
}

CLParser::~CLParser() {
}

bool CLParser::initialize(int argc, char** argv)
{
    m_pThisExecutable = NULL;
    m_pCLArguments = NULL;

    if (argc > 0)
        m_pThisExecutable = new std::string(argv[0]);

    m_pCLArguments = new std::map< std::string, std::string* >();

    std::map< std::string, std::string* >::iterator oIt;
    bool bAwaitsInput = false;

    for (int32_t i = 1; i < argc; ++i) {

        std::string* pArgument = new std::string(argv[i]);

        // does it start with - or -- ?
        if ((!bAwaitsInput) && ((pArgument->at(0) == '-') && ( pArgument->at(1) == '-'))) {



            std::string sArgument = pArgument->substr(2, pArgument->length());

            oIt = m_pCLArguments->insert(m_pCLArguments->end(), std::pair<std::string, std::string*>(sArgument, NULL));
            delete pArgument;

            bAwaitsInput = false;

            continue;
        }

        if ((!bAwaitsInput) && ((pArgument->at(0) == '-') && (pArgument->at(1) != '-'))) {

            std::string sArgument = pArgument->substr(1, pArgument->length());

            oIt = m_pCLArguments->insert(m_pCLArguments->end(), std::pair<std::string, std::string*>(sArgument, NULL));
            delete pArgument;

            bAwaitsInput = true;
            continue;
        }

        if (!bAwaitsInput) {
            this->handleParsingError(pArgument);

            return false;
        }

        // bAwaitsInput == true

        (*oIt).second = pArgument;
        bAwaitsInput = false;


    }

    return true;
}
