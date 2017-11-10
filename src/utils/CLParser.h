//
// Created by mjoppich on 11/1/17.
//

#ifndef TSXCOUNT_CLPARSER_H
#define TSXCOUNT_CLPARSER_H

#include <map>
#include <inttypes.h>
#include <iostream>
#include "Utils.h"


class CLParserException: public std::exception
{
public:
    CLParserException(std::string sText)
            : std::exception(), m_sText(sText)
    {
    }


    virtual const char* what() const throw()
    {
        return m_sText.c_str();
    }


protected:

    const std::string m_sText;

};

class CLParser {
public:
    CLParser(std::string sArgs);
    CLParser(int argc, char** argv);
    CLParser(const CLParser& orig);
    virtual ~CLParser();

    std::map< std::string, std::string* >* getArguments() {
        return m_pCLArguments;
    }

    std::string getArgumentByValue(std::string sArg) {

        std::map<std::string, std::string*>::iterator oIt = m_pCLArguments->find(sArg);

        if (oIt == m_pCLArguments->end())
        {
            throw new CLParserException("Does not have argument: " + sArg);
        }


        return std::string( *((*oIt).second) );

    }

    std::string* getArgument(std::string sArg) {

        std::map<std::string, std::string*>::iterator oIt = m_pCLArguments->find(sArg);

        if (oIt == m_pCLArguments->end())
        {
            return NULL;
        }


        return (*oIt).second;

    }

    std::string* getFileName(std::string sArg)
    {
        std::string* pFileName = this->getArgument(sArg);

        if (!Utils::file_exists(pFileName))
            return NULL;

        return pFileName;
    }

    bool isSet(std::string sArg) {
        std::map<std::string, std::string*>::iterator oIt = m_pCLArguments->find(sArg);

        if (oIt == m_pCLArguments->end())
            return false;

        return true;
    }

    void setArgument(std::string skey, std::string sval)
    {
        m_pCLArguments->insert(std::pair<std::string, std::string*>( skey, new std::string(sval) ));
    }

    void printAllArguments()
    {

        std::map< std::string, std::string*>::iterator oIt = m_pCLArguments->begin();

        for (; oIt != m_pCLArguments->end(); ++oIt) {

            std::cout << oIt->first << '\t' << ((oIt->second != NULL) ? (*oIt->second) : "") << std::endl;

        }

    }

    int getIntArgument(std::string arg)
    {
        std::string argval = this->getArgumentByValue(arg);

        int intVal = std::atoi(argval.c_str());

        return intVal;
    }

private:

    bool initialize(int argc, char** argv);

    void handleParsingError(std::string* pArgument) {
        std::cerr << "Error: invalid input arguments at " << std::endl;

        this->deleteArguments();

        delete pArgument;

        m_pCLArguments = NULL;
    }



    void deleteArguments() {

        std::map< std::string, std::string*>::iterator oIt = m_pCLArguments->begin();

        for (; oIt != m_pCLArguments->end(); ++oIt) {
            delete (*oIt).second;
        }

        m_pCLArguments->empty();

        delete m_pCLArguments;

    }

    bool m_bSuccess;

    std::map< std::string, std::string* >* m_pCLArguments;
    std::string* m_pThisExecutable;

};



#endif //TSXCOUNT_CLPARSER_H
