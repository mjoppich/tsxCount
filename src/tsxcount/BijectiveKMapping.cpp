//
// Created by mjopp on 11.06.2016.
//

#include "BijectiveKMapping.h"


const char* BijectiveMappingNoInverseException::getExceptText() const
{

    std::stringstream oSS;

    oSS << "No matrix with det(A) = +/- 1 could be found for k= " << m_pMapping->getK() << std::endl;

    oSS << m_pMapping->printMatrix() << std::endl;

    return oSS.str().c_str();

}

const char* BijectiveMapping1DeterminantException::getExceptText() const
{

    std::stringstream oSS;

    oSS << "No matrix with det(A) = +/- 1 could be found for k= " << m_pMapping->getK() << std::endl;

    oSS << m_pMapping->printMatrix() << std::endl;

    return oSS.str().c_str();

}