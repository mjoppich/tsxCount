//
// Created by joppich on 6/28/16.
//

#include "UBigInt.h"

const char *UBigIntNotLargeEnough::getExceptText() const {

    std::stringstream oSS;
    oSS << m_iBitsQueried << " bits were queried, but object only has " << m_pBigInt->getBitCount();

    return oSS.str().c_str();

}

const char *UBigIntNotImplemented::getExceptText() const {

    std::stringstream oSS;
    oSS << "Function/Method " << m_sMethod << " is not yet implemented. Object at " << m_pBigInt;

    return oSS.str().c_str();

}

std::ostream& operator<< (std::ostream &oOS, const UBigInt &oInt)
{

    const std::string sPrint = oInt.to_string();

    oOS << sPrint;
    return oOS;
}