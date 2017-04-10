//
// Created by mjopp on 11.06.2016.
//

#ifndef TSXCOUNT_IBIJECTIVEFUNCTION_H
#define TSXCOUNT_IBIJECTIVEFUNCTION_H

#include <stdint.h>

#include "TSXTypes.h"



/**
 * \class IBijectiveFunction
 *
 *
 * \brief interface for Bijective Functions.
 *
 * This class defines two operations on a bijective function. A bijective function must implement apply and inverse apply
 *
 *
 * \author Markus Joppich
 *
 */
class IBijectiveFunction
{

public:

    virtual UBigInt apply(UBigInt& iKmer) = 0;
    virtual UBigInt inv_apply(UBigInt& iKey) = 0;

};


#endif //TSXCOUNT_IBIJECTIVEFUNCTION_H
