//
// Created by mjopp on 15/03/2020.
//

#include <src/tsxcount/TSXHashMap.h>
#include "SBigInt.h"
#include "UBigInt.h"

SBIGINT::SBIGINT *SBIGINT::getEmptySBIGINT(uint16_t iKeyValBits) {

    std::div_t dfields = TSXHashMap::udiv( iKeyValBits, 8);

    uint8_t iFields = dfields.quot;

    if (dfields.rem > 0)
    {
        ++iFields;
    }

    return getEmptySBIGINT(iKeyValBits, iFields);

}

SBIGINT::SBIGINT *SBIGINT::fromClassToStruct(UBigInt *pElem, SBIGINT*ret) {


    for (uint8_t i = 0; i < pElem->m_iFields; ++i)
    {
        ret->pdata[i] = pElem->m_pArray[i];
    }

    return ret;

}

SBIGINT::SBIGINT *SBIGINT::fromClassToStruct(UBigInt *pElem) {

    SBIGINT* ret = getEmptySBIGINT(pElem->m_iBits, pElem->m_iFields);

    return fromClassToStruct(pElem, ret);

}

UBigInt SBIGINT::fromStructToClass(SBIGINT *pElem, MemoryPool<FIELDTYPE> *pPool) {
    UBigInt ret = UBigInt(pElem->iFields*pElem->iFieldSize, true, pPool);

    for (uint8_t i = 0; i < pElem->iFields; ++i)
    {
        ret.m_pArray[i] = pElem->pdata[i];
    }

    return ret;
}

SBIGINT::SBIGINT *SBIGINT::getEmptySBIGINT(uint16_t iKeyValBits, uint8_t iFields) {
    SBIGINT* ret = new SBIGINT();

    ret->iFieldSize = sizeof(FIELDTYPE) * 8;
    ret->iFields = iFields;
    ret->pdata = (FIELDTYPE*) malloc(sizeof(FIELDTYPE)* iFields);

    for (uint8_t i = 0; i < iFields; ++i)
    {
        ret->pdata[i] = 0;
    }

    return ret;
}

void SBIGINT::getFromMemory(SBIGINT *pRes, SBIGINT *pFields, uint8_t iOffset, uint16_t iKeyValBits,
                            FIELDTYPE *pData) {
    uint8_t totalFields = iKeyValBits / pRes->iFieldSize;
    FIELDTYPE iPartOne, iPartTwo;


    for (uint32_t i = 0;  i < totalFields; ++i)
    {

        pFields->pdata[i] = pData[i];
        iPartOne = pData[i] >> iOffset;
        iPartTwo = pData[i+1] << (pRes->iFieldSize-iOffset);

        pRes->pdata[i] = iPartOne | iPartTwo;

        iKeyValBits -= pRes->iFieldSize;

    }

    if (iKeyValBits > 0)
    {
        pFields->pdata[totalFields] = pData[totalFields];

        iPartOne = pData[totalFields];

        iPartOne = iPartOne << (pRes->iFieldSize-iOffset-iKeyValBits);
        iPartOne = iPartOne >> (pRes->iFieldSize-iOffset-iKeyValBits);
        iPartOne = iPartOne >> iOffset;
        pRes->pdata[totalFields] = iPartOne;
    }
}

void
SBIGINT::storeInMemory(FIELDTYPE *pSrc, FIELDTYPE *pDest, uint32_t iBitStart, uint16_t iBitCount, uint8_t iFieldSize) {

    //FIELDTYPE iThisVal, iOldVal, iPartRight, iPartLeft;
    FIELDTYPE values [4];
    // thisval = 0
    // oldval = 1
    // left = 2
    // right = 3

    uint8_t iOffset, iBitsRemaining, iFullFields, i;

    // how many bits can I copy into first, not fully occupied field?
    iOffset = iBitStart % iFieldSize;

    // TODO copy first field
    values[0] = pSrc[0] << iOffset;
    values[1] = pDest[0] << (iFieldSize - iOffset);
    values[1] = values[1] >> (iFieldSize - iOffset);

    pDest[0] = values[0] | values[1];

    // how many bits remain to be copied?
    iBitsRemaining = iBitCount - (iFieldSize - iOffset);
    iFullFields = iBitsRemaining / iFieldSize;

    //uint8_t iRemainFields = iBitsRemaining - iFullFields*iFieldSize;

    // copy full fields over

    for ( i = 0; i < iFullFields; ++i)
    {
        // TODO copy full fields
        values[3] = pSrc[i] >> (iFieldSize-iOffset);
        values[2] = pSrc[i+1] << iOffset;
        values[0] = values[2] | values[3];

        pDest[i+1] = values[0];

        iBitsRemaining -= iFieldSize;
    }


    if (iBitsRemaining > 0)
    {

        if (iBitsRemaining <= iOffset)
        {

            // we do not need a rest from the next field
            values[3] = pSrc[iFullFields] << (iFieldSize-iOffset);
            values[3] = values[3] << (iOffset-iBitsRemaining);
            values[3] = values[3] >> (iOffset-iBitsRemaining);
            values[3] = values[3] >> iFieldSize-iOffset;

            values[1] = pDest[iFullFields+1];
            values[1] = values[1] >> iBitsRemaining;
            values[1] = values[1] << iBitsRemaining;

            pDest[iFullFields+1] = values[1] | values[3];

        } else {

            // we need a rest from the next field
            values[3] = pSrc[iFullFields] >> (iFieldSize-iOffset);

            // values[0] is used to store iRemaing
            //iRemaining = iBitsRemaining - iOffset;
            values[0] = iBitsRemaining - iOffset;

            // need the rightmost iRemaining bits
            values[2] = pSrc[iFullFields+1] << (iFieldSize-values[0]); // clears all but iRemaining bits
            values[2] = values[2] >> (iFieldSize-values[0]);
            values[2] = values[2] << iOffset;

            values[0] = values[2] | values[3];

            values[1] = pDest[iFullFields+1];
            values[1] = values[1] >> iBitsRemaining;
            values[1] = values[1] << iBitsRemaining;

            pDest[iFullFields+1] = values[1] | values[0];

        }


        /*
        std::cout << "Having i: " << i << " and fields " << fields.quot << std::endl;
        std::cout << "Using the last bits: " << iBitsRemaining << std::endl;
        std::cout<< std::bitset<8>(pDest[fields.quot+1]) << " into field " << fields.quot+1 << " " << static_cast<void*>(pDest + fields.quot+1) << std::endl;
        */
    }
}

void SBIGINT::getFromMemory(SBIGINT *pRes, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE *pData) {
    uint8_t totalFields = iKeyValBits / pRes->iFieldSize;
    FIELDTYPE iPartOne, iPartTwo;


    for (uint32_t i = 0;  i < totalFields; ++i)
    {

        iPartOne = pData[i] >> iOffset;
        iPartTwo = pData[i+1] << (pRes->iFieldSize-iOffset);

        pRes->pdata[i] = iPartOne | iPartTwo;

        iKeyValBits -= pRes->iFieldSize;

    }

    if (iKeyValBits > 0)
    {
        iPartOne = pData[totalFields];
        iPartOne = iPartOne << (pRes->iFieldSize-iOffset-iKeyValBits);
        iPartOne = iPartOne >> (pRes->iFieldSize-iOffset-iKeyValBits);
        iPartOne = iPartOne >> iOffset;
        pRes->pdata[totalFields] = iPartOne;
    }
}

void SBIGINT::add_simpleSBIGINT(SBIGINT *pElem) {

    FIELDTYPE iCarriage = 1; // in order to add 1
    FIELDTYPE iMax = 0;

    uint32_t iFields = pElem->iFields;


    uint8_t i = 0;
    for (; i < iFields; ++i) {

        iMax = pElem->pdata[i];
        pElem->pdata[i] = pElem->pdata[i] + iCarriage;

        iCarriage = 0;

        if (pElem->pdata[i] < iMax)
            iCarriage = 1;

    }


    // technically we need to handle the case where more bits are used. But this is not needed here, because it is checked upfront

}

void SBIGINT::shiftLeft(SBIGINT *pElem, uint16_t iShift) {

    if (iShift >= pElem->iFieldSize)
    {

        uint16_t iFields = iShift / pElem->iFieldSize;

        uint32_t i;
        for (i = 0; i < pElem->iFields-iFields; ++i)
        {
            pElem->pdata[pElem->iFields-1 -i] = pElem->pdata[pElem->iFields-1-i-iFields];
        }

        for (i; i < pElem->iFields; ++i)
        {
            pElem->pdata[pElem->iFields-1 -i] = 0;
        }
    }

    iShift = iShift % pElem->iFieldSize;
    uint32_t iCounterShift = pElem->iFieldSize - iShift;

    FIELDTYPE iCurField, iNextFieldPart;

    for (uint32_t i = 0; i < pElem->iFields-1; ++i)
    {
        iCurField = (pElem->pdata[pElem->iFields-1-i] << iShift);
        iNextFieldPart = (pElem->pdata[pElem->iFields-1-i-1] >> iCounterShift);
        iCurField = iCurField | iNextFieldPart;

        pElem->pdata[pElem->iFields-1 -i] = iCurField;
    }

    pElem->pdata[0] = (pElem->pdata[0] << iShift);

}

void SBIGINT::shiftRight(SBIGINT *pElem, uint16_t iShift) {

    if (iShift >= pElem->iFieldSize)
    {
        uint64_t iFields = iShift / pElem->iFieldSize;

        uint32_t i;
        for (i = 0; i < pElem->iFields - iFields; ++i)
        {
            pElem->pdata[i] = pElem->pdata[i + iFields];
        }

        for (i; i < pElem->iFields; ++i)
        {
            pElem->pdata[i] = 0;
        }
    }

    iShift = iShift % pElem->iFieldSize;
    uint32_t iCounterShift = pElem->iFieldSize - iShift;

    for (uint32_t i = 0; i < pElem->iFields-1; ++i)
    {

        FIELDTYPE iNewRight = (pElem->pdata[i] >> iShift);
        FIELDTYPE iNewLeft = (pElem->pdata[i+1] << iCounterShift);

        FIELDTYPE iNewValue =  iNewLeft | iNewRight;
        pElem->pdata[i] = iNewValue;

    }

    pElem->pdata[pElem->iFields-1] = (pElem->pdata[pElem->iFields-1] >> iShift);
}

void SBIGINT::setOnes(SBIGINT *pRes, uint8_t iBits) {

    uint8_t iFullFields = iBits / pRes->iFieldSize;

    for (uint8_t i = 0; i < iFullFields; ++i)
    {
        pRes->pdata[i] = -1;
    }

    uint8_t remains = iBits-(iFullFields*pRes->iFieldSize);

    if (remains > 0)
    {
        for (uint8_t i = 0; i < remains; ++i)
        {
            pRes->pdata[pRes->iFields-1] = (pRes->pdata[pRes->iFields-1] & ~(1UL << i)) | (1 << i);
        }
    }

}

void SBIGINT::bitComplement(SBIGINT *pRes) {
    for (uint8_t i = 0; i < pRes->iFields; ++i)
    {
        pRes->pdata[i] = ~(pRes->pdata[i]);
    }
}

void SBIGINT::bitOr(SBIGINT *pRes, SBIGINT *pOther) {
    for (uint8_t i = 0; i < pRes->iFields; ++i)
    {
        pRes->pdata[i] |= (pOther->pdata[i]);
    }
}

bool SBIGINT::isZero(SBIGINT *pRes, uint8_t iBits) {
    div_t procs = div(iBits, pRes->iFieldSize);
    size_t checkIdx = 0;

    if (procs.quot > 0)
    {
        for (checkIdx = 0; checkIdx <= procs.quot-1; ++checkIdx)
        {
            if (pRes->pdata[checkIdx] != 0)
                return false;
        }
    }

    uint32_t iShift = pRes->iFieldSize - procs.rem;
    FIELDTYPE thisRemain = pRes->pdata[checkIdx] << iShift;
    bool result = (thisRemain == 0);

    return result;
}

void SBIGINT::clear(SBIGINT *pRes) {
    for (uint8_t i = 0; i < pRes->iFields; ++i)
    {
        pRes->pdata[i] = 0;
    }
}

void SBIGINT::bitAnd(SBIGINT *pRes, SBIGINT *pOther) {
    for (uint8_t i = 0; i < pRes->iFields; ++i)
    {
        pRes->pdata[i] &= (pOther->pdata[i]);
    }
}
