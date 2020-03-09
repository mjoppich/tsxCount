//
// Created by mjopp on 26/02/2020.
//

#ifndef TSXCOUNT_SBIGINT_H
#define TSXCOUNT_SBIGINT_H

#include <inttypes.h>
#include <src/tsxcount/commons.h>

namespace SBIGINT {


    struct SBIGINT{

        FIELDTYPE* pdata;
        uint8_t iFields;
        uint8_t iFieldSize;

    } ;



    void add_simpleSBIGINT(SBIGINT* pElem) {

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


    void shiftLeft( SBIGINT* pElem, uint16_t iShift)
    {

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

    void shiftRight( SBIGINT* pElem, uint16_t iShift)
    {

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

    void getFromMemory(SBIGINT* pRes, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE* pData )
    {
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

    void storeInMemory(FIELDTYPE* pSrc, FIELDTYPE* pDest,uint32_t iBitStart, uint16_t iBitCount, uint8_t iFieldSize) {

        // how many bits can I copy into first, not fully occupied field?
        uint8_t iOffset = iBitStart % iFieldSize;

        // TODO copy first field
        FIELDTYPE iThisVal = pSrc[0] << iOffset;
        FIELDTYPE iOldVal = pDest[0] << (iFieldSize - iOffset);
        iOldVal = iOldVal >> (iFieldSize - iOffset);

        pDest[0] = iThisVal | iOldVal;

        // how many bits remain to be copied?
        uint32_t iBitsRemaining = iBitCount - (iFieldSize - iOffset);
        // how many full fields are this?
        div_t fields = div(iBitsRemaining, iFieldSize);

        // copy full fields over
        FIELDTYPE iPartRight;
        FIELDTYPE iPartLeft;

        uint32_t i = 0;
        for ( i = 0; i < fields.quot; ++i)
        {
            // TODO copy full fields
            iPartRight = pSrc[i] >> (iFieldSize-iOffset);
            iPartLeft = pSrc[i+1] << iOffset;
            iThisVal = iPartLeft | iPartRight;

            pDest[i+1] = iThisVal;

            iBitsRemaining -= iFieldSize;
        }

        if (iBitsRemaining > 0)
        {
            if (iBitsRemaining <= iOffset)
            {

                // we do not need a rest from the next field
                iPartRight = pSrc[fields.quot] << (iFieldSize-iOffset);
                iPartRight = iPartRight << (iOffset-iBitsRemaining);
                iPartRight = iPartRight >> (iOffset-iBitsRemaining);
                iPartRight = iPartRight >> iFieldSize-iOffset;

                iOldVal = pDest[fields.quot+1];
                iOldVal = iOldVal >> iBitsRemaining;
                iOldVal = iOldVal << iBitsRemaining;

                pDest[fields.quot+1] = iOldVal | iPartRight;

            } else {

                // we need a rest from the next field
                iPartRight = pSrc[fields.quot] >> (iFieldSize-iOffset);

                uint32_t iRemaining = iBitsRemaining - iOffset;

                // need the rightmost iRemaining bits
                iPartLeft = pSrc[fields.quot+1] << (iFieldSize-iRemaining); // clears all but iRemaining bits
                iPartLeft = iPartLeft >> (iFieldSize-iRemaining);
                iPartLeft = iPartLeft << iOffset;

                iThisVal = iPartLeft | iPartRight;

                iOldVal = pDest[fields.quot+1];
                iOldVal = iOldVal >> iBitsRemaining;
                iOldVal = iOldVal << iBitsRemaining;

                pDest[fields.quot+1] = iOldVal | iThisVal;

            }
            /*
            std::cout << "Having i: " << i << " and fields " << fields.quot << std::endl;
            std::cout << "Using the last bits: " << iBitsRemaining << std::endl;
            std::cout<< std::bitset<8>(pDest[fields.quot+1]) << " into field " << fields.quot+1 << " " << static_cast<void*>(pDest + fields.quot+1) << std::endl;
            */
        }
    }

    void getFromMemory(SBIGINT* pRes, SBIGINT* pFields, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE* pData )
    {
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

    SBIGINT* getEmptySBIGINT(uint16_t iKeyValBits, uint8_t iFields)
    {
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


    SBIGINT* getEmptySBIGINT(uint16_t iKeyValBits)
    {

        std::div_t dfields = TSXHashMap::udiv( iKeyValBits, 8);

        uint8_t iFields = dfields.quot;

        if (dfields.rem > 0)
        {
            ++iFields;
        }

        return getEmptySBIGINT(iKeyValBits, iFields);

    }

    SBIGINT* fromClassToStruct(UBigInt* pElem)
    {

        SBIGINT* ret = getEmptySBIGINT(pElem->m_iBits, pElem->m_iFields);

        for (uint8_t i = 0; i < pElem->m_iFields; ++i)
        {
            ret->pdata[i] = pElem->m_pArray[i];
        }

        return ret;

    }

    UBigInt fromStructToClass(SBIGINT* pElem, MemoryPool<FIELDTYPE>* pPool)
    {
        UBigInt ret = UBigInt(pElem->iFields*pElem->iFieldSize, true, pPool);

        for (uint8_t i = 0; i < pElem->iFields; ++i)
        {
            ret.m_pArray[i] = pElem->pdata[i];
        }

        return ret;
    }

    void setOnes(SBIGINT* pRes, uint8_t iBits)
    {

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

    void bitComplement(SBIGINT* pRes)
    {
        for (uint8_t i = 0; i < pRes->iFields; ++i)
        {
            pRes->pdata[i] = ~(pRes->pdata[i]);
        }
    }

    void bitOr(SBIGINT* pRes, SBIGINT* pOther)
    {
        for (uint8_t i = 0; i < pRes->iFields; ++i)
        {
            pRes->pdata[i] |= (pOther->pdata[i]);
        }
    }

    void bitAnd(SBIGINT* pRes, SBIGINT* pOther)
    {
        for (uint8_t i = 0; i < pRes->iFields; ++i)
        {
            pRes->pdata[i] &= (pOther->pdata[i]);
        }
    }

    void clear(SBIGINT* pRes)
    {
        for (uint8_t i = 0; i < pRes->iFields; ++i)
        {
            pRes->pdata[i] = 0;
        }
    }

    bool isZero(SBIGINT* pRes, uint8_t iBits)
    {
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

}



#endif //TSXCOUNT_SBIGINT_H