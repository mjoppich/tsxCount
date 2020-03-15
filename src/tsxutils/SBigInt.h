//
// Created by mjopp on 26/02/2020.
//

#ifndef TSXCOUNT_SBIGINT_H
#define TSXCOUNT_SBIGINT_H

#include <src/tsxcount/commons.h>

class UBigInt;

namespace SBIGINT {


    struct SBIGINT{

        FIELDTYPE* pdata;
        uint8_t iFields;
        uint8_t iFieldSize;

    } ;



    void add_simpleSBIGINT(SBIGINT* pElem);


    void shiftLeft( SBIGINT* pElem, uint16_t iShift);

    void shiftRight( SBIGINT* pElem, uint16_t iShift);

    void getFromMemory(SBIGINT* pRes, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE* pData );

    void storeInMemory(FIELDTYPE* pSrc, FIELDTYPE* pDest,uint32_t iBitStart, uint16_t iBitCount, uint8_t iFieldSize);

    void getFromMemory(SBIGINT* pRes, SBIGINT* pFields, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE* pData );

    SBIGINT* getEmptySBIGINT(uint16_t iKeyValBits, uint8_t iFields);


    SBIGINT* getEmptySBIGINT(uint16_t iKeyValBits);


    SBIGINT* fromClassToStruct(UBigInt* pElem, SBIGINT* ret);

    SBIGINT* fromClassToStruct(UBigInt* pElem);


    UBigInt fromStructToClass(SBIGINT* pElem, MemoryPool<FIELDTYPE>* pPool);

    void setOnes(SBIGINT* pRes, uint8_t iBits);

    void bitComplement(SBIGINT* pRes);

    void bitOr(SBIGINT* pRes, SBIGINT* pOther);

    void bitAnd(SBIGINT* pRes, SBIGINT* pOther);

    void clear(SBIGINT* pRes);

    bool isZero(SBIGINT* pRes, uint8_t iBits);

}



#endif //TSXCOUNT_SBIGINT_H
