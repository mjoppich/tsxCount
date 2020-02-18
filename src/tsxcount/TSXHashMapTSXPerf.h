//
// Created by mjopp on 03/02/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPTSXPERF_H
#define TSXCOUNT_TSXHASHMAPTSXPERF_H


#include "TSXHashMap.h"
#include <pthread.h>
#include <stack>


#include <immintrin.h>
#include <unistd.h>
#include <rtmintrin.h>
#include <cmath>



TSX::tsx_kmer_t fromSequence(std::string& seq, MemoryPool<FIELDTYPE>* pPool)
{

    UBigInt oRet(seq.length()*2, false, pPool);
    bool hadN = false;

    for (size_t i = 0; i < seq.length(); ++i)
    {

        //size_t iBitPos = 2*(seq.length() -1-i);
        size_t iBitPos = 2*i;

        switch (seq.at(i))
        {
            case 'A': // 0 00

                oRet.setBit(iBitPos, 0);
                oRet.setBit(iBitPos+1, 0);

                break;

            case 'C': // 1 01

                oRet.setBit(iBitPos, 1);
                oRet.setBit(iBitPos+1, 0);

                break;

            case 'G': // 2 10

                oRet.setBit(iBitPos, 0);
                oRet.setBit(iBitPos+1, 1);

                break;

            case 'T': // 3 11
                oRet.setBit(iBitPos, 1);
                oRet.setBit(iBitPos+1, 1);
                break;

            default: // random e.g. N

                uint8_t iBit1 = rand() % 2;
                uint8_t iBit2 = rand() % 2;

                oRet.setBit(iBitPos, iBit1);
                oRet.setBit(iBitPos+1, iBit2);

                hadN = true;

                break;
        }
    }

    return oRet;

}



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

void bitComplement(SBIGINT* pRes)
{
    for (uint8_t i = 0; i < pRes->iFields; ++i)
    {
        pRes->pdata[i] = ~(pRes->pdata[i]);
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


FIELDTYPE volatile PREFETCH;


class TSXHashMapTSXPerf : public TSXHashMap {

public:

    TSXHashMapTSXPerf(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
    : TSXHashMap(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    virtual ~TSXHashMapTSXPerf()
    {
        free(m_pLocked);
    }
    size_t iAddCount = 0;
    size_t iAddKmerCount = 0;
    size_t iAborts = 0;
    size_t iTotalAborts = 0;


    /**
     *
     * @param kmer
     * @param verbose
     * @return
     */
    virtual bool addKmer_tsx(TSX::tsx_kmer_t& kmer, bool verbose=false)
    {
        bool bInserted = false;

        this->iAddKmerCount+=1;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        uint8_t iThreadID = omp_get_thread_num();
        uint8_t iRetries = 0;
        UBigInt uBigOne = UBigInt(1, m_pPool);

        UBigInt* p_mask_func_reprobe = &m_mask_func_reprobe;

        while ( iReprobes < m_iMaxReprobes ) {


            /*
             *
             * Let's assume the simple case: the field is empty
             *
             */

            uint64_t iPos = this->getPosition(basekey, iReprobes);

            TSX::tsx_keyval_t key_reprobe_shift = this->makeKey(basekey, iReprobes);
            key_reprobe_shift.resize(m_iKeyValBits);
            key_reprobe_shift = (key_reprobe_shift << m_iStorageBits);


            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t* pSavedKey = &savedkey;
            TSX::tsx_keyval_t* pKeyReprobeShiftUBIGINT = &key_reprobe_shift;

            uint64_t iBitsToPos = (iPos)*m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE)*8);
            uint32_t iStartOffset = iBitsToPos-((sizeof(FIELDTYPE)*8) * iStartPos);
            FIELDTYPE* pPos = m_pCounterArray + iStartPos;


            // THIS PREFETCH is necessary to avoid stupid status==0...
            PREFETCH = pPos[0];
            asm volatile("":::"memory");

            int iinc=0;

            uint status;
            if ((status = _xbegin ()) == _XBEGIN_STARTED) {

                //this->performIncrement(&incRet);
                // increment element

                bool elemEmpty=true;
                pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
                pSavedKey->bitAnd(m_mask_key_value);

                // check elements are equal
                for (uint8_t ei=0; ei < pSavedKey->m_iFields; ++ei)
                {
                    if (pSavedKey->m_pArray[ei] != 0)
                    {
                        elemEmpty = false;
                    }
                }

                if (!elemEmpty)
                {
                    _xabort(0xff);
                }

                pKeyReprobeShiftUBIGINT->m_pArray[0] = pKeyReprobeShiftUBIGINT->m_pArray[0] | 1;

                //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                (*pKeyReprobeShiftUBIGINT).copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);

                // transaction completes here
                _xend();
                asm volatile("":::"memory");


                // so we can find kmer start positions later without knowing the kmer
                m_iKmerStarts.setBit(iPos, 1);
                // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                m_setUsedPositions.insert(iPos);

                this->iAddCount += 1;
                bInserted = true;
                break;


            } else {

                ++iTotalAborts;
                if (_XABORT_CODE(status) == 0xff) {

                    /*
                     * The first transaction did not succeed because the position is not empty
                     */
                    // TODO where is the case kmer starts multiple positions later handled?
                    bool bIsKmerStart = m_iKmerStarts.getBit(iPos) == 1;
                    bool bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                    /*
                    if (bIsKmerStart) {
                        bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
                    } else {
                        bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
                    }
                     */

                    std::string sTestStr = "ACAGGACAAATACG";
                    UBigInt sTest = fromSequence(sTestStr, m_pPool);

                    bool isCurTest = kmer == sTest;

                    if (isCurTest)
                    {
                        std::cout << sTestStr << " " << bIsKmerStart << " " << bMatchesKey << std::endl;
                    }


                    if ((bIsKmerStart) && (bMatchesKey)) {

                        //std::cout << "going to increment " << std::endl;

                        uint8_t iAddStatus = 0;
                        uint8_t iOverflowStatus = 0;

                        while (iAddStatus == 0)
                        {
                            iAddStatus= this->incrementElement_tsx(kmer, iPos, basekey, iReprobes, false, verbose);
                        }

                        if (isCurTest)
                        {
                            std::cout << sTestStr << " after increment " << iAddStatus << " " << iOverflowStatus << std::endl;
                        }

                        if (iAddStatus == 2)
                        {
                            if (isCurTest)
                            {
                                std::cout << sTestStr << " before overflow" << std::endl;
                            }

                            while (iOverflowStatus == 0)
                            {
                                iOverflowStatus = this->handleOverflow_tsx(kmer, iPos, basekey, iReprobes, verbose);
                            }

                            if (isCurTest)
                            {
                                std::cout << sTestStr << " after  overflow increment " << iAddStatus << " " << iOverflowStatus << std::endl;
                            }
                        }

                        if (isCurTest)
                        {
                            std::cout << sTestStr << " before end " << iAddStatus << " " << iOverflowStatus << std::endl;
                        }


                        this->iAddCount += 1;
                        bInserted = true;

                        break;

                    } else {

                        //std::cout << "Skipping position " << iPos << ": kmer=" << bIsKmerStart << " match" << bMatchesKey << std::endl;

                        // transaction completed => incorrect position => reprobe
                        ++iReprobes;
                        continue;
                    }





                } else {
                    
                    ++iAborts;
                }

                if (iTotalAborts % 100000 == 0)
                {
                    std::cout << "aborts " << iTotalAborts << " inserts " << iAddCount << std::endl;
                }
            }

        }

        //std::cout << "Inserted element" << std::endl;

        if (!bInserted)
        {
            uint32_t maxPositions = 1 << m_iL;

            std::cerr << "Could not insert kmer " << kmer.to_string()<< std::endl;
            std::cerr << "Used fields: " << m_setUsedPositions.size() << std::endl;
            std::cerr << "Available fields: " << std::pow(2.0, m_iL) << std::endl;
            std::cerr << "k=" << m_iK << " l=" << (uint32_t) m_iL << " entry (key+value) bits=" << m_iKeyValBits << " storage bits=" << m_iStorageBits << std::endl;

            if (m_setUsedPositions.size() == maxPositions)
            {
                exit(42);
            }

            //this->addKmer(kmer);
        }

        if (!bInserted)
        {
            std::cout << "TSX ADDKMER RET NOT INSERTED" << std::endl;
        }

        return bInserted;

    }


    /**
     *
     * @param kmer
     * @param iPosition
     * @param key
     * @param reprobes
     * @param verbose
     * @return
     */
    uint8_t incrementElement_key_value(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false)
    {

        /**
         *
         * ONLY INCREMENT VALUE PART => ret 2 on overflow
         *
         */

        std::string stest = m_mask_value_key.to_string();
        stest = m_mask_key_value.to_string();

        //std::cout << "in inc key value" << std::endl;

        for (uint8_t iTSXRetries=0; iTSXRetries < 10; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            PREFETCH = pPos[0];

            TSX::tsx_val_t value = UBigInt(m_iStorageBits, true, this->m_pPool);

            SBIGINT* pZero = getEmptySBIGINT(m_iStorageBits);
            SBIGINT* pKeyVal = getEmptySBIGINT(m_iKeyValBits);
            SBIGINT* pTest = getEmptySBIGINT(m_iKeyValBits);
            SBIGINT* pValue = getEmptySBIGINT(m_iKeyValBits);


            getFromMemory(pTest, iStartOffset, m_iKeyValBits, pPos);

            UBigInt saved_b = fromStructToClass(pTest, m_pPool);
            //std::cout << "before inc position " << " " << iPosition << " " << saved_b.to_string() << std::endl;


            asm volatile("":: :"memory");

            bool elemEmpty = false;

            uint status = _xbegin();
            if (status == _XBEGIN_STARTED) {

                //this->performIncrement(&incRet);
                // increment element

                //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
                getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);
                getFromMemory(pTest, iStartOffset, m_iKeyValBits, pPos);

                //value = (*pSavedKey) &  m_mask_value_key;

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                {
                    pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_value_key.m_pArray[i];
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_key_value.m_pArray[i];
                }

                bitComplement(pValue);
                elemEmpty = isZero(pValue, m_iStorageBits);


                // check whether overflow occurs, or not
                // elemEmpty == true => OVERFLOW OCCURS!
                //elemEmpty = false;//(~value).isZero();

                if (elemEmpty) {

                    for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                    {
                        pValue->pdata[i] = 0;
                    }

                } else {
                    bitComplement(pValue);

                    add_simpleSBIGINT(pValue);

                    for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                    {
                        // keeps key (func+reprobe), nulls value
                        pKeyVal->pdata[i] = pKeyVal->pdata[i] | pValue->pdata[i];
                    }
                }

                //pSavedKey->bitAnd(m_mask_key_value);
                //pSavedKey->bitOr(value);

                //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                //(*pSavedKey).copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);
                storeInMemory(pKeyVal->pdata, pPos, iStartOffset, m_iKeyValBits, pKeyVal->iFieldSize);

                // transaction completes here
                _xend();
                asm volatile("":: :"memory");


                UBigInt value = fromStructToClass(pValue, m_pPool);
                UBigInt saved = fromStructToClass(pKeyVal, m_pPool);

                //std::cout << "after inc position " << elemEmpty << " " << iPosition << " " << value.to_string() << " " << saved.to_string() << std::endl;

                if (elemEmpty)
                {
                    //std::cout << "elemEmpty" << std::endl;
                    // OVERFLOW in VALUE part
                    return 2;
                }

                // NO OVERFLOW
                return 1;

            } else {

                ++iTotalAborts;

                /**
                 *
                 * TODO abort stats
                 */

            }


        }

        return 0;

    }

    /**
     *
     * @param kmer
     * @param iPosition
     * @param key
     * @param reprobes
     * @param verbose
     * @return
     */
    uint8_t incrementElement_func(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false)
    {

        /**
         *
         * ONLY INCREMENT FUNC PART => ret 2 on overflow
         *
         */

        for (uint8_t iTSXRetries=0; iTSXRetries < 10; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t* pSavedKey = &savedkey;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            PREFETCH = pPos[0];
            asm volatile("":: :"memory");

            SBIGINT* pZero = getEmptySBIGINT(m_iStorageBits);
            SBIGINT* pKeyVal = getEmptySBIGINT(m_iKeyValBits);
            SBIGINT* pValue = getEmptySBIGINT(m_iKeyValBits);

            m_mask_kv_reprobevalue;

            bool elemEmpty = false;

            uint status = _xbegin();
            if (status == _XBEGIN_STARTED) {

                //this->performIncrement(&incRet);
                // increment element

                getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);
                //value = (*pSavedKey) &  m_mask_value_key;

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                {
                    pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_kv_func.m_pArray[i];
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_kv_reprobevalue.m_pArray[i];
                }

                shiftRight(pValue, m_iL+m_iStorageBits);

                //pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
                //TSX::tsx_val_t funcValue = ((*pSavedKey) >> (m_iL + m_iStorageBits));

                // check whether overflow occurs, or not
                // elemEmpty == true => OVERFLOW OCCURS!
                //elemEmpty = (~funcValue).isZero();

                bitComplement(pValue);
                elemEmpty = isZero(pValue, 2*m_iK+m_iStorageBits);

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                {
                    // keeps key (value+reprobe), nulls value
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_kv_reprobevalue.m_pArray[i];
                }

                //pSavedKey->bitAnd(m_mask_kv_reprobevalue); // clears func bits

                if (elemEmpty) {
                    //is already null...

                    for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                    {
                        pValue->pdata[i] = 0;
                    }

                } else {
                    bitComplement(pValue);

                    add_simpleSBIGINT(pValue);
                    shiftLeft(pValue, m_iL+m_iStorageBits);

                }

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
                {
                    // keeps key (func+reprobe), nulls value
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] | pValue->pdata[i];
                }

                //pSavedKey->copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);
                storeInMemory(pKeyVal->pdata, pPos, iStartOffset, m_iKeyValBits, pKeyVal->iFieldSize);


                // transaction completes here
                _xend();
                asm volatile("":: :"memory");

                if (elemEmpty)
                {
                    // OVERFLOW in VALUE part
                    return 2;
                }

                // NO OVERFLOW
                return 1;

            } else {

                ++iTotalAborts;

                /**
                 *
                 * TODO abort stats
                 */

            }


        }

        return 0;

    }


    /**
     *
     * @param iPosition position in array
     * @param key key to look for
     * @param reprobes amount reprobes used
     * @param bKeyIsValue if true, the key part belongs to value
     * @return 1 if no overflow, 0 not added and 2 if overflow occured
     *
     * @details this function already handles overflows
     */
    uint8_t incrementElement_tsx(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool bKeyIsValue, bool verbose=false)
    {

        std::string sTestStr = "AACAGCGTTTCGTC";
        UBigInt sTest = fromSequence(sTestStr, m_pPool);

        bool isCurTest = kmer == sTest;
        bool handleFuncOverflow = false;

        //std::cout << "increment element tsx " << iPosition << std::endl;

        // try to increment value

        for (uint8_t iIncrementTries=0; iIncrementTries < 10; ++iIncrementTries)
        {
            uint8_t iIncrementState = this->incrementElement_key_value(kmer, iPosition, key, reprobes, verbose);

            if (iIncrementState == 0)
            {
                //std::cout << "simple increment state 0" << std::endl;
                continue;

            } else if (iIncrementState == 1)
            {

                if (isCurTest)
                {
                    std::cout << "return 1 after increment" << std::endl;
                }

                // value increment, everything done
                return 1;
            } else if (iIncrementState == 2)
            {
                // value incremented, but overflow occurred
                if (bKeyIsValue)
                {

                    if (isCurTest)
                    {
                        std::cout << "increment returns 2, break" << std::endl;
                    }

                    handleFuncOverflow = true;
                    break;

                } else {

                    if (isCurTest)
                    {
                        std::cout << "increment returns 2, no break" << std::endl;
                    }

                    //std::cout << "value overflow in pos " << iPosition << std::endl;
                    return 2;
                }
            }


        }




        // the value part had an overflow => maybe the func part can handle this => only if key is value
        if (handleFuncOverflow && bKeyIsValue)
        {

            if (isCurTest)
            {
                std::cout << "cur test overflow in increment" << std::endl;
            }
            // try to avoid overflow by incrementing func part
            for (uint8_t iFuncIncrementTries=0; iFuncIncrementTries < 10; ++iFuncIncrementTries) {

                uint8_t iIncrementFuncState = this->incrementElement_func(kmer, iPosition, key, reprobes, verbose);


                if (iIncrementFuncState == 0)
                {
                    continue;

                } else if (iIncrementFuncState == 1)
                {
                    // value increment, everything done
                    return 1;
                } else if (iIncrementFuncState == 2)
                {
                    // value incremented, but overflow occurred => this was already func => return overflow

                    uint8_t iOverflowReturn = this->handleOverflow_tsx(kmer, iPosition, key, reprobes, verbose);

                    if (iOverflowReturn == 0)
                    {
                        std::cout << "overflow return 0" << std::endl;
                    }

                    return 1;

                }


            }


        } else {
            if (isCurTest)
            {
                std::cout << "cur test no overflow in increment" << std::endl;
            }
        }

        exit(33);

        // should never reach this => not inserted
        return 0;
    }

    virtual uint8_t handleOverflow_tsx(TSX::tsx_kmer_t& kmer, uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, bool verbose=false)
    {
        uint32_t iPerformedReprobes = 0;
        uint8_t iThreadID = omp_get_thread_num();

        UBigInt uBigOne = UBigInt(1, m_pPool);

        std::string sTestStr = "AACAGCGTTTCGTC";
        UBigInt sTest = fromSequence(sTestStr, m_pPool);
        bool isCurTest = sTest == kmer;


        if (isCurTest)
        {
            std::cout << "overflow add before " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << " with reprobe " << iReprobe << std::endl;
        }

        while (iPerformedReprobes < m_iMaxReprobes)
        {

            iPerformedReprobes += 1;

            if (isCurTest)
            {
                std::cout << "overflow add before " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << " with reprobe " << iReprobe << " perf reprobe " << iPerformedReprobes << std::endl;
            }

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe + iPerformedReprobes);

            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t* pSavedKey = &savedkey;

            TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
            TSX::tsx_key_t updkey = TSX::tsx_key_t(2*m_iK, true, this->m_pPool);
            updkey = updkey | UBigInt(iPerformedReprobes, this->m_pPool);

            TSX::tsx_keyval_t oKeyVal(updkey, 2 * m_iK + m_iStorageBits);
            oKeyVal = (oKeyVal << m_iStorageBits);
            oKeyVal.resize(m_iKeyValBits);

            uint64_t iBitsToPos = (iPos)*m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE)*8);
            uint32_t iStartOffset = iBitsToPos-((sizeof(FIELDTYPE)*8) * iStartPos);
            FIELDTYPE* pPos = m_pCounterArray + iStartPos;


            // THIS PREFETCH is necessary to avoid stupid status==0...
            PREFETCH = pPos[0];
            asm volatile("":::"memory");

            uint status = _xbegin();
            if(status == _XBEGIN_STARTED) {

                //this->performIncrement(&incRet);
                // increment element

                bool elemEmpty = true;
                pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
                pSavedKey->bitAnd(m_mask_key_value);

                // check elements are equal
                for (uint8_t ei = 0; ei < pSavedKey->m_iFields; ++ei) {
                    if (pSavedKey->m_pArray[ei] != 0) {
                        elemEmpty = false;
                    }
                }

                if (!elemEmpty) {
                    _xabort(0xff);
                }

                oKeyVal.m_pArray[0] = oKeyVal.m_pArray[0] | 1;

                //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                oKeyVal.copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);

                // transaction completes here
                _xend();
                asm volatile("":::"memory");

                if (isCurTest)
                {
                    std::cout << "overflow add empty  after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;
                }

                //std::cout << "overflow add empty  after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;

                return 1;

            } else {

                ++iTotalAborts;

                if (_XABORT_CODE(status) == 0xff) {

                    /*
                     * The first transaction did not succeed because the position is not empty
                     */
                    // TODO where is the case kmer starts multiple positions later handled?
                    bool bIsKmerStart = m_iKmerStarts.getBit(iPos) == 1;

                    /*
                    if (bIsKmerStart)
                    {
                        continue;
                    }
                     */

                    // the reprobe part must match the number of reprobes back to the previous entry!
                    bool bMatchesKey = positionMatchesReprobe(iPos, basekey,
                                                              iPerformedReprobes); // was iReprobe+iPerformedReprobes

                    if ((bIsKmerStart) || (!bMatchesKey)) {
                        //std::cout << "overflow Skipping position " << iPos << ": kmer=" << bIsKmerStart << " reprobes=" << iReprobe << " match=" << bMatchesKey << std::endl;

                        continue;
                    }

                    //std::cout << "overflow add before inc " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                    uint8_t incremented = this->incrementElement_tsx(kmer, iPos, basekey, iPerformedReprobes, true, verbose); // was iReprobe+iPerformedReprobes

                    if (incremented == 0) {

                        // repeat -> hence undo reprobe
                        --iPerformedReprobes;

                        // redo increment ...
                        continue;
                    }
                    //std::cout << "overflow add after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                    return 1;



                } else {

                    //std::cout << "error code " << status << std::endl;
                    --iPerformedReprobes;

                    ++iAborts;
                }

                if (iTotalAborts % 10 == 0) {
                    std::cout << "aborts " << iTotalAborts << " inserts " << iAddCount << std::endl;
                }

            }
        }


        std::cout << "TSX OVERFLOW RET 0" << std::endl;

        return 0;
    }

    virtual bool addKmer(TSX::tsx_kmer_t& kmer, bool verbose=false)
    {
        return this->addKmer_tsx(kmer, verbose);
    }

protected:

    void initialiseLocks()
    {
    }



    /**
     *
     * @param iThreadID thread id for which the lock is acquired
     * @param iArrayPos array position for which the lock is acquired
     * @return True if lock successfully acquired, False otherwise
     */
    bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos)
    {
        return true;
    }


    /**
     *
     * unlocks all locked array positions for given thread
     *
     * @param iThreadID
     */
    void unlock_thread(uint8_t iThreadID)
    {
    }

    bool releaseLock(uint8_t iThreadID, uint64_t iPos)
    {
        return true;
    }

};


#endif //TSXCOUNT_TSXHASHMAPTSXPERF_H
