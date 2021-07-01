//
// Created by mjopp on 25/02/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPCAS_H
#define TSXCOUNT_TSXHASHMAPCAS_H

#include "TSXHashMap.h"
#include "TSXHashMapPerf.h"
#include <pthread.h>
#include <stack>

#include <stdbool.h>

#include <immintrin.h>
#include <unistd.h>
#include <atomic>
#include <cmath>
#include <bitset>
#include <src/tsxutils/SBigInt.h>

#ifdef __cplusplus
#include <atomic>
#else
#include <stdatomic.h>
#endif

void getFromMemory(SBIGINT::SBIGINT* pRes, SBIGINT::SBIGINT* pFields, uint8_t iOffset, uint16_t iKeyValBits, FIELDTYPE* pData )
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

void storeInMemoryOriginal(SBIGINT::SBIGINT* pSrcKeyVal, SBIGINT::SBIGINT* pOrigValues, FIELDTYPE* pDest,uint32_t iBitStart, uint16_t iBitCount, uint8_t iFieldSize) {

    // how many bits can I copy into first, not fully occupied field?
    uint8_t iOffset = iBitStart % iFieldSize;

    // TODO copy first field
    FIELDTYPE iThisVal = pSrcKeyVal->pdata[0] << iOffset;
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
        iPartRight = pSrcKeyVal->pdata[i] >> (iFieldSize-iOffset);
        iPartLeft = pSrcKeyVal->pdata[i+1] << iOffset;
        iThisVal = iPartLeft | iPartRight;

        pDest[i+1] = iThisVal;

        iBitsRemaining -= iFieldSize;
    }

    if (iBitsRemaining > 0)
    {
        if (iBitsRemaining <= iOffset)
        {

            // we do not need a rest from the next field
            iPartRight = pSrcKeyVal->pdata[fields.quot] << (iFieldSize-iOffset);
            iPartRight = iPartRight << (iOffset-iBitsRemaining);
            iPartRight = iPartRight >> (iOffset-iBitsRemaining);
            iPartRight = iPartRight >> iFieldSize-iOffset;

            iOldVal = pDest[fields.quot+1];
            iOldVal = iOldVal >> iBitsRemaining;
            iOldVal = iOldVal << iBitsRemaining;

            pDest[fields.quot+1] = iOldVal | iPartRight;

        } else {

            // we need a rest from the next field
            iPartRight = pSrcKeyVal->pdata[fields.quot] >> (iFieldSize-iOffset);

            uint32_t iRemaining = iBitsRemaining - iOffset;

            // need the rightmost iRemaining bits
            iPartLeft = pSrcKeyVal->pdata[fields.quot+1] << (iFieldSize-iRemaining); // clears all but iRemaining bits
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


uint8_t storeInMemoryCAS(SBIGINT::SBIGINT* pSrcKeyVal, SBIGINT::SBIGINT* pOrigValues,SBIGINT::SBIGINT* pEmpty, FIELDTYPE* pDest,uint32_t iBitStart, uint16_t iBitCount, uint8_t iFieldSize, uint32_t iKeyValueBits) {

    // how many bits can I copy into first, not fully occupied field?
    uint8_t iOffset = iBitStart % iFieldSize;

    //??000000 0000000 00000000 000000??

    SBIGINT::clear(pEmpty);
    //00001111 11111111 111111111 11111111
    SBIGINT::setOnes(pEmpty, iKeyValueBits);

    //00111111 11111111 111111111 11111100
    SBIGINT::shiftLeft(pEmpty, iOffset);

    //11000000 0000000 0000000 00000011
    SBIGINT::bitComplement(pEmpty);

    //____1111 11111111 1111111 11111111
    //__111111 11111111 1111111 11111100
    SBIGINT::shiftLeft(pSrcKeyVal, iOffset);

    //??000000 0000000 00000000 000000??
    SBIGINT::bitAnd(pEmpty, pOrigValues);

    //??111111 11111111 11111111 111111??
    SBIGINT::bitOr(pEmpty, pSrcKeyVal);

    uint32_t iRemainingBits = iKeyValueBits;

    for (uint8_t i = 0; i < pOrigValues->iFields; ++i)
    {

        if (iRemainingBits == 0)
        {
            break;
        }


        //std::bitset<8> bs_orig(pDest[i]);
        //std::bitset<8> bs_desired(pEmpty->pdata[i]);

        //std::cout << bs_orig << " , " << bs_desired << " --> ";

        //std::atomic_compare_exchange_weak_explicit(&oDest, pOrig, pEmpty->pdata[i], std::memory_order_release, std::memory_order_relaxed);
        bool CASresult = __atomic_compare_exchange_n(pDest+i, pOrigValues->pdata+i, pEmpty->pdata[i], false, __ATOMIC_RELAXED,__ATOMIC_RELAXED);

        //std::bitset<8> bs_orig2(pDest[i]);
        //std::bitset<8> bs_desired2(pEmpty->pdata[i]);

        //std::cout << bs_orig2 << " , " << bs_desired2 << std::endl;

        //bool casResult = __atomic_compare_exchange_n(pDest+i, pOrigValues->pdata+i, pEmpty->pdata[i], false, __ATOMIC_RELAXED,__ATOMIC_RELAXED);

        iRemainingBits -= (pOrigValues->iFieldSize);

        if (!CASresult)
        {
            if (i == 0)
            {
                // very first try did not succeed => retry
                return 0;
            }

            SBIGINT::clear(pEmpty);
            //00001111 11111111 111111111 11111111
            SBIGINT::setOnes(pEmpty, iKeyValueBits);

            //00111111 11111111 111111111 11111100
            SBIGINT::shiftLeft(pEmpty, iOffset);

            //11000000 0000000 0000000 00000011
            SBIGINT::bitComplement(pEmpty);

            for (uint8_t j = i; j > 1; --j)
            {
                // 111111?? -> 000000??
                pDest[j] = (pDest[j] & pEmpty->pdata[j]) | (pOrigValues->pdata[j] & ~pEmpty->pdata[j]);
            }

            std::cout << "rolled back " << (int) i  << " with remains " << (int) iRemainingBits << std::endl;

            // need to care about roleback, because most likely the neighbouring field has changed

            return 0;
        }


    }

    return 1;

}


class TSXHashMapCAS : public TSXHashMapPerf {

public:

    TSXHashMapCAS(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
            : TSXHashMapPerf(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);

    }


    virtual ~TSXHashMapCAS()
    {
        free(m_pLocked);


        free(m_pTMP_KEYVAL);
        free(m_pTMP_VALUE);
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
    virtual bool addKmer(TSX::tsx_kmer_t& kmer, bool verbose=false, bool noPrimaryAddition=false) override
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
            TSX::tsx_keyval_t* pKeyReprobeShiftUBIGINT = &key_reprobe_shift;

            uint64_t iBitsToPos = (iPos)*m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE)*8);
            uint32_t iStartOffset = iBitsToPos-((sizeof(FIELDTYPE)*8) * iStartPos);
            FIELDTYPE* pPos = m_pCounterArray + iStartPos;


            // THIS PREFETCH is necessary to avoid stupid status==0...
            asm volatile("":::"memory");

            SBIGINT::SBIGINT *pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT *pValue = m_pTMP_VALUE[iThreadID];
            SBIGINT::SBIGINT *pOrigValue = m_pTMP_ORIGINAL[iThreadID];
            SBIGINT::SBIGINT *pClear = m_pTMP_CLEAR[iThreadID];

            //this->performIncrement(&incRet);
            // increment element

            //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            SBIGINT::getFromMemory(pKeyVal, pOrigValue, iStartOffset, m_iKeyValBits, pPos);


            //this->performIncrement(&incRet);
            // increment element

            //pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            //pSavedKey->bitAnd(m_mask_key_value);


            for (uint8_t i = 0; i < pKeyVal->iFields; ++i) {
                pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_value_key.m_pArray[i];
                pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_key_value.m_pArray[i];
            }

            bool elemEmpty = SBIGINT::isZero(pKeyVal, m_iKeyValBits);

            if (elemEmpty) {
                for (uint8_t i = 0; i < pKeyReprobeShiftUBIGINT->m_iFields; ++i) {
                    pKeyVal->pdata[i] = pKeyReprobeShiftUBIGINT->m_pArray[i];
                }

                pKeyVal->pdata[0] = pKeyVal->pdata[0] | 1;

                //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                //(*pKeyReprobeShiftUBIGINT).copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);

                // transaction completes here
                uint8_t CASresult = storeInMemoryCAS(pKeyVal, pOrigValue, pClear, pPos, iStartOffset, m_iKeyValBits, sizeof(FIELDTYPE) * 8, m_iKeyValBits);

                if (CASresult == 0)
                {
                    ++iTotalAborts;

                    if (iTotalAborts % 100000 == 0)
                    {
                        std::cout << "aborts " << iTotalAborts << " inserts " << iAddCount << std::endl;
                    }

                    SBIGINT::getFromMemory(pKeyVal, pOrigValue, iStartOffset, m_iKeyValBits, pPos);

                    // try again
                    continue;

                } else {

                    //uint64_t iCount = this->getKmerCount(kmer, false).toUInt();
                    //std::cout << "Added Empty: " << kmer.to_string() << " in pos " << iPos << " " << iCount << std::endl;


                    // so we can find kmer start positions later without knowing the kmer
                    m_iKmerStarts.setBit(iPos, 1);
                    // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                    m_setUsedPositions.insert(iPos);

                    this->iAddCount += 1;
                    bInserted = true;
                    break;
                }


            } else {

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

                std::string sTestStr = toSequence(kmer);
                UBigInt sTest = fromSequence(sTestStr, m_pPool);

                bool isCurTest = kmer == sTest;
                isCurTest = false;

                if (isCurTest)
                {
                    std::cout << sTestStr << " " << bIsKmerStart << " " << bMatchesKey << std::endl;
                }


                if ((bIsKmerStart) && (bMatchesKey)) {

                    //std::cout << "going to increment " << std::endl;

                    uint8_t iAddStatus = 0;
                    uint8_t iOverflowStatus = 0;

                    uint64_t beforeCount = this->getKmerCount(kmer).toUInt();

                    if (isCurTest)
                    {
                        std::cout << sTestStr << " before increment " << (int) iAddStatus << " " << (int) iOverflowStatus << " count: " << beforeCount << std::endl;
                    }



                    while (iAddStatus == 0)
                    {
                        iAddStatus= this->incrementElement_new(kmer, iPos, basekey, iReprobes, iReprobes, false, verbose);
                    }

                    /*
                    uint64_t afterCount = this->getKmerCount(kmer).toUInt();

                    if (isCurTest)
                    {
                        std::cout << sTestStr << " after increment " << (int) iAddStatus << " " << (int) iOverflowStatus << " count: " << this->getKmerCount(kmer).toUInt() << std::endl;
                    }

                    if ((iAddStatus != 2) && (afterCount-beforeCount != 1))
                    {
                        std::cout << "ERROR JUST OCCURRED!" << std::endl;
                    }
                     */

                    if (iAddStatus == 2)
                    {
                        if (isCurTest)
                        {
                            std::cout << sTestStr << " before overflow" << std::endl;
                        }

                        while (iOverflowStatus == 0)
                        {
                            iOverflowStatus = this->handleOverflow(kmer, iPos, basekey, iReprobes, iReprobes, verbose);
                        }

                        if (isCurTest)
                        {
                            std::cout << sTestStr << " after overflow increment " << (int) iAddStatus << " " << (int) iOverflowStatus << " count: " << this->getKmerCount(kmer).toUInt() << std::endl;
                        }
                    }

                    if (isCurTest)
                    {
                        std::cout << sTestStr << " before end " << (int)iAddStatus << " " << (int)iOverflowStatus << std::endl;
                    }


                    this->iAddCount += 1;
                    bInserted = true;

                    break;

                }
            }

            ++iReprobes;

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
    virtual uint8_t incrementElement_key_value(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false) override {

        /**
         *
         * ONLY INCREMENT VALUE PART => ret 2 on overflow
         *
         */

        std::string stest = m_mask_value_key.to_string();
        stest = m_mask_key_value.to_string();

        uint8_t iThreadID = omp_get_thread_num();


        //std::cout << "in inc key value" << std::endl;

        for (uint8_t iTSXRetries = 0; iTSXRetries < 10; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            TSX::tsx_val_t value = UBigInt(m_iStorageBits, true, this->m_pPool);

            asm volatile("":: :"memory");

            bool elemEmpty = false;

            SBIGINT::SBIGINT *pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT *pValue = m_pTMP_VALUE[iThreadID];
            SBIGINT::SBIGINT *pOrigValue = m_pTMP_ORIGINAL[iThreadID];
            SBIGINT::SBIGINT *pClear = m_pTMP_CLEAR[iThreadID];
            //this->performIncrement(&incRet);
            // increment element

            //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            SBIGINT::getFromMemory(pKeyVal, pOrigValue, iStartOffset, m_iKeyValBits, pPos);

            //value = (*pSavedKey) &  m_mask_value_key;

            for (uint8_t i = 0; i < pKeyVal->iFields; ++i) {
                pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_value_key.m_pArray[i];
                pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_key_value.m_pArray[i];
            }

            bitComplement(pValue);
            elemEmpty = isZero(pValue, m_iStorageBits);


            // check whether overflow occurs, or not
            // elemEmpty == true => OVERFLOW OCCURS!
            //elemEmpty = false;//(~value).isZero();

            if (elemEmpty) {

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i) {
                    pValue->pdata[i] = 0;
                }

            } else {
                bitComplement(pValue);

                add_simpleSBIGINT(pValue);

                for (uint8_t i = 0; i < pKeyVal->iFields; ++i) {
                    // keeps key (func+reprobe), nulls value
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] | pValue->pdata[i];
                }
            }

            //pSavedKey->bitAnd(m_mask_key_value);
            //pSavedKey->bitOr(value);

            //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
            //(*pSavedKey).copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);
            uint8_t CASReturn = storeInMemoryCAS(pKeyVal, pOrigValue, pClear, pPos, iStartOffset, m_iKeyValBits,
                                                 pKeyVal->iFieldSize, m_iKeyValBits);

            if (CASReturn == 1) {
                // CAS Worked - yes!

                //UBigInt value = fromStructToClass(pValue, m_pPool);
                //UBigInt saved = fromStructToClass(pKeyVal, m_pPool);

                //std::cout << "after inc position " << elemEmpty << " " << iPosition << " " << value.to_string() << " " << saved.to_string() << std::endl;

                if (elemEmpty) {
                    //std::cout << "elemEmpty" << std::endl;
                    // OVERFLOW in VALUE part
                    return 2;
                }

                // NO OVERFLOW
                return 1;
            } else {

                std::cout << "INC KEY VAL CAS 0" << std::endl;

                //CAS DID NOT WORK => REPEAT!
                ++iTotalAborts;

                continue;
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
    virtual uint8_t incrementElement_func(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false) override
    {

        /**
         *
         * ONLY INCREMENT FUNC PART => ret 2 on overflow
         *
         */
        uint8_t iThreadID = omp_get_thread_num();

        for (uint8_t iTSXRetries=0; iTSXRetries < 10; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t* pSavedKey = &savedkey;

            asm volatile("":: :"memory");

            SBIGINT::SBIGINT* pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT* pValue = m_pTMP_VALUE[iThreadID];
            SBIGINT::SBIGINT *pOrigValue = m_pTMP_ORIGINAL[iThreadID];
            SBIGINT::SBIGINT*pClear = m_pTMP_CLEAR[iThreadID];

            m_mask_kv_reprobevalue;

            bool elemEmpty = false;



            //this->performIncrement(&incRet);
            // increment element

            SBIGINT::getFromMemory(pKeyVal, pOrigValue, iStartOffset, m_iKeyValBits, pPos);
            //value = (*pSavedKey) &  m_mask_value_key;

            for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
            {
                pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_kv_func.m_pArray[i];
                pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_kv_reprobevalue.m_pArray[i];

            }


            //10______________________________
            shiftRight(pValue, m_iL+m_iStorageBits);
            //00000000000000000000000000000010

            //pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            //TSX::tsx_val_t funcValue = ((*pSavedKey) >> (m_iL + m_iStorageBits));

            // check whether overflow occurs, or not
            // elemEmpty == true => OVERFLOW OCCURS!
            //elemEmpty = (~funcValue).isZero();


            //00000000000000000000000000000010
            bitComplement(pValue);
            //11111111111111111111111111111101
            elemEmpty = isZero(pValue, 2*m_iK-m_iL);
            //11111111111111111111111111111101
            bitComplement(pValue);
            //00000000000000000000000000000010

            //pSavedKey->bitAnd(m_mask_kv_reprobevalue); // clears func bits

            if (elemEmpty) {

                for (uint8_t i = 0; i < pValue->iFields; ++i) {
                    pValue->pdata[i] = 0;
                }

                //is already null...
                //00000000000000000000000000000000
            } else {

                //00000000000000000000000000000010
                add_simpleSBIGINT(pValue);
                //00000000000000000000000000000011
                shiftLeft(pValue, m_iL+m_iStorageBits);
                //11000000000000000000000000000000
            }

            for (uint8_t i = 0; i < pKeyVal->iFields; ++i)
            {
                // keeps key (func+reprobe), nulls value
                pKeyVal->pdata[i] = pKeyVal->pdata[i] | pValue->pdata[i];
            }

            //pSavedKey->copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);


            uint8_t CASresult = storeInMemoryCAS(pKeyVal, pOrigValue, pClear, pPos, iStartOffset, m_iKeyValBits, pKeyVal->iFieldSize, m_iKeyValBits);

            if (CASresult != 0)
            {
                /*
                UBigInt oVal =fromStructToClass(pValue, m_pPool);
                UBigInt oKV =fromStructToClass(pKeyVal, m_pPool);
                UBigInt oOV = fromStructToClass(pOrigValue, m_pPool);


                std::cout << "FUNC INC " << elemEmpty << " " << oVal.to_string() << " pos " << iPosition << " " << "" << std::endl;
                std::cout << "FUNC INC " << elemEmpty << " " << oKV.to_string() << " pos " << iPosition << " " << "" << std::endl;
                std::cout << "FUNC INC " << elemEmpty << " " << oOV.to_string() << " pos " << iPosition << " " << "" << std::endl;
                 */

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

                continue;
            }

        }

        return 0;

    }

/*

    */
/**
     *
     * @param iPosition position in array
     * @param key key to look for
     * @param reprobes amount reprobes used
     * @param bKeyIsValue if true, the key part belongs to value
     * @return 1 if no overflow, 0 not added and 2 if overflow occured
     *
     * @details this function already handles overflows
     *//*

    uint8_t incrementElement_tsx(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, uint32_t iInitialReprobes, bool bKeyIsValue, bool verbose=false)
    {

        std::string sTestStr = "TTCCATTCCATTCC";
        UBigInt sTest = fromSequence(sTestStr, m_pPool);

        bool isCurTest = kmer == sTest;
        isCurTest = false;
        bool handleFuncOverflow = false;

        //std::cout << "increment element tsx " << iPosition << std::endl;

        // try to increment value
        uint8_t iIncrementTries = 0;
        uint8_t iIncrementState = 0;

        //for (iIncrementTries=0; iIncrementTries < 100; ++iIncrementTries)
        while(iIncrementState == 0)
        {
            iIncrementState = this->incrementElement_key_value(kmer, iPosition, key, reprobes, verbose);

            if (iIncrementState == 0)
            {
                if (isCurTest)
                {
                    std::cout << "simple increment return 0, pos="<< iPosition << std::endl;
                }
                continue;

            } else if (iIncrementState == 1)
            {

                if (isCurTest)
                {
                    std::cout << "return 1 after increment, pos="<< iPosition << std::endl;
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
                        std::cout << "increment returns 2, break, pos="<< iPosition << std::endl;
                    }

                    handleFuncOverflow = true;
                    break;

                } else {

                    if (isCurTest)
                    {
                        std::cout << "increment returns 2, no break, pos="<< iPosition << std::endl;
                    }

                    //std::cout << "value overflow in pos " << iPosition << std::endl;
                    return 2;
                }
            }


        }



        uint8_t iFuncIncrementTries =0;
        uint8_t iIncrementFuncState = -1;
        // the value part had an overflow => maybe the func part can handle this => only if key is value
        if (handleFuncOverflow && bKeyIsValue)
        {

            if (isCurTest)
            {
                std::cout << "cur test overflow in increment" << " pos="<< iPosition<< std::endl;
            }
            // try to avoid overflow by incrementing func part
            for (iFuncIncrementTries=0; iFuncIncrementTries < 10; ++iFuncIncrementTries) {

                iIncrementFuncState = this->incrementElement_func(kmer, iPosition, key, reprobes, verbose);

                if (isCurTest)
                {
                    std::cout << "cur test after overflow in func " << (int) iIncrementFuncState << " pos="<< iPosition<< std::endl;
                }

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

                    uint8_t iOverflowReturn = this->handleOverflow_tsx(kmer, iPosition, key, reprobes, iInitialReprobes, verbose);

                    if (iOverflowReturn == 0)
                    {
                        std::cout << "overflow return 0" << " pos="<< iPosition<< std::endl;
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


        std::cout << handleFuncOverflow << " " << bKeyIsValue << " " << (int) iIncrementTries << " " << (int) iIncrementState << " " << (int)iFuncIncrementTries << " " << (int)iIncrementFuncState << std::endl;
        this->print_stats();

        exit(33);

        // should never reach this => not inserted
        return 0;
    }
*/

    /*
    virtual UBigInt findOverflowCounts(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, bool verbose)
    {
        bool bHandled = false;
        uint32_t iPerformedReprobes = 0;

        UBigInt oReturn(0, m_pPool);
        uint32_t iRequiredBits = 0;
        uint64_t origPos = iPos;

        uint32_t iInitialReprobes = iReprobe;

        while (iPerformedReprobes < 10)
        {

            iPerformedReprobes += 1;

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe+iPerformedReprobes);

            //std::cout << "OVFL Looking at pos " << iPos << " " << iReprobe << std::endl;

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (!bEmpty)
            {

                if (m_iKmerStarts.getBit(iPos) == 1)
                    continue;
                //std::cout << "OVFL B Match Key " << iPos << " " << iReprobe << std::endl;

                // the reprobe part must match the number of reprobes back to the previous entry!
                bool bMatchesKey = positionMatchesOverflowReprobe(iPos, basekey, iReprobe, iPerformedReprobes); // was iReprobe

                if (!bMatchesKey)
                {

                    if (verbose)
                    {
                        TSX::tsx_keyval_t elem = this->getElement(iPos);

                        std::cout << "OVFL Unmatched Key; orig pos " << origPos << " test pos " << iPos << " initial reprobes " << iInitialReprobes << " reprobe " << iReprobe << " perf reprobe " << iPerformedReprobes << " all_reprobes " << iPerformedReprobes+iReprobe << std::endl;
                        std::cout << elem.to_string() << std::endl;
                    }

                    continue;
                }



                if (verbose)
                {
                    std::cout << "OVFL A Match Key " << iPos << " init_reprobe " << iInitialReprobes << " reprobes " << iReprobe << " p_reprobes " << iPerformedReprobes << " " << iPerformedReprobes+iReprobe << std::endl;
                }

                TSX::tsx_keyval_t elem = this->getElement(iPos);

                UBigInt posValue = this->getFuncValFromKeyVal(elem);

                //std::cout << "PV " << posValue.to_string() << std::endl;

                uint32_t iOldRequired = iRequiredBits;
                iRequiredBits += 2*m_iK - m_iL + m_iStorageBits;
                posValue.resize(iRequiredBits);

                oReturn = (posValue << ( iOldRequired )) | oReturn;


                iReprobe += iPerformedReprobes;
                // reset performed reprobes as this should indicate number of reprobes needed!
                // now we can try to find further matching positions :)
                iPerformedReprobes = 0;

            }

        }

        // no more matching position as max number of reprobes reached without finding a matching one
        return oReturn;

    }
     */

    virtual uint8_t handleOverflow(TSX::tsx_kmer_t& kmer, uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, uint32_t iInitialReprobes, bool verbose=false) override
    {
        uint32_t iPerformedReprobes = 0;
        uint8_t iThreadID = omp_get_thread_num();

/*        UBigInt uBigOne = UBigInt(1, m_pPool);

        std::string sTestStr = "AACAGCGTTTCGTC";
        UBigInt sTest = fromSequence(sTestStr, m_pPool);
        bool isCurTest = sTest == kmer;
        isCurTest = false;

        sTestStr = toSequence(kmer);


        if (isCurTest)
        {
            std::cout << "overflow add before " << sTestStr << " " << kmer.to_string() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << " with reprobe " << iReprobe << std::endl;
        }*/

        while (iPerformedReprobes < m_iMaxReprobes)
        {

            iPerformedReprobes += 1;

/*
            if (isCurTest)
            {
                std::cout << "overflow add before " << sTestStr << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << " with reprobe " << iReprobe << " perf reprobe " << iPerformedReprobes << std::endl;
            }
*/

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe + iPerformedReprobes);

            TSX::tsx_key_t reprobePart = this->makeOverflowReprobe(iReprobe, iPerformedReprobes);

/*            if (isCurTest)
            {
                std::cout << "adding overflow kmer  " << sTestStr << " in " << iPos << " with reprobe " << iReprobe << " perf reprobe " << iPerformedReprobes << std::endl;
            }*/

            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t* pSavedKey = &savedkey;

            TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
            TSX::tsx_key_t updkey = TSX::tsx_key_t(2*m_iK, true, this->m_pPool);
            updkey = updkey | reprobePart; //UBigInt(iPerformedReprobes+iReprobe-iInitialReprobes, this->m_pPool);

            TSX::tsx_keyval_t oKeyVal(updkey, 2 * m_iK + m_iStorageBits);
            oKeyVal = (oKeyVal << m_iStorageBits);
            oKeyVal.resize(m_iKeyValBits);

            uint64_t iBitsToPos = (iPos)*m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE)*8);
            uint32_t iStartOffset = iBitsToPos-((sizeof(FIELDTYPE)*8) * iStartPos);
            FIELDTYPE* pPos = m_pCounterArray + iStartPos;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            asm volatile("":::"memory");


            SBIGINT::SBIGINT *pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT *pValue = m_pTMP_VALUE[iThreadID];
            SBIGINT::SBIGINT *pOrigValue = m_pTMP_ORIGINAL[iThreadID];
            SBIGINT::SBIGINT *pClear = m_pTMP_CLEAR[iThreadID];
            //this->performIncrement(&incRet);
            // increment element

            //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            SBIGINT::getFromMemory(pKeyVal, pOrigValue, iStartOffset, m_iKeyValBits, pPos);

            bool elemEmpty = true;
            for (uint8_t i = 0; i < pKeyVal->iFields; ++i) {
                pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_value_key.m_pArray[i];
                pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_key_value.m_pArray[i];

                if (pKeyVal->pdata[i])
                {
                    elemEmpty = false;
                }
            }

            //bool elemEmpty = SBIGINT::isZero(pKeyVal, m_iKeyValBits);

            if (elemEmpty) {
                for (uint8_t i = 0; i < oKeyVal.m_iFields; ++i) {
                    pKeyVal->pdata[i] = oKeyVal.m_pArray[i];
                }

                pKeyVal->pdata[0] = pKeyVal->pdata[0] | 1;
                // transaction completes here
                uint8_t CASresult = storeInMemoryCAS(pKeyVal, pOrigValue, pClear, pPos, iStartOffset, m_iKeyValBits, sizeof(FIELDTYPE) * 8, m_iKeyValBits);

                if (CASresult == 1) {
                    m_setUsedPositions.insert(iPos);

                    /*if (isCurTest)
                    {
                        std::cout << "overflow add empty after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;
                        std::cout << "overflow add empty " << sTestStr  << " " << iPos << " " << oKeyVal.to_string() << std::endl;
                        std::cout << "overflow add empty initial reprobes " << iInitialReprobes << " reprobes " << iReprobe << " perf reprobes " << iPerformedReprobes << std::endl;
                    }*/

                    //std::cout << "overflow add empty  after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;

                    return 1;

                } else if (CASresult == 0) {
                    ++iTotalAborts;
                    --iPerformedReprobes;
                    ++iAborts;



                    if (iTotalAborts % 10 == 0) {
                        std::cout << "aborts " << iTotalAborts << " inserts " << iAddCount << std::endl;
                    }

                    // try again
                    continue;

                }


            } else {


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
                bool bMatchesKey = positionMatchesOverflowReprobe(iPos, basekey, iReprobe, iPerformedReprobes); // was iReprobe+iPerformedReprobes

                if ((bIsKmerStart) || (!bMatchesKey)) {
                    //std::cout << "overflow Skipping position " << iPos << ": kmer=" << bIsKmerStart << " reprobes=" << iReprobe << " match=" << bMatchesKey << std::endl;

                    continue;
                }

                //std::cout << "overflow add before inc " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                uint8_t incremented = this->incrementElement_new(kmer, iPos, basekey, iReprobe+iPerformedReprobes, iInitialReprobes, true, verbose); // was iReprobe+iPerformedReprobes

                if (incremented == 0) {

                    // repeat -> hence undo reprobe
                    --iPerformedReprobes;

                    // redo increment ...
                    continue;
                } else if (incremented == 2){
                    this->handleOverflow(kmer, iPos, basekey, iReprobe + iPerformedReprobes, iInitialReprobes, verbose);
                }
                //std::cout << "overflow add after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                return 1;

            }
        }


        std::cout << "TSX OVERFLOW RET 0" << std::endl;

        return 0;
    }

protected:

    SBIGINT::SBIGINT** m_pTMP_ORIGINAL;
    SBIGINT::SBIGINT** m_pTMP_CLEAR;

    virtual void initialiseLocks()
    {
        TSXHashMapPerf::initialiseLocks();

        m_pTMP_ORIGINAL = (SBIGINT::SBIGINT**) malloc(sizeof(SBIGINT::SBIGINT*) * m_iThreads);
        m_pTMP_CLEAR = (SBIGINT::SBIGINT**) malloc(sizeof(SBIGINT::SBIGINT*) * m_iThreads);


        for (uint8_t i = 0; i < m_iThreads; ++i)
        {
            m_pTMP_ORIGINAL[i] = SBIGINT::getEmptySBIGINT(m_iKeyValBits + sizeof(FIELDTYPE)*8); // at max one field more
            m_pTMP_CLEAR[i] = SBIGINT::getEmptySBIGINT(m_iKeyValBits + sizeof(FIELDTYPE)*8); // at max one field more
        }


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



#endif //TSXCOUNT_TSXHASHMAPCAS_H
