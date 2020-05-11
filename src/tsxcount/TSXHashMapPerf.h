//
// Created by mjopp on 06/03/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPPERF_H
#define TSXCOUNT_TSXHASHMAPPERF_H


#include "TSXHashMap.h"
#include <src/tsxutils/SBigInt.h>

#include <pthread.h>
#include <stack>

#include <stdbool.h>

#include <immintrin.h>
#include <unistd.h>
#include <atomic>
#include <cmath>
#include <bitset>

#ifdef __cplusplus
#include <atomic>
#else
#include <stdatomic.h>
#endif




class TSXHashMapPerf : public TSXHashMap {

public:

    TSXHashMapPerf(uint8_t iL, uint32_t iStorageBits, uint16_t iK)
            : TSXHashMap(iL, iStorageBits, iK) {

        this->setThreads(1);

    }


    virtual ~TSXHashMapPerf() {
        free(m_pLocked);
        free(m_pTMP_KEYVAL);
        free(m_pTMP_VALUE);
    }

    size_t iAddCount = 0;
    size_t iAddKmerCount = 0;
    size_t iAborts = 0;
    size_t iTotalAborts = 0;


    virtual bool addKmer(TSX::tsx_kmer_t& kmer, bool verbose=false) override
    {

        bool bInserted = false;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        uint8_t iThreadID = omp_get_thread_num();

        while ( iReprobes < m_iMaxReprobes )
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            bool lockAcquired = this->acquireLock(iThreadID, iPos);

            if (!lockAcquired)
            {
                continue;
            }

            TSX::tsx_keyval_t key_reprobe_shift = this->makeKey(basekey, iReprobes);
            key_reprobe_shift.resize(m_iKeyValBits);
            key_reprobe_shift.shiftLeft(m_iStorageBits);// (key_reprobe_shift << m_iStorageBits);,
            key_reprobe_shift.m_pArray[0] = key_reprobe_shift.m_pArray[0] | 1;

            uint64_t iBitsToPos = (iPos) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            SBIGINT::SBIGINT *pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT *pKeyReprobeShift = m_pTMP_VALUE[iThreadID];

            SBIGINT::fromClassToStruct(&key_reprobe_shift, pKeyReprobeShift);

            int iinc = 0;

            uint8_t i = 0;

            //volatile FIELDTYPE prefVal;
            //for (i = 0; i < pKeyVal->iFields; ++i) {
            //    __atomic_fetch_or (pPos+i, pPos[i], __ATOMIC_RELAXED);
            //    __atomic_store(pPos + i, pPos + i, __ATOMIC_RELAXED);
            //    prefVal = pKeyVal->pdata[i];
            //    prefVal = pKeyReprobeShift->pdata[i];
            //}
            //prefVal = pPos[pKeyVal->iFields];
            //__atomic_store(pPos + i, pPos + i, __ATOMIC_RELAXED);


            uint16_t iKeyValBits = m_iKeyValBits;
            asm volatile("":: :"memory");


            // increment element
            SBIGINT::getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);

            bool elemEmpty = SBIGINT::isZero(pKeyVal, m_iKeyValBits);


            if (elemEmpty) {
                //SBIGINT::getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);
                SBIGINT::storeInMemory(pKeyReprobeShift->pdata, pPos, iStartOffset, m_iKeyValBits, sizeof(FIELDTYPE) * 8);

                // transaction completes here
                asm volatile("":: :"memory");

                //std::cout << "after add and transaction" << std::endl;

                // so we can find kmer start positions later without knowing the kmer
                m_iKmerStarts.setBit(iPos, 1);
                // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                m_setUsedPositions.insert(iPos);

                this->releaseLock(iThreadID, iPos);
                bInserted = true;
                break;

            } else {

                // TODO where is the case kmer starts multiple positions later handled?

                //std::cout << "managed abort" << std::endl;
                // TODO where is the case kmer starts multiple positions later handled?
                bool bIsKmerStart = m_iKmerStarts.getBit(iPos) == 1;
                bool bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);



                if ((bIsKmerStart) && (bMatchesKey))
                {
                    uint8_t iAddStatus = 0;
                    uint8_t iOverflowStatus = 0;

                    while (iAddStatus == 0) {
                        iAddStatus = this->incrementElement_new(kmer, iPos, basekey, iReprobes, iReprobes, false,
                                                                verbose);
                    }

                    this->releaseLock(iThreadID, iPos);

                    if (iAddStatus == 2) {
                        /*
                        if (isCurTest) {
                            std::cout << sTestStr << " before overflow" << std::endl;
                        }
                         */

                        while (iOverflowStatus == 0) {
                            iOverflowStatus = this->handleOverflow(kmer, iPos, basekey, iReprobes, iReprobes,
                                                                   verbose);
                        }
                        /*
                        if (isCurTest) {
                            std::cout << sTestStr << " after overflow increment " << (int) iAddStatus << " "
                                      << (int) iOverflowStatus << std::endl;
                        }
                         */
                    }

                    bInserted = true;

                    break;

                    /*
                    CIncrementElement incRet = this->incrementElement(iPos, basekey, iReprobes, 1, false);
                    std::vector<CIncrementElement> vOPS;
                    vOPS.insert(vOPS.end(), incRet);

                    if (incRet.iOverflow == 1)
                    {
                        uint8_t iOvflw = this->handleOverflow(iPos, basekey, iReprobes, &vOPS, &kmer);

                        if (iOvflw == 2)
                        {
                            // TRY AGAIN
                            continue;
                        }

                        //this->performIncrement(&incRet);
                    }

                    std::vector<size_t> usedPos;
                    for (auto rit = vOPS.rbegin(); rit != vOPS.rend(); ++rit)
                    {
                        CIncrementElement inc = (*rit);

                        //std::cout << inc.iPosition;
                        //inc.keyval.print_string();
                        this->performIncrement(&inc);
                        this->releaseLock(iThreadID, inc.iPosition);

                    }
                     */

                    //this->releaseLock(iThreadID, iPos);


                } else {

                    //bool bMatchesKey2 = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
                    ++iReprobes;

                    this->releaseLock(iThreadID, iPos);

                    continue;
                }


            }


        }

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

            this->addKmer(kmer);
        }

        return bInserted;

    }

    volatile FIELDTYPE TSXHASHMAPPREFETCH;

    /**
     *
     * @param kmer
     * @param iPosition
     * @param key
     * @param reprobes
     * @param verbose
     * @return
     */
    virtual uint8_t incrementElement_key_value(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false)
    {

        /**
         *
         * ONLY INCREMENT VALUE PART => ret 2 on overflow
         *
         */

        std::string stest = m_mask_value_key.to_string();
        stest = m_mask_key_value.to_string();

        uint8_t iThreadID = omp_get_thread_num();


        //std::cout << "in inc key value" << std::endl;

        for (uint8_t iTSXRetries=0; iTSXRetries < 10; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            TSX::tsx_val_t value = UBigInt(m_iStorageBits, true, this->m_pPool);

            SBIGINT::SBIGINT* pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT* pValue = m_pTMP_VALUE[iThreadID];
            uint8_t i;
            bool elemEmpty = false;

            asm volatile("":: :"memory");



            //this->performIncrement(&incRet);
            // increment element

            //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
            SBIGINT::getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);

            //value = (*pSavedKey) &  m_mask_value_key;

            for (i = 0; i < pKeyVal->iFields; ++i)
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

                for (i = 0; i < pKeyVal->iFields; ++i)
                {
                    pValue->pdata[i] = 0;
                }

            } else {
                bitComplement(pValue);

                add_simpleSBIGINT(pValue);

                for (i = 0; i < pKeyVal->iFields; ++i)
                {
                    // keeps key (func+reprobe), nulls value
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] | pValue->pdata[i];
                }
            }

            //pSavedKey->bitAnd(m_mask_key_value);
            //pSavedKey->bitOr(value);

            //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
            //(*pSavedKey).copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);
            SBIGINT::storeInMemory(pKeyVal->pdata, pPos, iStartOffset, m_iKeyValBits, pKeyVal->iFieldSize);

            // transaction completes here

            //UBigInt value = fromStructToClass(pValue, m_pPool);
            //UBigInt saved = fromStructToClass(pKeyVal, m_pPool);

            //std::cout << "after inc position " << elemEmpty << " " << iPosition << " " << value.to_string() << " " << saved.to_string() << std::endl;

            if (elemEmpty)
            {
                //std::cout << "elemEmpty" << std::endl;
                // OVERFLOW in VALUE part
                return 2;
            }

            // NO OVERFLOW
            return 1;

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
    virtual uint8_t incrementElement_func(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool verbose=false)
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

            // THIS PREFETCH is necessary to avoid stupid status==0...
            TSXHASHMAPPREFETCH = pPos[0];
            asm volatile("":: :"memory");

            SBIGINT::SBIGINT* pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT* pValue = m_pTMP_VALUE[iThreadID];

            m_mask_kv_reprobevalue;

            bool elemEmpty = false;


            //this->performIncrement(&incRet);
            // increment element

            SBIGINT::getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);
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
            SBIGINT::storeInMemory(pKeyVal->pdata, pPos, iStartOffset, m_iKeyValBits, pKeyVal->iFieldSize);


            // transaction completes here


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


        }

        return 0;

    }

    TSX::tsx_key_t makeOverflowReprobe(uint32_t iReprobe, uint32_t iPerfReprobe)
    {

        uint32_t iReprobeLength = m_iL;
        uint32_t iPerfReprobeBits = (uint32_t) std::floor(iReprobeLength/2.0);

        TSX::tsx_key_t oRet = UBigInt(iReprobeLength, true, m_pPool);
        oRet.resize(iReprobeLength);
        TSX::tsx_key_t oReprobe = UBigInt(iReprobe, m_pPool);
        oReprobe.resize(iReprobeLength);
        TSX::tsx_key_t oPerfReprobe = UBigInt(iPerfReprobe, m_pPool);
        oPerfReprobe.resize(iReprobeLength);

        oRet = oRet | oReprobe;
        oRet = oRet << iReprobeLength-iPerfReprobeBits;
        oRet = oRet | oPerfReprobe;

        return oRet;

    }

    bool positionMatchesOverflowReprobe(uint64_t pos, TSX::tsx_key_t& key, uint32_t iReprobe, uint32_t iPerfReprobe)
    {

        // get key from position
        TSX::tsx_keyval_t oKeyVal = getElement(pos);
        TSX::tsx_key_t oKey = this->getKeyFromKeyVal(oKeyVal);

        TSX::tsx_key_t oReprobe = this->makeOverflowReprobe(iReprobe, iPerfReprobe);//iReprobe, this->m_pPool);
        oReprobe.resize(2*m_iK);

        // extract position reprobe and compare with reprobe value
        bool reprobeMatch = (oKey & m_mask_reprobe_func) == oReprobe;

        return reprobeMatch;

    }

    void print_stats()
    {
        std::cerr << "Used fields: " << m_setUsedPositions.size() << std::endl;
        std::cerr << "Available fields: " << std::pow(2.0, m_iL) << std::endl;
        std::cerr << "k=" << m_iK << " l=" << (uint32_t) m_iL << " entry (key+value) bits=" << m_iKeyValBits << " storage bits=" << m_iStorageBits << std::endl;
    }

    std::vector<uint64_t> getKmerPositions(TSX::tsx_kmer_t& kmer)
    {
        // TODO this is a copy of getKmerCount and should be removed
        bool bFound = false;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        UBigInt oResult(64, true, this->m_pPool);

        std::vector<uint64_t> vPositions;

        while ((!bFound) && ( iReprobes < m_iMaxReprobes))
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (!bEmpty)
            {
                TSX::tsx_keyval_t elem = this->getElement(iPos);
                bool bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                if (bMatchesKey)
                {
                    bFound = true;


                    oResult = this->getValFromKeyVal(elem);

                    // TODO find possible remaining entries!
                    UBigInt oOverflows = findOverflowCounts(iPos, basekey, iReprobes, false);

                    vPositions.push_back(iPos);


                    uint32_t iUsedOverflowBits = oOverflows.getBitCount();
                    if (iUsedOverflowBits > 0)
                    {
                        oOverflows.resize( iUsedOverflowBits + m_iStorageBits );

                        oOverflows = oOverflows << m_iStorageBits;
                        oOverflows = oOverflows | oResult;

                        oResult = oOverflows;
                    }

                }

            } else {
                // why would the insertion skip an empty place?
                break;
            }
            ++iReprobes;
        }

        if (!bFound)
        {
            //std::cerr << "Kmer " << kmer.to_string() << " not in hash" << std::endl;
        }
        return vPositions;
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
    uint8_t incrementElement_new(TSX::tsx_kmer_t& kmer, uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, uint32_t iInitialReprobes, bool bKeyIsValue, bool verbose=false)
    {

        //std::string sTestStr = "TTCCATTCCATTCC";
        //UBigInt sTest = fromSequence(sTestStr, m_pPool);

        //bool isCurTest = kmer == sTest;
        //bool isCurTest = false;
        bool handleFuncOverflow = false;

        //std::cout << "increment element tsx " << iPosition << std::endl;

        // try to increment value
        uint8_t iIncrementTries = 0;
        uint8_t iIncrementState = 0;

        bool singleVerbose = false;

        //for (iIncrementTries=0; iIncrementTries < 100; ++iIncrementTries)
        while(iIncrementState == 0)
        {
            iIncrementState = this->incrementElement_key_value(kmer, iPosition, key, reprobes, verbose | singleVerbose);

            if (iIncrementState == 0)
            {
                //std::cout << "simple increment state 0" << std::endl;
                singleVerbose = true;
                continue;

            } else if (iIncrementState == 1)
            {

                /*
                if (isCurTest)
                {
                    std::cout << "return 1 after increment, pos="<< iPosition << std::endl;
                }
                 */

                // value increment, everything done
                return 1;
            } else if (iIncrementState == 2)
            {
                // value incremented, but overflow occurred
                if (bKeyIsValue)
                {

                    /*
                    if (isCurTest)
                    {
                        std::cout << "increment returns 2, break, pos="<< iPosition << std::endl;
                    }
                    */

                    handleFuncOverflow = true;
                    break;

                } else {

                    /*
                    if (isCurTest)
                    {
                        std::cout << "increment returns 2, no break, pos="<< iPosition << std::endl;
                    }
                    */

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

            /*
            if (isCurTest)
            {
                std::cout << "cur test overflow in increment" << " pos="<< iPosition<< std::endl;
            }
             */
            // try to avoid overflow by incrementing func part
            //for (iFuncIncrementTries=0; iFuncIncrementTries < 10; ++iFuncIncrementTries) {
            while(true) {

                iIncrementFuncState = this->incrementElement_func(kmer, iPosition, key, reprobes, verbose);

                /*
                if (isCurTest)
                {
                    std::cout << "cur test after overflow in func " << (int) iIncrementFuncState << " pos="<< iPosition<< std::endl;
                }
                 */

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

                    uint8_t iOverflowReturn = this->handleOverflow(kmer, iPosition, key, reprobes, iInitialReprobes, verbose);

                    if (iOverflowReturn == 0)
                    {
                        std::cout << "overflow return 0" << " pos="<< iPosition<< std::endl;
                    }

                    return 1;

                }


            }


        } else {
            /*
            if (isCurTest)
            {
                std::cout << "cur test no overflow in increment" << std::endl;
            }
             */
        }


        std::cout << handleFuncOverflow << " " << bKeyIsValue << " " << (int) iIncrementTries << " " << (int) iIncrementState << " " << (int)iFuncIncrementTries << " " << (int)iIncrementFuncState << std::endl;
        this->print_stats();

        exit(33);

        // should never reach this => not inserted
        return 0;
    }

    virtual uint8_t handleOverflow(TSX::tsx_kmer_t& kmer, uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, uint32_t iInitialReprobes, bool verbose=false)
    {
        uint32_t iPerformedReprobes = 0;
        uint8_t iThreadID = omp_get_thread_num();

        /*
        UBigInt uBigOne = UBigInt(1, m_pPool);
        std::string sTestStr = "AACAGCGTTTCGTC";
        UBigInt sTest = fromSequence(sTestStr, m_pPool);
        bool isCurTest = sTest == kmer;
        isCurTest = false;
        sTestStr = toSequence(kmer);

        if (isCurTest)
        {
            std::cout << "overflow add before " << sTestStr << " " << kmer.to_string() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << " with reprobe " << iReprobe << std::endl;
        }
         */

        while (iPerformedReprobes < m_iMaxReprobes) {

            iPerformedReprobes += 1;

            /*
            if (isCurTest) {
                std::cout << "overflow add before " << sTestStr << " " << this->getKmerCount(kmer).toUInt()
                          << " in pos " << iPos << " with reprobe " << iReprobe << " perf reprobe "
                          << iPerformedReprobes << std::endl;
            }
            */

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe + iPerformedReprobes);

            bool lockAcquired = this->acquireLock(iThreadID, iPos);

            if (!lockAcquired)
            {
                return 2;
            }

            TSX::tsx_key_t reprobePart = this->makeOverflowReprobe(iReprobe, iPerformedReprobes);

            /*
            if (isCurTest) {
                std::cout << "adding overflow kmer  " << sTestStr << " in " << iPos << " with reprobe " << iReprobe
                          << " perf reprobe " << iPerformedReprobes << std::endl;
            }
            */

            TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
            TSX::tsx_keyval_t *pSavedKey = &savedkey;

            TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
            TSX::tsx_key_t updkey = TSX::tsx_key_t(2 * m_iK, true, this->m_pPool);
            updkey = updkey | reprobePart; //UBigInt(iPerformedReprobes+iReprobe-iInitialReprobes, this->m_pPool);

            TSX::tsx_keyval_t oKeyVal(updkey, 2 * m_iK + m_iStorageBits);
            oKeyVal = (oKeyVal << m_iStorageBits);
            oKeyVal.resize(m_iKeyValBits);

            uint64_t iBitsToPos = (iPos) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;


            // THIS PREFETCH is necessary to avoid stupid status==0...
            //TSXHASHMAPPREFETCH = pPos[0];
            asm volatile("":: :"memory");

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

            if (elemEmpty) {


                oKeyVal.m_pArray[0] = oKeyVal.m_pArray[0] | 1;

                //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                oKeyVal.copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);

                // transaction completes here
                asm volatile("":: :"memory");

                m_setUsedPositions.insert(iPos);

                this->releaseLock(iThreadID, iPos);

                /*
                if (isCurTest) {
                    std::cout << "overflow add empty after " << kmer.to_debug() << " "
                              << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;
                    std::cout << "overflow add empty " << sTestStr << " " << iPos << " " << oKeyVal.to_string()
                              << std::endl;
                    std::cout << "overflow add empty initial reprobes " << iInitialReprobes << " reprobes " << iReprobe
                              << " perf reprobes " << iPerformedReprobes << std::endl;

                }
                 */

                //std::cout << "overflow add empty  after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " " << iPos << std::endl;

                return 1;

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
                    this->releaseLock(iThreadID, iPos);

                    continue;
                }

                //std::cout << "overflow add before inc " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                uint8_t incremented = this->incrementElement_new(kmer, iPos, basekey, iReprobe + iPerformedReprobes,
                                                                 iInitialReprobes, true,
                                                                 verbose); // was iReprobe+iPerformedReprobes

                this->releaseLock(iThreadID, iPos);

                if (incremented == 0) {

                    // repeat -> hence undo reprobe
                    --iPerformedReprobes;

                    // redo increment ...
                    continue;
                } else if (incremented == 2) {
                    this->handleOverflow(kmer, iPos, basekey, iReprobe + iPerformedReprobes, iInitialReprobes,
                                         verbose);
                }
                //std::cout << "overflow add after " << kmer.to_debug() << " " << this->getKmerCount(kmer).toUInt() << " in pos " << iPos << std::endl;

                return 1;


            }
        }

        std::cout << "TSX OVERFLOW RET 0" << std::endl;

        return 0;
    }

protected:

    SBIGINT::SBIGINT** m_pTMP_KEYVAL;
    SBIGINT::SBIGINT** m_pTMP_VALUE;

    virtual void initialiseLocks()
    {

        m_pTMP_KEYVAL = (SBIGINT::SBIGINT**) malloc(sizeof(SBIGINT::SBIGINT*) * m_iThreads);
        m_pTMP_VALUE = (SBIGINT::SBIGINT**) malloc(sizeof(SBIGINT::SBIGINT*) * m_iThreads);


        for (uint8_t i = 0; i < m_iThreads; ++i)
        {
            m_pTMP_KEYVAL[i] = SBIGINT::getEmptySBIGINT(m_iKeyValBits);
            m_pTMP_VALUE[i] = SBIGINT::getEmptySBIGINT(m_iKeyValBits);
        }

        std::cout << "keyval and value initialized" << std::endl;

    }

};


#endif //TSXCOUNT_TSXHASHMAPPERF_H
