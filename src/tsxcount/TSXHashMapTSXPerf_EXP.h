//
// Created by mjopp on 03/02/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPTSXPERFEXP_H
#define TSXCOUNT_TSXHASHMAPTSXPERFEXP_H


#include "TSXHashMapTSXPerf.h"
#include "commons.h"
#include <pthread.h>
#include <stack>
#include <map>

#include <xmmintrin.h>
#include <immintrin.h>
#include <unistd.h>
#include <cmath>
#include <src/tsxutils/SBigInt.h>

class TSXHashMapTSXPerfExperimental : public TSXHashMapTSXPerf {

public:

    TSXHashMapTSXPerfExperimental(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
    : TSXHashMapTSXPerf(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);

    }

    virtual ~TSXHashMapTSXPerfExperimental()
    {

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

        uint8_t iThreadID = omp_get_thread_num();
        uint8_t allowedTSXRetries = 10;

        //std::cout << "in inc key value" << std::endl;

        for (uint8_t iTSXRetries=0; iTSXRetries < allowedTSXRetries; ++iTSXRetries) {


            uint64_t iBitsToPos = (iPosition) * m_iKeyValBits;
            uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE) * 8);
            uint32_t iStartOffset = iBitsToPos - ((sizeof(FIELDTYPE) * 8) * iStartPos);
            FIELDTYPE *pPos = m_pCounterArray + iStartPos;

            // THIS PREFETCH is necessary to avoid stupid status==0...
            // TSX::tsx_val_t value = UBigInt(m_iStorageBits, true, this->m_pPool);

            SBIGINT::SBIGINT* pKeyVal = m_pTMP_KEYVAL[iThreadID];
            SBIGINT::SBIGINT* pValue = m_pTMP_VALUE[iThreadID];
            SBIGINT::clear(pKeyVal); // ADD 210621


            // THIS PREFETCH is necessary to avoid stupid status==0...
            volatile FIELDTYPE prefVal;
            uint8_t i = 0;
            for (i = 0; i < pKeyVal->iFields; ++i)
            {
                __builtin_prefetch(pPos+i, 1, 3);
                __atomic_compare_exchange_n(pPos+i, pPos+i, *(pPos+i), true, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
                prefVal = pPos[i];
                prefVal = pKeyVal->pdata[i];
            }
            __atomic_compare_exchange_n(pPos+i, pPos+i, *(pPos+i), true, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
            __builtin_prefetch(pPos+i, 1, 3);

            bool elemEmpty = false;

            uint status = _xbegin();
            if (status == _XBEGIN_STARTED) {

                //this->performIncrement(&incRet);
                // increment element

                //(*pSavedKey).copy_content_bits(pPos, iStartOffset, m_iKeyValBits);
                SBIGINT::getFromMemory(pKeyVal, iStartOffset, m_iKeyValBits, pPos);

                //value = (*pSavedKey) &  m_mask_value_key;

                for (i = 0; i < pKeyVal->iFields; ++i)
                {
                    // 11111111111111111111111111110000
                    pValue->pdata[i] = pKeyVal->pdata[i] & m_mask_value_key.m_pArray[i];
                    pKeyVal->pdata[i] = pKeyVal->pdata[i] & m_mask_key_value.m_pArray[i];
                }

                //00000000000000000000000000001110
                SBIGINT::bitComplement(pValue);
                //11111111111111111111111111110001
                elemEmpty = SBIGINT::isZero(pValue, m_iStorageBits);


                // check whether overflow occurs, or not
                // elemEmpty == true => OVERFLOW OCCURS!
                //elemEmpty = false;//(~value).isZero();

                if (elemEmpty) {

                    // 210618 not required b/c already 0 and not used thereafter!
                    //for (i = 0; i < pKeyVal->iFields; ++i)
                    //{
                    //    pValue->pdata[i] = 0;
                    //}

                } else {
                    SBIGINT::bitComplement(pValue);
                    //00000000000000000000000000001110
                    SBIGINT::add_simpleSBIGINT(pValue);
                    //00000000000000000000000000001111

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
                _xend();
                //asm volatile("":: :"memory");


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

            } else {

                ++iTotalAborts;
                pCounter[iThreadID][status].increment();

            }


        }

        return 0;

    }



};


#endif //TSXCOUNT_TSXHASHMAPTSXPERFEXP_H
