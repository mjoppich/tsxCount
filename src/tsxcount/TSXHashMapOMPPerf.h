//
// Created by mjopp on 05/03/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPOMPPERF_H
#define TSXCOUNT_TSXHASHMAPOMPPERF_H
#include "TSXHashMap.h"
#include "TSXHashMapOMP.h"
#include "TSXHashMapPerf.h"
#include <pthread.h>
#include <stack>

#include <stdbool.h>
#include <omp.h>


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




class TSXHashMapOMPPerf : public TSXHashMapPerf {

public:

    TSXHashMapOMPPerf(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
            : TSXHashMapPerf(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);

    }

    void initialiseLocks() override
    {
        TSXHashMapPerf::initialiseLocks();

        m_pLocked = (std::vector<uint64_t>*) calloc(m_iThreads, sizeof(std::vector<uint64_t>));

        for (uint8_t i = 0; i < m_iThreads; ++i)
        {
            m_pLocked[i] = std::vector<uint64_t>();
        }


        omp_lock_hint_t lockHint = omp_lock_hint_speculative;
        omp_init_lock_with_hint(&m_oOMPLock, lockHint );
        //omp_init_lock(&m_oOMPLock);

    }

    /**
     *
     * @param iThreadID thread id for which the lock is acquired
     * @param iArrayPos array position for which the lock is acquired
     * @return True if lock successfully acquired, False otherwise
     */
    bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos) override
    {

        omp_set_lock(&m_oOMPLock);
        bool success = false;

        if (this->canAcquireLock(iThreadID, iArrayPos)) {
            m_pLocked[iThreadID].insert(m_pLocked[iThreadID].end(), iArrayPos);
            success = true;
        }

        omp_unset_lock(&m_oOMPLock);
        return success;

    }


    /**
     *
     * unlocks all locked array positions for given thread
     *
     * @param iThreadID
     */
    void unlock_thread(uint8_t iThreadID) override
    {
        omp_set_lock(&m_oOMPLock);
        m_pLocked[iThreadID].clear();
        omp_unset_lock(&m_oOMPLock);

    }

    bool releaseLock(uint8_t iThreadID, uint64_t iPos) override
    {
        omp_set_lock(&m_oOMPLock);
        std::vector<uint64_t>::iterator oIt = std::find(m_pLocked[iThreadID].begin(), m_pLocked[iThreadID].end(), iPos);

        bool bRetVal = false;

        if (oIt != m_pLocked[iThreadID].end())
        {
            bRetVal = true;

            m_pLocked[iThreadID].erase(oIt);
        }

        omp_unset_lock(&m_oOMPLock);

        return bRetVal;
    }





protected:


};


#endif //TSXCOUNT_TSXHASHMAPOMPPERF_H
