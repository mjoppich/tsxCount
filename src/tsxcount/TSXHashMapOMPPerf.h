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
    virtual bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos)
    {

        omp_set_lock(&m_oOMPLock);
        bool success = false;

        if (this->canAcquireLock(iThreadID, iArrayPos)) {

            //std::cerr << "Locked Pos " << iArrayPos << " for thread " << (int) iThreadID << std::endl;


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
    virtual void unlock_thread(uint8_t iThreadID)
    {
        omp_set_lock(&m_oOMPLock);
        m_pLocked[iThreadID].clear();
        omp_unset_lock(&m_oOMPLock);

    }

    virtual bool releaseLock(uint8_t iThreadID, uint64_t iPos)
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


    virtual ~TSXHashMapOMPPerf()
    {
        free(m_pLocked);


        free(m_pTMP_KEYVAL);
        free(m_pTMP_VALUE);
    }
    size_t iAddCount = 0;
    size_t iAddKmerCount = 0;
    size_t iAborts = 0;
    size_t iTotalAborts = 0;


protected:


};


#endif //TSXCOUNT_TSXHASHMAPOMPPERF_H
