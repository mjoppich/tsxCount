//
// Created by mjopp on 15/03/2020.
//

#ifndef TSXCOUNT_TSXHASHMAPPTHREADPERF_H
#define TSXCOUNT_TSXHASHMAPPTHREADPERF_H

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




class TSXHashMapPThreadPerf : public TSXHashMapPerf {

public:

    TSXHashMapPThreadPerf(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
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

        pthread_mutex_init(&m_oLockMutex, NULL);

    }


    /**
     *
     * @param iArrayPos checks whether iArrayPos is locked
     * @return threadID of thread who locks iArrayPos or -1
     */
    uint8_t position_locked(uint64_t iArrayPos)
    {

        for (uint8_t i = 0; i < m_iThreads; ++i)
        {

            std::vector<uint64_t>::iterator iPos = std::find(m_pLocked[i].begin(), m_pLocked[i].end(), iArrayPos);

            if (iPos != m_pLocked[i].end())
            {
                return i+1;
            }
        }

        return 0;
    }

    bool canAcquireLock(uint8_t iThreadID, uint64_t iArrayPos) override
    {
        uint8_t iPosLocked = this->position_locked(iArrayPos);
        return iPosLocked == 0 or iThreadID+1 == iPosLocked;
    }

    /**
     *
     * @param iThreadID thread id for which the lock is acquired
     * @param iArrayPos array position for which the lock is acquired
     * @return True if lock successfully acquired, False otherwise
     */
    bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos) override
    {

        pthread_mutex_lock(&m_oLockMutex);
        bool success = false;

        if (this->canAcquireLock(iThreadID, iArrayPos)) {
            m_pLocked[iThreadID].insert(m_pLocked[iThreadID].end(), iArrayPos);
            success = true;
        }

        pthread_mutex_unlock(&m_oLockMutex);
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
        pthread_mutex_lock(&m_oLockMutex);
        m_pLocked[iThreadID].clear();
        pthread_mutex_unlock(&m_oLockMutex);

    }

    bool releaseLock(uint8_t iThreadID, uint64_t iPos) override
    {
        pthread_mutex_lock(&m_oLockMutex);
        std::vector<uint64_t>::iterator oIt = std::find(m_pLocked[iThreadID].begin(), m_pLocked[iThreadID].end(), iPos);

        bool bRetVal = false;

        if (oIt != m_pLocked[iThreadID].end())
        {
            bRetVal = true;

            m_pLocked[iThreadID].erase(oIt);
        }

        pthread_mutex_unlock(&m_oLockMutex);

        if (bRetVal == false)
        {
            std::cerr << "error releasing lock" << std::endl;
        }

        return bRetVal;
    }

protected:


};



#endif //TSXCOUNT_TSXHASHMAPPTHREADPERF_H
