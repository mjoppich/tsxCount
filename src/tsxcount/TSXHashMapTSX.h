//
// Created by mjopp on 21/10/2018.
//

#ifndef TSXCOUNT_TSXHASHMAPTSX_H
#define TSXCOUNT_TSXHASHMAPTSX_H

#include "TSXHashMap.h"
#include <pthread.h>
#include <stack>

#include <immintrin.h>

class TSXHashMapTSX : public TSXHashMap {

public:

    TSXHashMapTSX(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
            : TSXHashMap(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    ~TSXHashMapTSX()
    {
        free(m_pLocked);
    }

    bool countKmer(TSX::tsx_kmer_t& kmer)
    {
        return this->addKmer(kmer);
    }

protected:

    void setThreads(uint8_t iThreads)
    {
        m_iThreads = iThreads;
    }

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
        _xbegin();
    }


    /**
     *
     * unlocks all locked array positions for given thread
     *
     * @param iThreadID
     */
    void unlock_thread(uint8_t iThreadID)
    {
        pthread_mutex_lock(&m_oLockMutex);
        m_pLocked[iThreadID].clear();
        pthread_mutex_unlock(&m_oLockMutex);

    }

    bool releaseLock(uint8_t iThreadID, uint64_t iPos)
    {

    }


private:


};


#endif //TSXCOUNT_TSXHASHMAPPTHREAD_H

