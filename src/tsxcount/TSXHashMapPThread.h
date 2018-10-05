//
// Created by mjoppich on 11/15/17.
//

#ifndef TSXCOUNT_TSXHASHMAPPTHREAD_H
#define TSXCOUNT_TSXHASHMAPPTHREAD_H

#include "TSXHashMap.h"
#include <pthread.h>
#include <stack>

class TSXHashMapPThread : public TSXHashMap {

public:

    TSXHashMapPThread(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
    : TSXHashMap(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    ~TSXHashMapPThread()
    {
        free(m_pLocked);
    }

protected:

    void setThreads(uint8_t iThreads)
    {
        m_iThreads = iThreads;
    }

    void initialiseLocks()
    {

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
                return i;
            }
        }

        return -1;
    }

    bool canAcquireLock(uint8_t iThreadID, uint64_t iArrayPos)
    {
        uint8_t iPosLocked = this->position_locked(iArrayPos);
        return iPosLocked == -1 or iThreadID == iPosLocked;
    }

    /**
     *
     * @param iThreadID thread id for which the lock is acquired
     * @param iArrayPos array position for which the lock is acquired
     * @return True if lock successfully acquired, False otherwise
     */
    bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos)
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
    void unlock_thread(uint8_t iThreadID)
    {
        pthread_mutex_lock(&m_oLockMutex);
        m_pLocked[iThreadID].clear();
        pthread_mutex_unlock(&m_oLockMutex);

    }

    bool releaseLock(uint8_t iThreadID, uint64_t iPos)
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

        return bRetVal;
    }



    bool addKmer(TSX::tsx_kmer_t& kmer)
    {
        bool bInserted = false;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        uint8_t iThreadID = omp_get_thread_num();

        while ( iReprobes < m_iMaxReprobes )
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            this->acquireLock( iThreadID, iPos );

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
                TSX::tsx_key_t updkey = this->makeKey(key, iReprobes);

                this->incrementElement(iPos, key, iReprobes, false);

                // so we can find kmer start positions later without knowing the kmer
                m_iKmerStarts.setBit(iPos, 1);
                // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                m_setUsedPositions.insert(iPos);

                bInserted = true;

                this->releaseLock(iThreadID, iPos);

                break;

            } else {

                // TODO where is the case kmer starts multiple positions later handled?
                bool bIsKmerStart = m_iKmerStarts.getBit(iPos) == 1;
                bool bMatchesKey;// = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                if (bIsKmerStart)
                {
                    bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
                } else {
                    bMatchesKey = false;
                }

                if ((bIsKmerStart) && (bMatchesKey))
                {

                    this->acquireLock(iThreadID, iPos);
                    bool bOverflow = this->incrementElement(iPos, basekey, iReprobes, false);

                    if (bOverflow)
                    {
                        handleOverflow(iPos, basekey, iReprobes);
                    }

                    this->releaseLock(iThreadID, iPos);

                    bInserted = true;
                    break;

                } else {

                    bool bMatchesKey2 = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                    ++iReprobes;
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



private:


};


#endif //TSXCOUNT_TSXHASHMAPPTHREAD_H
