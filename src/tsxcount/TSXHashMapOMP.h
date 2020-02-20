//
// Created by mjopp on 22/10/2018.
//

#ifndef TSXCOUNT_TSXHASHMAPOMP_H
#define TSXCOUNT_TSXHASHMAPOMP_H


#include <omp.h>
#include "TSXHashMapPThread.h"

class TSXHashMapOMP : public TSXHashMapPThread{

public:

    TSXHashMapOMP(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
    : TSXHashMapPThread(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    virtual UBigInt findOverflowCounts(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, bool verbose)
    {
        bool bHandled = false;
        uint32_t iPerformedReprobes = 0;

        UBigInt oReturn(0, m_pPool);
        uint32_t iRequiredBits = 0;

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
                bool bMatchesKey = positionMatchesReprobe(iPos, basekey, iPerformedReprobes); // was iReprobe

                if (!bMatchesKey)
                {

                    if (verbose)
                    {
                        TSX::tsx_keyval_t elem = this->getElement(iPos);

                        std::cout << "OVFL Unmatched Key " << iPos << " " << iReprobe << " " << iPerformedReprobes << " " << iPerformedReprobes+iReprobe << std::endl;
                        std::cout << elem.to_string() << std::endl;
                    }

                    continue;
                }



                if (verbose)
                {
                    std::cout << "OVFL A Match Key " << iPos << " " << iReprobe << " " << iPerformedReprobes << " " << iPerformedReprobes+iReprobe << std::endl;
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


    virtual void initialiseLocks()
    {

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

};


#endif //TSXCOUNT_TSXHASHMAPOMP_H
