//
// Created by mjopp on 21/10/2018.
//

#ifndef TSXCOUNT_TSXHASHMAPTSX_H
#define TSXCOUNT_TSXHASHMAPTSX_H

#include "TSXHashMap.h"
#include <pthread.h>
#include <stack>

#include <immintrin.h>
#include <unistd.h>

#define _XA_EXPLICIT		0
#define _XA_RETRY		1
#define _XA_CONFLICT		2
#define _XA_CAPACITY		3

#ifdef ABORT_COUNT
__thread unsigned _abrt[4] __attribute__((aligned(L1DSZ)));
#define ABRT_COUNT(type, status)					\
do {									\
	if (status & (1 << type))					\
		_abrt[type]++;						\
} while (0)
#else
#define ABRT_COUNT(...)
#endif

class TSXHashMapTSX : public TSXHashMap {

public:

    TSXHashMapTSX(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
            : TSXHashMap(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    virtual ~TSXHashMapTSX()
    {
        free(m_pLocked);
    }
    size_t iAddCount = 0;
    size_t iAddKmerCount = 0;
    size_t iAborts = 0;

    virtual bool addKmer(TSX::tsx_kmer_t& kmer)
    {
        bool bInserted = false;

        this->iAddKmerCount+=1;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        uint8_t iThreadID = omp_get_thread_num();
        uint8_t iRetries = 0;

        while ( iReprobes < m_iMaxReprobes )
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            // start transaction here
            //unsigned status = _xbegin();

            bool bEmpty = false;


            bEmpty = positionEmpty(iPos);

            if (bEmpty) {

                TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
                TSX::tsx_key_t updkey = this->makeKey(key, iReprobes);

                // no overflow can happen here ...
                CIncrementElement incRet = this->incrementElement(iPos, key, iReprobes, false);

                uint status = _xbegin();
                if(status == _XBEGIN_STARTED) {
                    this->performIncrement(&incRet);

                    // compare elements
                    bool elemsEqual = this->getElement(incRet.iPosition) == incRet.keyval;

                    if (!elemsEqual)
                    {
                        _xabort(0xff);
                    }

                    // transaction completes here
                    _xend();

                    // so we can find kmer start positions later without knowing the kmer
                    m_iKmerStarts.setBit(iPos, 1);
                    // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                    m_setUsedPositions.insert(iPos);

                    this->iAddCount += 1;
                    bInserted = true;
                    break;
                }




                if (status != _XBEGIN_STARTED) {
                    //std::cout << "STATUS ERROR " << (uint) status << std::endl;
                    if (_XABORT_CODE(status)==0xff)
                    {
                        ++iAborts;
                    }
                    // transaction aborted
                    // => retry
                    continue;
                }


            } else {

                // TODO where is the case kmer starts multiple positions later handled?
                bool bIsKmerStart = m_iKmerStarts.getBit(iPos) == 1;
                bool bMatchesKey;// = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                if (bIsKmerStart) {
                    bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
                } else {
                    bMatchesKey = false;
                }

                if ((bIsKmerStart) && (bMatchesKey)) {

                    CIncrementElement incRet = this->incrementElement(iPos, basekey, iReprobes, false);
                    std::vector<CIncrementElement> vOPS;
                    vOPS.insert(vOPS.end(), incRet);

                    if (incRet.iOverflow == 1) {

                        this->handleOverflow(iPos, basekey, iReprobes, &vOPS);

                    }

                    /*
                    if (vOPS.size() > 1)
                    {
                        std::cout << "OVERFLOW INCOMMING" << std::endl;


                        for (size_t i = 0; i < vOPS.size(); ++i) {

                            CIncrementElement* pINC = &(vOPS.data()[i]);
                            std::cout << pINC->keyval.to_string() << " overflow " << (int) pINC->iOverflow << " position " << pINC->iPosition << " " << this->getElement(pINC->iPosition).to_string() << std::endl;
                        }
                    }
                     */

                    int iAbort = 0;
                    uint status = _xbegin();

                    if(status == _XBEGIN_STARTED) {

                        for (size_t i = 0; i < vOPS.size(); ++i) {

                            CIncrementElement* pINC = &(vOPS.data()[i]);

                            bool elemsEqual = this->getElement(pINC->iPosition) == pINC->original;

                            if (!elemsEqual)
                            {
                                //iAbort = 1;
                                _xabort(0xff);
                            }

                            // we want to get the iPosition-th entry
                            std::div_t oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);

                            //UBigInt oReturn(m_iKeyValBits);
                            //UBigInt::storeIntoMemory(m_pCounterArray, (uint32_t) oStartPos.quot, (uint32_t) oStartPos.rem, pINC->keyval, m_iKeyValBits);
                            FIELDTYPE* pPos = m_pCounterArray + oStartPos.quot;
                            pINC->keyval.copy_content_to_array(pPos, oStartPos.rem, m_iKeyValBits);

                            // compare elements
                            elemsEqual = this->getElement(pINC->iPosition) == pINC->keyval;

                            if (!elemsEqual)
                            {
                                //iAbort = 1;
                                _xabort(0xff);
                            }
                        }

                        _xend();



                        this->iAddCount += 1;
                        bInserted = true;
                        break;

                    }

                    if (status != _XBEGIN_STARTED) {

                        if (_XABORT_CODE(status)==0xff)
                        {
                            ++iAborts;
                        }
                        // transaction aborted
                        // => retry

                        /*
                        std::cout << "STATUS ERROR NE " << (uint) status << " " << (int) _XABORT_CODE(status) << " " << vOPS.size() << std::endl;


                        for (size_t i = 0; i < vOPS.size(); ++i) {

                            CIncrementElement* pINC = &(vOPS.data()[i]);
                            std::cout << pINC->keyval.to_string() << " overflow " << (int) pINC->iOverflow << " position " << pINC->iPosition << " " << this->getElement(pINC->iPosition).to_string() << std::endl;
                        }
                         */

                        continue;
                    }



                } else {

                    // @TODO why was that there anyhow?
                    //bool bMatchesKey2 = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                    // transaction completed => incorrect position => reprobe
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

            //this->addKmer(kmer);
        }

        return bInserted;

    }

    bool countKmer(TSX::tsx_kmer_t& kmer)
    {
        return this->addKmer(kmer);
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


private:


};


#endif //TSXCOUNT_TSXHASHMAPPTHREAD_H

