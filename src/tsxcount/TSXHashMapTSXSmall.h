//
// Created by mjopp on 16/12/2018.
//

#ifndef TSXCOUNT_TSXHASHMAPTSXSMALL_H
#define TSXCOUNT_TSXHASHMAPTSXSMALL_H

#include <assert.h>
#include "TSXHashMapTSX.h"

class TSXHashMapTSXSmall : public TSXHashMapTSX {

    TSXHashMapTSXSmall(uint8_t iL, uint32_t iStorageBits, uint16_t iK, uint8_t iThreads=2)
            : TSXHashMapTSX(iL, iStorageBits, iK)
    {

        this->setThreads(iThreads);
        this->initialiseLocks();

    }

    virtual ~TSXHashMapTSXSmall()
    {
        free(m_pLocked);
    }
    size_t iAddCount = 0;
    size_t iAddKmerCount = 0;
    size_t iAborts = 0;

    bool incrementPosition(CIncrementElement* pINC)
    {

        TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
        TSX::tsx_keyval_t* pSavedKey = &savedkey;

        int iAbort = 0;
        uint status = _xbegin();

        if(status == _XBEGIN_STARTED) {


            std::div_t oStartPos;
            FIELDTYPE* pPos;
            bool elemsEqual;


            // calc pos
            uint32_t iBitsPerField = sizeof(FIELDTYPE) * 8;
            uint64_t iDivPos = pINC->iPosition * m_iKeyValBits;
            oStartPos = udiv( iDivPos, iBitsPerField);
            pPos = m_pCounterArray + oStartPos.quot;

            // get element
            pSavedKey->copy_content_bits(pPos, (uint32_t) oStartPos.rem, m_iKeyValBits);

            // compare with original value
            elemsEqual = pSavedKey->isEqual(pSavedKey, &(pINC->original), m_iKeyValBits);

            if (!elemsEqual)
            {
                //iAbort = 1;
                _xabort(0xff);
            }
            // increment element
            // we want to get the iPosition-th entry
            pINC->keyval.copy_content_to_array(pPos, oStartPos.rem, m_iKeyValBits);

            _xend();

            return true;

        }

        if (status != _XBEGIN_STARTED) {

            if (_XABORT_CODE(status)==0xff)
            {
                ++iAborts;
            } else {
                //std::cout << "unknown abort filled " << (int) status << std::endl;
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

            return false;
        }

        return false;
    }

    /**
     *
     * @param iPosition position in array
     * @param key key to look for
     * @param reprobes amount reprobes used
     * @param bKeyIsValue if true, the key part belongs to value
     * @return 1 if overflow, 0 if no overflow, 2 if locking error
     */
    CIncrementElement incrementElement(uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool bKeyIsValue)
    {

        TSX::tsx_keyval_t keyval = getElement(iPosition);
        TSX::tsx_val_t value = this->getValFromKeyVal(keyval);


        // while it did not work ....
        while (true)
        {
            CIncrementElement oRet;

            oRet.original = keyval;

            bool isZero = (~value).isZero();

            if (isZero) // value == 1111111
            {
                // WHY IS THIS VALUE == 1 ?
                value = 0;
                value.resize(m_iStorageBits);

                // maybe we can prevent the overflow?
                if (bKeyIsValue)
                {

                    TSX::tsx_func_t funcpart = this->getFuncFromKeyVal(keyval);

                    if ((~funcpart).isZero())
                    {
                        // overflow in func key part => set funcpart = 0 and propagate overflow further
                        funcpart = 0;
                        TSX::tsx_key_t updkey = (UBigInt(funcpart, 2*m_iK) << m_iL) | UBigInt(reprobes, this->m_pPool);

                        //storeElement(iPosition, updkey, value);

                        TSX::tsx_keyval_t oKeyVal(updkey, 2*m_iK + m_iStorageBits);
                        oKeyVal = (oKeyVal << m_iStorageBits) | value;
                        oKeyVal.resize(m_iKeyValBits);

                        oRet.iOverflow = 1;
                        oRet.iPosition = iPosition;
                        oRet.keyval = oKeyVal;

                        bool posIncrement = this->incrementPosition(&oRet);

                        if (posIncrement)
                        {
                            this->handleOverflow(iPosition, key, reprobes, nullptr);
                            return oRet;

                        }


                        //retry
                        continue;
                        //return 1;

                    } else {

                        ++funcpart;

                        TSX::tsx_key_t updkey = (UBigInt(funcpart, 2*m_iK) << m_iL) | UBigInt(reprobes, this->m_pPool);
                        value = 0;
                        value.resize(m_iStorageBits);

                        //storeElement(iPosition, updkey, value);
                        TSX::tsx_keyval_t oKeyVal(updkey, 2*m_iK + m_iStorageBits);
                        oKeyVal = (oKeyVal << m_iStorageBits) | value;
                        oKeyVal.resize(m_iKeyValBits);

                        oRet.iOverflow = 0;
                        oRet.iPosition = iPosition;
                        oRet.keyval = oKeyVal;

                        bool posIncrement = this->incrementPosition(&oRet);

                        if (posIncrement)
                            return oRet;

                        continue;
                    }

                }

                // also set value == 0
                TSX::tsx_key_t updkey = (key & m_mask_func_reprobe) | UBigInt(reprobes, this->m_pPool);
                value = 0;
                value.resize(m_iStorageBits);

                //storeElement(iPosition, updkey, value);

                // overflow occurred!
                TSX::tsx_keyval_t oKeyVal(updkey, 2*m_iK + m_iStorageBits);
                oKeyVal = (oKeyVal << m_iStorageBits) | value;
                oKeyVal.resize(m_iKeyValBits);

                oRet.iOverflow = 1;
                oRet.iPosition = iPosition;
                oRet.keyval = oKeyVal;

                bool posIncrement = this->incrementPosition(&oRet);

                if (posIncrement)
                {
                    this->handleOverflow(iPosition, key, reprobes, nullptr);

                    return oRet;
                }

                //retry
                continue;

            } else {

                bool bEmpty = this->positionEmpty( iPosition );

                ++value;

                TSX::tsx_key_t updkey;

                if (bKeyIsValue)
                {

                    if (bEmpty)
                    {

                        updkey = TSX::tsx_key_t(2*m_iK, true, this->m_pPool);
                        updkey = updkey | UBigInt(reprobes, this->m_pPool);
                        //storeElement(iPosition, updkey, value);

                    } else {

                        updkey = this->getKeyFromKeyVal(keyval);

                        //storeElement(iPosition, oldkey, value);
                    }


                } else {
                    updkey = this->makeKey(key, reprobes);
                    //storeElement(iPosition, updkey, value);
                }

                //std::cerr << value.to_string() << std::endl;

                TSX::tsx_keyval_t oKeyVal(updkey, 2*m_iK + m_iStorageBits);
                oKeyVal = (oKeyVal << m_iStorageBits) | value;
                oKeyVal.resize(m_iKeyValBits);

                oRet.iOverflow = 0;
                oRet.iPosition = iPosition;
                oRet.keyval = oKeyVal;

                //std::cerr << oRet.oval.to_string() << std::endl;

                bool posIncrement = this->incrementPosition(&oRet);

                if (posIncrement)
                    return oRet;

                //retry
                continue;

            }
        }

        //
        assert (false);
    }


    /**
  *
  * @param iPos
  * @param basekey
  * @param iReprobe
  * @return 1 if succeeded, 0 if not succeeded, 2 if blocked
  */
    virtual uint8_t handleOverflow(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, UBigInt* kmer=NULL)
    {
        uint32_t iPerformedReprobes = 0;
        uint8_t iThreadID = omp_get_thread_num();


        while (iPerformedReprobes < m_iMaxReprobes)
        {

            iPerformedReprobes += 1;

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe+iPerformedReprobes);

            bool lockAcquired = this->acquireLock(iThreadID, iPos);

            if (!lockAcquired)
            {
                return 2;
            }


            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {
                //std::cerr << "New field: POS " << std::to_string(iPos) << " for basekey " << basekey.to_string() << " with reprobes " << std::to_string(iPerformedReprobes) << std::endl;

                // the reprobe part must match the number of reprobes back to the previous entry!
                CIncrementElement oRet = this->incrementElement(iPos, basekey, iReprobe+iPerformedReprobes, true);

                continue;

                //this->releaseLock(iThreadID, iPos);
                /*
                if (kmer != NULL)
                {
                    std::cout << "Kmer " << kmer->to_string() << std::endl;
                    std::cout << "BaseKey " << basekey.to_string() << " " << iReprobe+iPerformedReprobes << " " << this->getPosition(basekey, iReprobe+iPerformedReprobes) << std::endl;
                    std::cout << "New Position " << iPos << " reprobe " << iReprobe << " iperfre " << iPerformedReprobes << std::endl << std::endl;
                }
                 */


                // there can not be an overflow ;)
                return 1;

            } else {

                if (this->m_iKmerStarts.getBit(iPos) == 1)
                {
                    this->releaseLock(iThreadID, iPos);
                    continue;
                }

                // the reprobe part must match the number of reprobes back to the previous entry!
                bool bMatchesKey = positionMatchesReprobe(iPos, basekey, iReprobe+iPerformedReprobes);

                if (!bMatchesKey)
                {
                    this->releaseLock(iThreadID, iPos);
                    continue;
                }

                /*
                if (iReprobe >= 1)
                {
                    if ((kmer != NULL) && (false))
                    {
                        std::cout << "Kmer " << kmer->to_string() << std::endl;
                        std::cout << "Accepted Position " << iPos  << " reprobe " << iReprobe << " iperfre " << iPerformedReprobes<< std::endl << std::endl;
                    }
                }
                 */

                CIncrementElement incRet = this->incrementElement(iPos, basekey, iReprobe+iPerformedReprobes, true);

                if (incRet.iOverflow == 0)
                {
                    //no overflow
                    //this->performIncrement(&incRet);
                    //this->releaseLock(iThreadID, iPos);

                    // no overflow occurred, we are all fine :)
                    return 1;
                } else {

                    // get next position lock
                    TSX::tsx_key_t okey = this->getKeyFromKeyVal(incRet.keyval);
                    uint8_t iOverflowHandled = this->handleOverflow(incRet.iOverflow, basekey, iReprobe+iPerformedReprobes, kmer);

                    if (iOverflowHandled == 2)
                    {
                        return 2;
                    }

                    // perform increment
                    //this->performIncrement(&incRet);
                    return 1;

                }

            }

        }

        return 0;
    }


    /**
     *
     * @param kmer which should be added to counter
     * @return
     */
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

            //std::cout << "Empty position " << (int) bEmpty << std::endl;

            if (bEmpty) {


                // we can directly add the element to this position
                TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
                TSX::tsx_key_t updkey = this->makeKey(key, iReprobes);

                // no overflow can happen here ...
                CIncrementElement oINC = this->incrementElement(iPos, key, iReprobes, false);
                CIncrementElement* pINC = &(oINC);
                TSX::tsx_keyval_t savedkey = UBigInt(m_iKeyValBits, true, this->m_pPool);
                TSX::tsx_keyval_t* pSavedKey = &savedkey;

                int iinc=0;

                uint status = _xbegin();
                if(status == _XBEGIN_STARTED) {

                    //this->performIncrement(&incRet);
                    // increment element

                    bool elemsEqual=true;
                    uint64_t iBitsToPos = (pINC->iPosition)*m_iKeyValBits;
                    uint32_t iStartPos = iBitsToPos / (sizeof(FIELDTYPE)*8);
                    uint32_t iStartOffset = iBitsToPos-((sizeof(FIELDTYPE)*8) * iStartPos);


                    FIELDTYPE* pPos = m_pCounterArray + iStartPos;
                    pSavedKey->copy_content_bits(pPos, iStartOffset, m_iKeyValBits);

                    // check elements are equal
                    elemsEqual = pSavedKey->isEqual(pSavedKey, &(pINC->original), m_iKeyValBits);

                    if (!elemsEqual)
                    {
                        _xabort(0xff);
                    }

                    //oStartPos = udiv( (uint32_t) pINC->iPosition*m_iKeyValBits, sizeof(uint8_t) * 8);
                    pINC->keyval.copy_content_to_array(pPos, iStartOffset, m_iKeyValBits);

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
                    } else {
                        //std::cout << "unknown abort empty " << (int) status << std::endl;
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


                    // we can add the element ot this current position, but must handle overflows, etc.
                    // transactions are within incrementElement
                    this->incrementElement(iPos, basekey, iReprobes, false);

                } else {

                    // @TODO why was that there anyhow?
                    //bool bMatchesKey2 = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                    // transaction completed => incorrect position => reprobe
                    ++iReprobes;
                    continue;
                }


            }


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

        return bInserted;

    }

};


#endif //TSXCOUNT_TSXHASHMAPTSXSMALL_H
