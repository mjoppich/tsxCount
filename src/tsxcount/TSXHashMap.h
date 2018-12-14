//
// Created by mjopp on 11.06.2016.
//

#ifndef TSXCOUNT_TSXHASHMAP_H
#define TSXCOUNT_TSXHASHMAP_H

#include <cinttypes>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include "IBijectiveFunction.h"
#include "TSXTypes.h"
#include "BijectiveKMapping.h"

#include <iostream>
#include <bitset>
#include <cstdlib>
#include <vector>
#include <set>
#include <omp.h>

using namespace TSX;



class TSXException: public std::exception
{
public:
    TSXException(std::string sText)
            : std::exception(), m_sText(sText)
    {

    }

    virtual const char* what() const throw()
    {
        return m_sText.c_str();
    }


protected:

    const std::string m_sText;

};

/**
 * \class TSXHashMap
 *
 *
 * \brief implements the main hashmap for tsxcount:
 *
 * This class implements a HashMap for storing key/value pairs.
 *
 * \note Performance check this
 *
 * \author Markus Joppich
 *
 */
class TSXHashMap {
public:

    /**
     *
     * \param iL the actual number of bins will be 2^iExpLength
     * \param iStorageBits the bits used for each value
     * \param iK the actual K to measure
     * \note 2*iK > m_iL
     *
     */
    TSXHashMap(uint8_t iL, uint32_t iStorageBits, uint16_t iK)
            : m_pPool( new MemoryPool<FIELDTYPE>(1) ),
              m_iL( iL ),
              m_iLength((uint32_t) (std::pow(2, iL))),
              m_iStorageBits( iStorageBits ), //(uint32_t) std::ceil(std::log(iStorageBits))
              m_iKeyValBits( 2 * iK + iStorageBits ),
              m_iK( iK ),
              m_iMaxReprobes( (uint32_t) (1 << iL)-1 ),
              m_iMapSize( UBigInt::createFromBitShift(iL, iL, this->m_pPool))
    {

        if (2 * m_iK <= m_iL)
        {
            throw TSXException("Invalid lengths for hashmap size and value of k");
        }

        size_t iElements = std::pow(2, m_iL);

        // 2*k = key, m_iStorageBits = value
        size_t iBitsPerElement = 2 * m_iK + m_iStorageBits;
        size_t iBytes = std::ceil( (float)(iElements * iBitsPerElement) / 8.0f);

        std::cerr << "Creating array with " << std::to_string(iBytes) << " bytes for " << std::to_string(iElements) << " places." << std::endl;
        m_pCounterArray = (FIELDTYPE*) calloc( sizeof(FIELDTYPE), iBytes);

        m_iKmerStarts = UBigInt(0, m_pPool);
        m_iKmerStarts.resize(iElements);


        m_mask_value_key = UBigInt(2*m_iK + m_iStorageBits, true, this->m_pPool);
        m_mask_value_key.setBit(m_iStorageBits, 1);
        m_mask_value_key = m_mask_value_key-1;

        m_mask_key_value = ~m_mask_value_key;


        m_mask_reprobe_func = UBigInt(2*m_iK, true, this->m_pPool);
        m_mask_reprobe_func.setBit(m_iL, 1);

        m_mask_reprobe_func = m_mask_reprobe_func-1;
        m_mask_func_reprobe = ~m_mask_reprobe_func;


        std::cerr << "Value Mask (value = 1) m_mask_value_key" << std::endl;
        std::cerr << m_mask_value_key.to_string() << std::endl;

        std::cerr << "Key Mask (key = 1) m_mask_key_value" << std::endl;
        std::cerr << m_mask_key_value.to_string() << std::endl;

        std::cerr << "Reprobe Mask (reprobe = 1) m_mask_reprobe_func" << std::endl;
        std::cerr << m_mask_reprobe_func.to_string() << std::endl;

        std::cerr << "Func Mask (func = 1) m_mask_func_reprobe" << std::endl;
        std::cerr << m_mask_func_reprobe.to_string() << std::endl;


        m_pHashingFunction = new BijectiveKMapping(m_iK, this->m_pPool);


        this->setThreads(2);

    }

    virtual ~TSXHashMap()
    {
        delete m_pHashingFunction;
        free(m_pCounterArray);
    }


    MemoryPool<FIELDTYPE>* getMemoryPool()
    {
        return m_pPool;
    }

    uint32_t getK()
    {
        return m_iK;
    }

    size_t getUsedPositions()
    {
        return m_setUsedPositions.size();
    }

    virtual bool addKmer(TSX::tsx_kmer_t& kmer)
    {
        bool bInserted = false;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        uint8_t iThreadID = omp_get_thread_num();

        while ( iReprobes < m_iMaxReprobes )
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            // does this position match to key?

            bool lockAcquired = this->acquireLock(iThreadID, iPos);

            if (!lockAcquired)
            {
                continue;
            }

            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
                TSX::tsx_key_t updkey = this->makeKey(key, iReprobes);

                // no overflow can happen here ...
                CIncrementElement incRet = this->incrementElement(iPos, key, iReprobes, false);
                this->performIncrement(&incRet);

                // so we can find kmer start positions later without knowing the kmer
                m_iKmerStarts.setBit(iPos, 1);
                // TODO this slows down inserting, but is a nice measure ifdef out verbose?
                m_setUsedPositions.insert(iPos);

                this->releaseLock(iThreadID, iPos);

                bInserted = true;
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

                    CIncrementElement incRet = this->incrementElement(iPos, basekey, iReprobes, false);
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

                    /*
                    if (vOPS.size() > 2)
                    {
                        std::cout << "kmer " <<  kmer.to_string() << std::endl;
                    }
                    */

                    std::vector<size_t> usedPos;
                    for (auto rit = vOPS.rbegin(); rit != vOPS.rend(); ++rit)
                    {
                        CIncrementElement inc = (*rit);

                        /*
                        if (vOPS.size() > 2)
                        {
                            std::cout << "B" <<  inc.iPosition << " " << this->getElement(inc.iPosition).to_string() << " " << inc.iOverflow << std::endl;
                        }
                        */

                        //std::cout << inc.iPosition;
                        //inc.keyval.print_string();
                        this->performIncrement(&inc);
                        this->releaseLock(iThreadID, inc.iPosition);

                        /*
                        if (vOPS.size() > 2)
                        {
                            std::cout << "A" << inc.iPosition << " " << inc.keyval.to_string() << std::endl;
                        }
                         */

                    }

                    /*
                    if (vOPS.size() > 2)
                    {

                        std::cout << "BaseKey" << basekey.to_string() << std::endl;
                        for (int ii=1; ii < 6; ++ii)
                        {
                            std::cout << ii << " " << this->getPosition(basekey, ii) << std::endl;
                        }


                        std::cout << std::endl;
                    }
                     */

                    //this->releaseLock(iThreadID, iPos);

                    bInserted = true;
                    break;

                } else {

                    bool bMatchesKey2 = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);
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

    UBigInt getKmerCount(TSX::tsx_kmer_t& kmer)
    {
        bool bFound = false;


        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        UBigInt oResult(64, true, this->m_pPool);

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
                    UBigInt oOverflows = findOverflowCounts(iPos, basekey, iReprobes);

                    //std::cout << "Looking at pos " << iPos << std::endl;
                    //std::cout << oResult.to_string() << std::endl;
                    //std::cout << oOverflows.to_string() << std::endl;

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
            std::cerr << "Kmer " << kmer.to_string() << " not in hash" << std::endl;
        }

        return oResult;
    }

    uint64_t getKmerCount()
    {
        return this->m_iKmerStarts.sumBits();
    }

    std::vector<TSX::tsx_kmer_t> getAllKmers()
    {

        std::vector<TSX::tsx_kmer_t> oKmers;

        for (uint64_t i = 0; i < m_iKmerStarts.getBitCount(); ++i)
        {

            uint8_t iKmerStarts = m_iKmerStarts.getBit(i);

            bool isStart = m_setUsedPositions.find(i) != m_setUsedPositions.end();

            if (!iKmerStarts)
                continue;

            //std::cerr << "Position: " << i << std::endl;

            TSX::tsx_keyval_t oKeyVal = this->getElement(i);
            //std::cerr << "KeyVal: " << oKeyVal.to_string() << std::endl;

            TSX::tsx_key_t oKey = this->getKeyFromKeyVal(oKeyVal);
            //std::cerr << "Key: " << oKey.to_string() << std::endl;

            TSX::tsx_reprobe_t oReprobe = this->getReprobeFromKey(oKey);
            uint64_t iReprobeLength = oReprobe.getBitCount();

            uint64_t iReprobe = oReprobe.toUInt();
            oReprobe = this->reprobe(iReprobe);


            UBigInt oMissingPart = UBigInt(i, this->m_pPool);
            oMissingPart.resize(iReprobeLength);

            //std::cerr << "Position: " << oMissingPart.to_string() << std::endl;
            //std::cerr << "Reprobe: " << oReprobe.to_string() << std::endl;

            oMissingPart = oMissingPart-oReprobe;

            //std::cerr << "pos-reprobe: " << oMissingPart.to_string() << std::endl;

            oMissingPart = oMissingPart.mod2(m_iL);
            oMissingPart.resize(m_iL);

            //std::cerr << "Missing: " << oMissingPart.to_string() << " " << oMissingPart.getBitCount() << std::endl;
            //std::cerr << "Key: " << oKey.to_string() << " " << oKey.getBitCount() << std::endl;


            UBigInt oHKmer = oKey & ~m_mask_reprobe_func;
            oHKmer = oHKmer | oMissingPart;

            //std::cerr << "HKmer: " << oHKmer.to_string() << " " << oHKmer.getBitCount() << std::endl;

            TSX::tsx_kmer_t oKmer = this->m_pHashingFunction->inv_apply(oHKmer);

            //std::cerr << "Kmer: " << oKmer.to_string() << std::endl;

            oKmers.push_back(oKmer);

        }

        return oKmers;

    }

    void testHashFunction() {

        UBigInt oTest = UBigInt::fromString("11001001000000001011", this->m_pPool);

        TSX::tsx_key_t oKey = m_pHashingFunction->apply(oTest);
        UBigInt oInvKey = m_pHashingFunction->inv_apply(oKey);

        std::cout << oTest.to_string() << std::endl;
        std::cout << oInvKey.to_string() << std::endl;

        std::cout << std::endl;
    }

    int getThreads() {
        return m_iThreads;
    }

    void setThreads(uint8_t iThreads)
    {
        m_iThreads = iThreads;
        this->m_pPool->setThreads(iThreads);
        this->initialiseLocks();
    }

protected:



    virtual void initialiseLocks()
    {
    }


    /**
     *
     * @param iArrayPos checks whether iArrayPos is locked
     * @return threadID of thread who locks iArrayPos or -1
     */
    virtual uint8_t position_locked(uint64_t iArrayPos)
    {
        return omp_get_thread_num();
    }

    virtual bool canAcquireLock(uint8_t iThreadID, uint64_t iArrayPos)
    {
        return true;
    }

    /**
     *
     * @param iThreadID thread id for which the lock is acquired
     * @param iArrayPos array position for which the lock is acquired
     * @return True if lock successfully acquired, False otherwise
     */
    virtual bool acquireLock(uint8_t iThreadID, uint64_t iArrayPos)
    {
        return true;
    }


    /**
     *
     * unlocks all locked array positions for given thread
     *
     * @param iThreadID
     */
    virtual void unlock_thread(uint8_t iThreadID)
    {
    }

    virtual bool releaseLock(uint8_t iThreadID, uint64_t iPos)
    {
        return true;
    }

    UBigInt findOverflowCounts(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe)
    {
        bool bHandled = false;
        uint32_t iPerformedReprobes = 0;

        UBigInt oReturn(0, m_pPool);
        uint32_t iRequiredBits = 0;

        while (iPerformedReprobes < 10)
        {

            iPerformedReprobes += 1;
            iReprobe += 1;

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe);

            //std::cout << "OVFL Looking at pos " << iPos << " " << iReprobe << std::endl;

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (!bEmpty)
            {

                if (m_iKmerStarts.getBit(iPos) == 1)
                    continue;
                //std::cout << "OVFL B Match Key " << iPos << " " << iReprobe << std::endl;

                // the reprobe part must match the number of reprobes back to the previous entry!
                bool bMatchesKey = positionMatchesReprobe(iPos, basekey, iReprobe);

                if (!bMatchesKey)
                    continue;

                //std::cout << "OVFL A Match Key " << iPos << " " << iReprobe << std::endl;
                TSX::tsx_keyval_t elem = this->getElement(iPos);

                UBigInt posValue = this->getFuncValFromKeyVal(elem);

                //std::cout << "PV " << posValue.to_string() << std::endl;

                uint32_t iOldRequired = iRequiredBits;
                iRequiredBits += 2*m_iK - m_iL + m_iStorageBits;
                posValue.resize(iRequiredBits);

                oReturn = (posValue << ( iOldRequired )) | oReturn;


                // reset performed reprobes as this should indicate number of reprobes needed!
                // now we can try to find further matching positions :)
                iPerformedReprobes = 0;
            }

        }

        // no more matching position as max number of reprobes reached without finding a matching one
        return oReturn;

    }

    void printBinary(int x)
    {
        std::cout << std::bitset<32>(x) << std::endl;
    }

    TSX::tsx_reprobe_t reprobe(uint32_t i)
    {
        uint32_t j = i * (i+1) / 2;

        TSX::tsx_reprobe_t ret(j, this->m_pPool);
        ret.resize(m_iL);

        return ret;
    }

    TSX::tsx_key_t makeKey(TSX::tsx_key_t& basekey, uint32_t reprobe)
    {

        TSX::tsx_reprobe_t orepr( reprobe, this->m_pPool );
        orepr.resize(m_iL);

        return this->makeKey(basekey, orepr);
    }

    TSX::tsx_key_t makeKey(TSX::tsx_key_t& basekey, TSX::tsx_reprobe_t& reprobe)
    {

        TSX::tsx_key_t oRet = (basekey & m_mask_func_reprobe) | reprobe;
        oRet.resize(2*m_iK);
        return oRet;

    }

    inline uint64_t getPosition(TSX::tsx_key_t& basekey, uint32_t iReprobes)
    {


        UBigInt hashedKey(basekey);
        
        // Uk is still in range [0, 4^k -1]
        hashedKey = hashedKey + this->reprobe(iReprobes);

        //std::cerr << "kmer key: " << hashedKey.to_string() << std::endl;

        UBigInt mod2 = hashedKey.mod2( m_iL );

        //std::cerr << "kmer key % 2: " << mod2.to_string() << std::endl;

        // Uk.mod2 now is in range [0, M-1] = [0, 2^l -1]
        uint64_t iPos = mod2.toUInt();

        return iPos;
    }

    bool positionMatchesKeyAndReprobe(uint64_t pos, TSX::tsx_key_t& key, uint32_t iReprobe)
    {

        TSX::tsx_keyval_t oKeyVal = getElement(pos);

        TSX::tsx_func_t oKeyValFunc = (oKeyVal >> (m_iL + m_iStorageBits));
        TSX::tsx_func_t oKeyFunc = (key >> (m_iL));

        // the key matches
        bool keyPartMatch = ( oKeyValFunc == oKeyFunc );

        // and the reprobe matches
        bool reprobePartMatch = this->positionMatchesReprobe(pos, key, iReprobe);

        if (keyPartMatch && reprobePartMatch)
        {
            return true;
        } else {

            //std::cerr << key.to_string() << std::endl;
            //std::cerr << oKeyVal.to_string() << std::endl;

            return false;
        }
    }

    bool positionMatchesReprobe(uint64_t pos, TSX::tsx_key_t& key, uint32_t iReprobe)
    {

        // get key from position
        TSX::tsx_keyval_t oKeyVal = getElement(pos);
        TSX::tsx_key_t oKey = this->getKeyFromKeyVal(oKeyVal);

        TSX::tsx_key_t oReprobe(iReprobe, this->m_pPool);
        oReprobe.resize(2*m_iK);

        // extract position reprobe and compare with reprobe value
        bool reprobeMatch = (oKey & m_mask_reprobe_func) == oReprobe;

        return reprobeMatch;

    }

    /**
     * checks if position has been yet used or is empty
     * @param iPos
     * @return true if position was never used
     */
    bool positionEmpty(uint64_t iPos)
    {
        tsx_keyval_t element = this->getElement(iPos);

        return element.isZero();
    }

    inline std::div_t udiv(uint32_t a, uint32_t b)
    {

        std::div_t ret;

        ret.quot = a / b;
        ret.rem = a % b;

        return ret;
    }


    /**
     *
     * @param pos
     * @return keyval of element at position pos in array
     */
    inline TSX::tsx_keyval_t getElement(uint64_t pos)
    {

        uint8_t iArray = 0;//pos % 4;
        //pos = pos >> 4;

        // we want to get the iPosition-th entry
        uint32_t iBitsPerField = sizeof(uint8_t) * 8;
        uint32_t iDivPos = pos * m_iKeyValBits;

        std::div_t oStartPos = udiv( iDivPos, iBitsPerField);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt oReturn = UBigInt::createFromMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, m_iKeyValBits, this->m_pPool);

        //std::cerr << "keyval fetched: " << oReturn.to_string() << " from array " << std::to_string(iArray) << " in pos " << std::to_string(pos) << std::endl;

        return oReturn;

    }

    inline void getElement(uint64_t pos, TSX::tsx_keyval_t* pReturn)
    {
        // we want to get the iPosition-th entry
        uint32_t iBitsPerField = sizeof(FIELDTYPE) * 8;
        uint64_t iDivPos = pos * m_iKeyValBits;

        std::div_t oStartPos = udiv( iDivPos, iBitsPerField);

        FIELDTYPE* pPos = m_pCounterArray + oStartPos.quot;
        pReturn->copy_content_bits(pPos, (uint32_t) oStartPos.rem, m_iKeyValBits);

    }

    inline void getElementTest(uint64_t pos, TSX::tsx_keyval_t* pReturn)
    {
        // we want to get the iPosition-th entry
        uint32_t iBitsPerField = sizeof(FIELDTYPE) * 8;
        uint64_t iDivPos = pos * m_iKeyValBits;

        std::div_t oStartPos = udiv( iDivPos, iBitsPerField);

        FIELDTYPE* pPos = m_pCounterArray + oStartPos.quot;
        pReturn->copy_content_bits(pPos, (uint32_t) oStartPos.rem, m_iKeyValBits);

    }

    inline
    void storeKeyValElement(uint64_t pos, TSX::tsx_keyval_t& oKeyVal)
    {
        // we want to get the iPosition-th entry
        std::div_t oStartPos = udiv( (uint32_t) pos*m_iKeyValBits, sizeof(uint8_t) * 8);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt::storeIntoMemory(m_pCounterArray, (uint32_t) oStartPos.quot, (uint32_t) oStartPos.rem, oKeyVal, m_iKeyValBits);
    }


    void storeElement(uint64_t pos, TSX::tsx_key_t& key, TSX::tsx_val_t& val)
    {

        TSX::tsx_keyval_t oKeyVal(key, 2*m_iK + m_iStorageBits);
        oKeyVal = (oKeyVal << m_iStorageBits) | val;

        uint8_t iArray = 0;//pos % 4;
        //pos = pos >> 4;

        // we want to get the iPosition-th entry
        uint32_t iBitsPerField = sizeof(uint8_t) * 8;
        std::div_t oStartPos = udiv( pos*m_iKeyValBits, iBitsPerField);
        oKeyVal.resize(m_iKeyValBits);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt::storeIntoMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, oKeyVal, m_iKeyValBits);


        /*
        TSX::tsx_keyval_t oStoredKeyVal = this->getElement(pos);
        bool bIsStoredCorrectly = (oStoredKeyVal == oKeyVal);

        if (!bIsStoredCorrectly)
        {
            exit(111);

            std::cerr << "Should be: " << oKeyVal.to_string() << std::endl;
            std::cerr << "but is   : " << oStoredKeyVal.to_string() << std::endl;
            UBigInt::storeIntoMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, oKeyVal, m_iKeyValBits);
            oStoredKeyVal = this->getElement(pos);

        }

        */

        //std::cerr << "keyval stored: " << oKeyVal.to_string() << " into array " << std::to_string(iArray) << " in pos " << std::to_string(pos) << std::endl;

    }

    TSX::tsx_val_t getValFromKeyVal(TSX::tsx_keyval_t& keyval)
    {

        TSX::tsx_val_t ret = keyval &  m_mask_value_key ;
        ret.resize(m_iStorageBits);

        return ret;

    }

    TSX::tsx_key_t getKeyFromKeyVal(TSX::tsx_keyval_t& keyval)
    {
        TSX::tsx_key_t ret = keyval >> m_iStorageBits;
        ret.resize(2*m_iK);

        return ret;
    }

    TSX::tsx_func_t getFuncFromKey(TSX::tsx_key_t& key)
    {

        TSX::tsx_func_t ret = key & m_mask_func_reprobe;
        ret.resize(2*m_iK - m_iL);

        return ret;

    }

    TSX::tsx_reprobe_t getReprobeFromKey(TSX::tsx_key_t& key)
    {

        TSX::tsx_reprobe_t ret = key & m_mask_reprobe_func;
        ret.resize(m_iL);

        return ret;

    }

    TSX::tsx_func_t getFuncValFromKeyVal(TSX::tsx_keyval_t& keyval)
    {

        TSX::tsx_func_t ret = (keyval >> (m_iL + m_iStorageBits));
        ret = ret << ( m_iStorageBits );
        ret = ret | this->getValFromKeyVal(keyval);

        ret.resize(2*m_iK-m_iL + m_iStorageBits);

        return ret;
    }

    TSX::tsx_func_t getFuncFromKeyVal(TSX::tsx_keyval_t& keyval)
    {

        TSX::tsx_func_t ret = (keyval >> (m_iL + m_iStorageBits));
        ret.resize(2*m_iK-m_iL);

        return ret;
    }


    struct CIncrementElement
    {
        uint8_t iOverflow;
        uint64_t iPosition;

        //TSX::tsx_key_t okey;
        //TSX::tsx_val_t oval;

        TSX::tsx_keyval_t keyval;
        TSX::tsx_keyval_t original;
    };

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

                    return oRet;
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

                    return oRet;
                    //return 0;
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

            return oRet;
            //return 1;
        } else {

            bool bEmpty = this->positionEmpty( iPosition );

            ++value;

            //std::cerr << "Storing for key: " << key.to_string() << " value " << value.to_string() << std::endl;

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

            return oRet;
            //return 0;

        }
    }


    void performIncrement(CIncrementElement* pElement)
    {
        //storeElement(pElement->iPosition, pElement->okey, pElement->oval);
        this->storeKeyValElement(pElement->iPosition, pElement->keyval);
    }



    /**
     *
     * @param iPos
     * @param basekey
     * @param iReprobe
     * @return 1 if succeeded, 0 if not succeeded, 2 if blocked
     */
    virtual uint8_t handleOverflow(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe, std::vector<CIncrementElement>* pOPS, UBigInt* kmer=NULL)
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
                CIncrementElement incRet = this->incrementElement(iPos, basekey, iReprobe+iPerformedReprobes, true);

                // no overflow can happen!
                //this->performIncrement(&incRet);
                pOPS->insert(pOPS->end(), incRet);

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
                pOPS->insert(pOPS->end(), incRet);

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
                    uint8_t iOverflowHandled = this->handleOverflow(incRet.iOverflow, basekey, iReprobe+iPerformedReprobes, pOPS, kmer);

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


    /*
     *
     * Variables
     *
     */


    FIELDTYPE* m_pCounterArray;
    MemoryPool<FIELDTYPE>* m_pPool;

    /*
     *
     * Variables needed for the reading/writing of key,value pairs
     *
     */
    const uint8_t m_iL;
    const uint32_t m_iLength;
    const uint32_t m_iStorageBits;
    const uint32_t m_iK;
    const uint32_t m_iKeyValBits;
    const uint32_t m_iMaxReprobes;

    const UBigInt m_iMapSize;

    UBigInt m_iKmerStarts;

    IBijectiveFunction* m_pHashingFunction;

    std::set<uint64_t> m_setUsedPositions;


    uint8_t m_iThreads;

    std::vector<uint64_t>* m_pLocked;
    pthread_mutex_t m_oLockMutex;
    omp_lock_t m_oOMPLock;


    TSX::tsx_keyval_t m_mask_key_value = UBigInt(0, m_pPool);
    TSX::tsx_keyval_t m_mask_value_key = UBigInt(0, m_pPool);

    TSX::tsx_key_t m_mask_reprobe_func = UBigInt(0, m_pPool);
    TSX::tsx_key_t m_mask_func_reprobe = UBigInt(0, m_pPool);

};


#endif //TSXCOUNT_TSXHASHMAP_H
