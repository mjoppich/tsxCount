//
// Created by mjopp on 11.06.2016.
//

#ifndef TSXCOUNT_TSXHASHMAP_H
#define TSXCOUNT_TSXHASHMAP_H

#include <inttypes.h>
#include <stdint.h>
#include <cmath>
#include <cstdlib>
#include "IBijectiveFunction.h"
#include "TSXTypes.h"
#include "BijectiveKMapping.h"

#include <iostream>
#include <bitset>

using namespace TSX;



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
     * \param iExpLength the actual number of bins will be 2^iExpLength
     * \param iStorage the bits used for each value
     * \param iK the actual K to measure
     *
     */
    TSXHashMap(uint8_t iL, uint32_t iStorageBits, uint16_t iK)
            : m_iL( iL ),
              m_iLength((uint32_t) (std::pow(2, iL))),
              m_iStorageBits( iStorageBits ), //(uint32_t) std::ceil(std::log(iStorageBits))
              m_iK( iK ),
              m_iMapSize( UBigInt::createFromBitShift(iL, iL) ),
              m_iKeyValBits( 2 * iK + iStorageBits ),
              m_iMaxReprobes( (uint32_t) (1 << iL)-1 )
    {

        if (2 * m_iK <= m_iL)
        {
            throw "Invalid lengths for hashmap size and value of k";
        }

        size_t iElements = std::pow(2, m_iL);

        // 2*k = key, m_iStorageBits = value
        size_t iBitsPerElement = 2 * m_iK + m_iStorageBits;
        size_t iBytes = std::ceil( (float)(iElements * iBitsPerElement) / 8.0f);

        std::cerr << "Creating array with " << std::to_string(iBytes) << " bytes for " << std::to_string(iElements) << " places." << std::endl;
        m_pCounterArray = (uint8_t*) calloc( sizeof(uint8_t), iBytes);
        m_iKmerStarts.resize(iElements);


        m_mask_value_key = UBigInt(2*m_iK + m_iStorageBits, true);
        m_mask_value_key.setBit(m_iStorageBits, 1);
        m_mask_value_key = m_mask_value_key-1;

        m_mask_key_value = ~m_mask_value_key;


        m_mask_reprobe_func = UBigInt(2*m_iK, true);
        m_mask_reprobe_func.setBit(m_iL, 1);

        std::cerr << m_mask_reprobe_func.to_string() << std::endl;
        m_mask_reprobe_func = m_mask_reprobe_func-1;
        std::cerr << m_mask_reprobe_func.to_string() << std::endl;


        m_mask_func_reprobe = ~m_mask_reprobe_func;


        std::cerr << "Value Mask (value = 1)" << std::endl;
        std::cerr << m_mask_value_key.to_string() << std::endl;

        std::cerr << "Key Mask (key = 1)" << std::endl;
        std::cerr << m_mask_key_value.to_string() << std::endl;

        std::cerr << "Reprobe Mask (reprobe = 1)" << std::endl;
        std::cerr << m_mask_reprobe_func.to_string() << std::endl;

        std::cerr << "Func Mask (func = 1)" << std::endl;
        std::cerr << m_mask_func_reprobe.to_string() << std::endl;


        m_pHashingFunction = new BijectiveKMapping(m_iK);

    }

    ~TSXHashMap()
    {
        delete m_pHashingFunction;
        free(m_pCounterArray);
    }

    bool addKmer(TSX::tsx_kmer_t& kmer)
    {
        bool bInserted = false;

        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        std::cerr << "new iKmer " << kmer.to_string() << std::endl;
        std::cerr << "new basekey " << basekey.to_string() << std::endl;


        while ( iReprobes < m_iMaxReprobes )
        {

            // get possible position
            uint64_t iPos = this->getPosition( basekey, iReprobes);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                TSX::tsx_key_t key = (basekey & m_mask_func_reprobe);
                std::cerr << "key and inv " << key.to_string() << std::endl;

                this->incrementElement(iPos, key, iReprobes, false);

                // so we can find kmers later without knowing the kmer
                m_iKmerStarts.setBit(iPos, 1);

                bInserted = true;
                break;

            } else {

                bool bIsOtherKmerStart = m_iKmerStarts.getBit(iPos) == 1;
                bool bMatchesKey = positionMatchesKeyAndReprobe(iPos, basekey, iReprobes);

                if ((!bIsOtherKmerStart) && (bMatchesKey))
                {

                    bool bOverflow = this->incrementElement(iPos, basekey, iReprobes, false);

                    if (bOverflow)
                    {
                        handleOverflow(iPos, basekey, iReprobes);
                    }


                    bInserted = true;
                    break;

                } else {

                    ++iReprobes;
                    continue;
                }


            }


        }

        return bInserted;

    }

    UBigInt getKmerCount(TSX::tsx_kmer_t& kmer)
    {
        bool bFound = false;


        uint32_t iReprobes = 1;
        TSX::tsx_key_t basekey = m_pHashingFunction->apply( kmer );

        UBigInt oResult(m_iStorageBits, true);

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

                    uint32_t iUsedOverflowBits = oOverflows.getBitCount();
                    if (iUsedOverflowBits > 0)
                    {
                        oOverflows.resize( iUsedOverflowBits + m_iStorageBits );

                        oOverflows = oOverflows << m_iStorageBits;
                        oOverflows = oOverflows | oResult;

                        oResult = oOverflows;
                    }

                }

            }

            ++iReprobes;


        }

        if (!bFound)
        {
            std::cerr << "Kmer " << kmer.to_string() << " not in hash" << std::endl;
        }

        return oResult;
    }

protected:

    UBigInt findOverflowCounts(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe)
    {
        bool bHandled = false;
        uint32_t iPerformedReprobes = 0;

        UBigInt oReturn = 0;
        uint32_t iAddedBits = 0;

        while (iPerformedReprobes < m_iMaxReprobes)
        {

            iPerformedReprobes += 1;
            iReprobe += 1;

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (!bEmpty)
            {

                // the reprobe part must match the number of reprobes back to the previous entry!
                bool bMatchesKey = positionMatchesReprobe(iPos, basekey, iPerformedReprobes);

                if (!bMatchesKey)
                    continue;

                TSX::tsx_keyval_t elem = this->getElement(iPos);

                UBigInt posValue = this->getFuncValFromKeyVal(elem);

                oReturn = (posValue << iAddedBits) | oReturn;

                iAddedBits += 2*m_iK - m_iL + m_iStorageBits;

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

        TSX::tsx_reprobe_t ret(j);
        ret.resize(m_iL);

        return ret;
    }

    TSX::tsx_key_t makeKey(TSX::tsx_key_t& basekey, uint32_t reprobe)
    {

        TSX::tsx_reprobe_t orepr( reprobe );
        orepr.resize(m_iL);

        return this->makeKey(basekey, orepr);
    }

    TSX::tsx_key_t makeKey(TSX::tsx_key_t& basekey, TSX::tsx_reprobe_t& reprobe)
    {

        TSX::tsx_key_t oRet = (basekey & m_mask_func_reprobe) | reprobe;
        oRet.resize(2*m_iK);

        return oRet;

    }

    inline uint64_t getPosition(TSX::tsx_key_t& hashedKey, uint32_t iReprobes)
    {

        // Uk is still in range [0, 4^k -1]
        hashedKey = hashedKey + this->reprobe(iReprobes);

        std::cerr << "kmer key: " << hashedKey.to_string() << std::endl;

        UBigInt mod2 = hashedKey.mod2( m_iL );

        std::cerr << "kmer key % 2: " << mod2.to_string() << std::endl;

        // Uk.mod2 now is in range [0, M-1] = [0, 2^l -1]
        uint64_t iPos = mod2.toUInt();

        return iPos;
    }

    bool positionMatchesKeyAndReprobe(uint64_t pos, TSX::tsx_key_t& key, uint32_t iReprobe)
    {

        TSX::tsx_keyval_t oKeyVal = getElement(pos);

        // the key matches
        bool keyPartMatch = (oKeyVal >> (m_iL + m_iStorageBits)) == (key >> (m_iL));

        // and the reprobe matches
        bool reprobePartMatch = this->positionMatchesReprobe(pos, key, iReprobe);

        if (keyPartMatch && reprobePartMatch)
        {
            return true;
        } else {

            std::cerr << key.to_string() << std::endl;
            std::cerr << oKeyVal.to_string() << std::endl;

            return false;
        }
    }

    bool positionMatchesReprobe(uint64_t pos, TSX::tsx_key_t& key, uint32_t iReprobe)
    {

        // get key from position
        TSX::tsx_keyval_t oKeyVal = getElement(pos);
        TSX::tsx_key_t oKey = this->getKeyFromKeyVal(oKeyVal);

        TSX::tsx_key_t oReprobe = UBigInt::fromUint32(iReprobe);
        oReprobe.resize(2*m_iK);

        // extract position reprobe and compare with reprobe value
        bool reprobeMatch = (oKey & m_mask_reprobe_func) == oReprobe;

        if (reprobeMatch)
            return true;

        return false;

    }

    /**
     * checks if position has been yet used or is empty
     * @param iPos
     * @return true if position was never used
     */
    bool positionEmpty(uint64_t iPos)
    {

        tsx_keyval_t element = this->getElement(iPos);

        return element == 0;
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
        std::div_t oStartPos = div( pos*m_iKeyValBits, iBitsPerField);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt oReturn = UBigInt::createFromMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, m_iKeyValBits);

        std::cerr << "keyval fetched: " << oReturn.to_string() << " from array " << std::to_string(iArray) << " in pos " << std::to_string(pos) << std::endl;

        return oReturn;

    }

    void storeElement(uint64_t pos, TSX::tsx_key_t& key, TSX::tsx_val_t& val)
    {

        TSX::tsx_keyval_t oKeyVal(key, 2*m_iK + m_iStorageBits);
        oKeyVal = (oKeyVal << m_iStorageBits) | val;



        uint8_t iArray = 0;//pos % 4;
        //pos = pos >> 4;


        // we want to get the iPosition-th entry
        uint32_t iBitsPerField = sizeof(uint8_t) * 8;
        std::div_t oStartPos = div( pos*m_iKeyValBits, iBitsPerField);
        oKeyVal.resize(m_iKeyValBits);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt::storeIntoMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, oKeyVal, m_iKeyValBits);

        TSX::tsx_keyval_t oStoredKeyVal = this->getElement(pos);

        bool bIsStoredCorrectly = (oStoredKeyVal == oKeyVal);

        if (!bIsStoredCorrectly)
        {
            std::cerr << "Should be: " << oKeyVal.to_string() << std::endl;
            std::cerr << "but is   : " << oStoredKeyVal.to_string() << std::endl;
            UBigInt::storeIntoMemory(m_pCounterArray, oStartPos.quot, oStartPos.rem, oKeyVal, m_iKeyValBits);
            oStoredKeyVal = this->getElement(pos);
        }

        std::cerr << "keyval stored: " << oKeyVal.to_string() << " into array " << std::to_string(iArray) << " in pos " << std::to_string(pos);
        std::cerr << std::endl;

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
        ret = ret << (m_iL + m_iStorageBits);
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

    /**
     *
     * @param iPosition position in array
     * @param key key to look for
     * @param reprobes amount reprobes used
     * @param bKeyIsValue if true, the key part belongs to value
     * @return
     */
    bool incrementElement(uint64_t iPosition, TSX::tsx_key_t& key, uint32_t reprobes, bool bKeyIsValue)
    {

        TSX::tsx_keyval_t keyval = getElement(iPosition);
        TSX::tsx_val_t value = this->getValFromKeyVal(keyval);


        if (~value == 0)
        {
            // an overflow is going to happen!
            if (~value == 0)
            {
                std::cerr << "is zero" << std::endl;
            }

            value = 1;
            value.resize(m_iStorageBits);

            // maybe we can prevent the overflow?
            if (bKeyIsValue)
            {

                TSX::tsx_func_t funcpart = this->getFuncFromKeyVal(keyval);

                if (~funcpart == 0)
                {
                    // overflow in func key part => set funcpart = 0 and propagate overflow further
                    funcpart = 0;
                    TSX::tsx_key_t updkey = (UBigInt(funcpart, 2*m_iK) << m_iL) | reprobes;

                    storeElement(iPosition, updkey, value);

                    return true;

                } else {

                    funcpart += 1;

                    TSX::tsx_key_t updkey = (UBigInt(funcpart, 2*m_iK) << m_iL) | reprobes;
                    value = 0;
                    value.resize(m_iStorageBits);

                    storeElement(iPosition, updkey, value);
                    return false;
                }

            }

            // also set value == 0
            TSX::tsx_key_t updkey = (key & m_mask_func_reprobe) | reprobes;
            value = 0;
            value.resize(m_iStorageBits);

            storeElement(iPosition, updkey, value);

            // overflow occurred!
            return true;
        } else {

            value += 1;

            std::cerr << "Storing for key: " << key.to_string() << " value " << value.to_string() << std::endl;

            if (bKeyIsValue)
            {
                TSX::tsx_key_t updkey(2*m_iK, true);
                updkey = updkey | reprobes;
                storeElement(iPosition, updkey, value);
            } else {
                TSX::tsx_key_t updkey = this->makeKey(key, reprobes);
                storeElement(iPosition, updkey, value);
            }


            return false;

        }
    }

    bool handleOverflow(uint64_t iPos, TSX::tsx_key_t& basekey, uint32_t iReprobe)
    {
        bool bHandled = false;
        uint32_t iPerformedReprobes = 0;

        while (iPerformedReprobes < m_iMaxReprobes)
        {

            iPerformedReprobes += 1;
            iReprobe += 1;

            // this fetches the element using the global reprobe!
            uint64_t iPos = this->getPosition(basekey, iReprobe);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                // the reprobe part must match the number of reprobes back to the previous entry!
                this->incrementElement(iPos, basekey, iPerformedReprobes, true);

                // there can not be an overflow ;)
                return true;

            } else {

                // the reprobe part must match the number of reprobes back to the previous entry!
                bool bMatchesKey = positionMatchesReprobe(iPos, basekey, iPerformedReprobes);

                if (!bMatchesKey)
                    continue;

                bool bOverflow = this->incrementElement(iPos, basekey, iPerformedReprobes, true);

                if (!bOverflow)
                {
                    // no overflow occurred, we are all fine :)
                    return true;
                } else {

                    // reset performed reprobes as this should indicate number of reprobes needed!
                    iPerformedReprobes = 0;
                }

            }

        }

        return false;

    }


    /*
     *
     * Variables
     *
     */


    uint8_t* m_pCounterArray;

    TSX::tsx_keyval_t m_mask_key_value = 0;
    TSX::tsx_keyval_t m_mask_value_key = 0;

    TSX::tsx_key_t m_mask_reprobe_func = 0;
    TSX::tsx_key_t m_mask_func_reprobe = 0;


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


};


#endif //TSXCOUNT_TSXHASHMAP_H
