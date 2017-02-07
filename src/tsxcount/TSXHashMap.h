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

    std::string toBinStr(uint64_t val)
    {
        // mask has only the leftmost bit set.
        uint64_t iShift = std::numeric_limits<uint64_t>::digits-1;
        uint64_t mask = ((uint64_t)1) << (iShift) ;

        // skip leading bits that are not set.
        /*
        while ( 0 == (val & mask) && mask != 0 )
            mask >>= 1 ; // shift all bits to the right by 1
            */

        std::string binStr ;
        binStr.reserve(std::numeric_limits<uint64_t>::digits+1) ;

        do
        {
            // add a '1' or '0' depending the current bit.
            uint8_t iVal = static_cast<char>(val & mask);
            binStr += iVal + '0' ;

        } while ( (mask>>=1) != 0 ) ; // next bit, when mask is 0 we've processed all bits

        return binStr ;
    }

    /**
     *
     * \param iExpLength the actual number of bins will be 2^iExpLength
     * \param iStorage the bits used for each value
     * \param iK the actual K to measure
     *
     */
    TSXHashMap(uint8_t iExpLength, uint32_t iStorage, uint16_t iK)
            : m_iLengthExp( iExpLength ),
              m_iLength((uint32_t) (std::pow(2, iExpLength))),
              m_iStorageBits( (uint32_t) std::ceil(std::log(iStorage)) ),
              m_iK( iK ),
              m_iMapSize( UBigInt::createFromBitShift(iExpLength, iExpLength) ),
              m_iKeyValBits( 2 * iK - iExpLength + iStorage ),
              m_iReprobeBits( iExpLength ), m_iMaxReprobes( (uint32_t) (1 << iExpLength)-1 )
    {

        if (2 * m_iK <= iExpLength)
        {
            throw "Invalid lengths for hashmap size and value of k";
        }

        // this should always be a power of 2 !
        m_iCounterArrays = 4;
        m_pCounterArray= (uint8_t**) calloc( sizeof(uint8_t*), m_iCounterArrays );

        for (int i = 0; i < m_iCounterArrays; ++i)
        {

            size_t iElements = std::pow(2, m_iLengthExp-2);

            // 2*k = key, m_iStorageBits = value
            size_t iBitsPerElement = 2 * m_iK + m_iStorageBits;
            size_t iBytes = std::ceil( (float)(iElements * iBitsPerElement) / 8.0f);

            std::cerr << "Creating array " << std::to_string(i) << " with " << std::to_string(iBytes) << " bytes for " << std::to_string(iElements) << " places." << std::endl;

            m_pCounterArray[i] = (uint8_t*) calloc( sizeof(uint8_t), iBytes);

        }



        uint64_t iBytesUsed = sizeof(tsx_key_t);

        for (size_t i = 0; i < iBytesUsed*8; ++i)
        {

            if (i < (2*m_iK - m_iLengthExp))
            {
                m_iKeyMask = m_iKeyMask | 1;
            } else {
                m_iKeyMask = m_iKeyMask;
            }

            m_iKeyMask = m_iKeyMask << 1;


        }


        m_iReprobeMask = UBigInt(2*m_iK, true);
        m_iKeyMask = UBigInt(2*m_iK, true);
        m_iKeyMask = (~m_iKeyMask >> m_iReprobeBits) << m_iReprobeBits;
        m_iReprobeMask = ~m_iKeyMask;

        std::cerr << "Key Mask" << std::endl;
        std::cerr << m_iKeyMask.to_string() << std::endl;

        std::cerr << "Reprobe Mask" << std::endl;
        std::cerr << m_iReprobeMask.to_string() << std::endl;

        m_pHashingFunction = new BijectiveKMapping(m_iK);

    }

    ~TSXHashMap()
    {
        delete m_pHashingFunction;
        free(m_pCounterArray);
    }

    bool addKmer(TSX::tsx_kmer_t& kmer)
    {
        return this->insert_kmer(kmer);
    }

protected:

    void printBinary(int x)
    {
        std::cout << std::bitset<32>(x) << std::endl;
    }

    tsx_reprobe_t reprobe(uint32_t i)
    {
        return i*(i+1)/2;
    }

    inline uint64_t getPosition(TSX::tsx_kmer_t& iKmer, TSX::tsx_reprobe_t iReprobe)
    {

        // Uk is still in range [0, 4^k -1]
        TSX::tsx_kmer_t Uk = m_pHashingFunction->apply(iKmer) + iReprobe;

        std::cerr << "kmer key: " << Uk.to_string() << std::endl;

        UBigInt mod2 = Uk.mod2( m_iLengthExp );

        std::cerr << "kmer key % 2: " << mod2.to_string() << std::endl;

        // Uk.mod2 now is in range [0, M-1] = [0, 2^l -1]
        uint64_t iPos = mod2.toUInt();

        return iPos;
    }

    TSX::tsx_key_t extractKeyNoReprobe(TSX::tsx_keyval_t& keyval)
    {
        return keyval >> (m_iStorageBits + m_iReprobeBits);
    }

    TSX::tsx_key_t extractKey(TSX::tsx_keyval_t& keyval)
    {
        return keyval >> m_iStorageBits;
    }


    bool positionMatchesKey(uint64_t pos, TSX::tsx_key_t& key, TSX::tsx_reprobe_t iReprobe)
    {

        TSX::tsx_keyval_t oKeyVal = getElement(pos);

        // the key matches
        bool keyPartMatch = (oKeyVal >> (m_iReprobeBits + m_iStorageBits)) == (key >> (m_iReprobeBits + m_iStorageBits));
        // and the reprobe matches
        bool reprobePartMatch = this->positionMatchesReprobe(pos, key, iReprobe);

        return keyPartMatch && reprobePartMatch;
    }

    bool positionMatchesReprobe(uint64_t pos, TSX::tsx_key_t& key, TSX::tsx_reprobe_t iReprobe)
    {

        TSX::tsx_keyval_t oKeyVal = getElement(pos);
        TSX::tsx_key_t oKey = this->extractKey(oKeyVal);

        // only check the reprobe part matches
        bool reprobeMatch = (oKey & m_iReprobeMask) == iReprobe;

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
     * given a keyval, constructs the value part
     * @param keyval
     * @return value part of keyval entry
     */
    TSX::tsx_val_t extractVal(TSX::tsx_keyval_t& keyval)
    {
        return keyval & (~m_iKeyMask);
    }

    /**
     *
     * @param pos
     * @return keyval of element at position pos in array
     */
    inline TSX::tsx_keyval_t getElement(uint64_t pos)
    {

        uint8_t iArray = pos % 4;
        pos = pos >> 4;

        uint8_t* pArray = m_pCounterArray[iArray];

        // we want to get the iPosition-th entry

        uint64_t iStartBit = pos * m_iKeyValBits;
        uint64_t iStartByte = (iStartBit >> 3); // / 8
        uint8_t iStartOffset = (iStartBit & 0x3); // mod 8

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt oReturn = UBigInt::createFromMemory(pArray, iStartByte, iStartOffset, m_iKeyValBits);

        std::cerr << "keyval fetched: " << oReturn.to_string() << " from array " << std::to_string(iArray) << " in pos " << std::to_string(pos) << std::endl;

        return oReturn;

    }

    void storeElement(uint64_t pos, TSX::tsx_key_t& key, TSX::tsx_val_t& val)
    {

        TSX::tsx_keyval_t oKeyVal = (key << m_iStorageBits) | val;

        uint8_t iArray = pos % 4;
        pos = pos >> 4;

        uint8_t* pArray = m_pCounterArray[iArray];

        // we want to get the iPosition-th entry

        uint64_t iStartBit = pos * m_iKeyValBits;
        uint64_t iStartByte = (iStartBit >> 3); // / 8
        uint8_t iStartOffset = (iStartBit & 0x3); // mod 8

        oKeyVal.resize(m_iKeyValBits);

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt::storeIntoMemory(pArray, iStartByte, iStartOffset, oKeyVal, m_iKeyValBits);

        std::cerr << "keyval stored: " << oKeyVal.to_string() << " into array " << std::to_string(iArray) << " in pos " << std::to_string(pos) << std::endl;

    }


    /**
     *
     * @param iPosition position in array
     * @param key key to look for
     * @param reprobes amount reprobes used
     * @param bKeyIsValue if true, the key part belongs to value
     * @return
     */
    bool incrementElement(uint64_t iPosition, TSX::tsx_key_t& key, TSX::tsx_reprobe_t reprobes, bool bKeyIsValue)
    {

        TSX::tsx_keyval_t keyval = getElement(iPosition);
        TSX::tsx_val_t value = this->extractVal(keyval);

        if (~value == 0)
        {
            // an overflow is going to happen!
            value = 1;

            // maybe we can prevent the overflow?
            if (bKeyIsValue)
            {

                TSX::tsx_key_t keyVal = this->extractKeyNoReprobe(keyval);

                if (~keyVal == 0)
                {
                    // overflow in key part
                    keyVal = 0;

                    keyVal = (keyVal << m_iReprobeBits) | reprobes;

                    storeElement(iPosition, keyVal, value);

                    return true;

                } else {
                    key += 1;

                    keyVal = (keyVal << m_iReprobeBits) | reprobes;

                    storeElement(iPosition, keyVal, value);
                    return false;
                }

            }

            // overflow occurred!
            return true;
        } else {

            value += 1;

            // key does not change here
            storeElement(iPosition, key, value);
            return false;

        }
    }

    bool handleOverflow(uint64_t iPos, TSX::tsx_key_t& iKey, TSX::tsx_reprobe_t iReprobe)
    {
        bool bHandled = false;

        iReprobe += 1;
        tsx_reprobe_t iPerformedReprobes = 1;

        while (iPerformedReprobes < m_iMaxReprobes)
        {

            ++iPerformedReprobes;
            uint64_t iPos = this->getPosition(iKey, iReprobe);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                // there can not be an overflow ;)
                this->incrementElement(iPos, iKey, iPerformedReprobes, true);

                return true;

            } else {

                bool bMatchesKey = positionMatchesReprobe(iPos, iKey, iPerformedReprobes);

                bool bOverflow = this->incrementElement(iPos, iKey, iPerformedReprobes, true);

                if (!bOverflow)
                {
                    return true;
                }

            }

        }

        return false;

    }

    bool insert_kmer(TSX::tsx_kmer_t& iKmer, uint8_t iCount = 1)
    {

        bool bInserted = false;
        TSX::tsx_reprobe_t iReprobe = 1;

        TSX::tsx_key_t iBaseKey = m_pHashingFunction->apply(iKmer);

        std::cerr << "new iKmer " << iKmer.to_string() << std::endl;
        std::cerr << "new basekey " << iBaseKey.to_string() << std::endl;

        while (!bInserted)
        {

            // get possible position
            uint64_t iPos = this->getPosition(iKmer, iReprobe);

            // does this position match to key?
            bool bEmpty = positionEmpty(iPos);

            if (bEmpty)
            {

                TSX::tsx_key_t invReprobe = ~m_iReprobeMask;
                std::cerr << "reprobe " << invReprobe.to_string() << std::endl;

                TSX::tsx_key_t key = (iBaseKey & invReprobe);
                std::cerr << "key and inv " << key.to_string() << std::endl;
                key = key | iReprobe;

                std::cerr << "create keyval: " << key.to_string() << " for key " << key.to_string() << " and reprobes: " << iReprobe << std::endl;

                this->incrementElement(iPos, key, iReprobe, false);

                bInserted = true;
                break;

            } else {

                bool bMatchesKey = positionMatchesKey(iPos, iBaseKey, iReprobe);

                if (bMatchesKey)
                {

                    bool bOverflow = this->incrementElement(iPos, iBaseKey, iReprobe, 1);

                    if (bOverflow)
                        handleOverflow(iPos, iBaseKey, iReprobe);

                    bInserted = true;
                    break;

                } else {
                    ++iReprobe;
                    continue;
                }


            }


        }

    }


    /*
     *
     * Variables
     *
     */


    uint8_t m_iCounterArrays = -1;
    uint8_t** m_pCounterArray;

    TSX::tsx_key_t m_iKeyMask = 0;
    TSX::tsx_val_t m_iValMask = 0;
    TSX::tsx_key_t m_iReprobeMask = 0;
//    tsx_kmer_t m_iReprobeMask = 0;


    /*
     *
     * Variables needed for the reading/writing of key,value pairs
     *
     */
    const uint8_t m_iLengthExp;
    const uint32_t m_iLength;
    const uint32_t m_iStorageBits;
    const uint32_t m_iK;
    const uint32_t m_iKeyValBits;

    const uint32_t m_iReprobeBits;
    const uint32_t m_iMaxReprobes;

    const UBigInt m_iMapSize;

    IBijectiveFunction* m_pHashingFunction;


};


#endif //TSXCOUNT_TSXHASHMAP_H
