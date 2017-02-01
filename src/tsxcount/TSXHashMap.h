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

#include <fastxtools/SequenceUtils.h>
#include <fastxtools/Sequence.h>
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
              m_iKeyValBits( 2 * iK - iExpLength + iStorage )
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
        m_iKeyMask = (~m_iKeyMask >> iExpLength);
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

    bool addKmer(tsx_kmer_t kmer)
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

    inline uint64_t getPosition(tsx_kmer_t iKmer, tsx_reprobe_t iReprobe)
    {

        // Uk is still in range [0, 4^k -1]
        tsx_kmer_t Uk = m_pHashingFunction->apply(iKmer) + iReprobe;

        // Uk.mod2 now is in range [0, M-1] = [0, 2^l -1]
        uint64_t iPos = Uk.mod2( m_iLengthExp ).toUInt();

        return iPos;
    }

    inline TSX::tsx_keyval_t getElement(uint64_t iPosition)
    {

        uint8_t iArray = iPosition % 4;
        iPosition = iPosition >> 4;

        uint8_t* pArray = m_pCounterArray[iArray];

        // we want to get the iPosition-th entry

        uint64_t iStartBit = iPosition * m_iKeyValBits;
        uint64_t iStartByte = (iStartBit >> 3); // / 8
        uint8_t iStartOffset = (iStartBit & 0x3); // mod 8

        //UBigInt oReturn(m_iKeyValBits);
        UBigInt oReturn = UBigInt::createFromMemory(pArray, iStartByte, iStartOffset, m_iKeyValBits);

        return oReturn;

    }

    TSX::tsx_key_t extractKey(TSX::tsx_keyval_t& keyval)
    {

    }

    bool positionMatchesKey(uint64_t pos, tsx_key_t& key)
    {

        TSX::tsx_keyval_t oKeyVal = getElement(pos);
        TSX::tsx_key_t oKey = this->extractKey(oKeyVal);

        return (oKey == key);

    }

    bool insert_kmer(tsx_kmer_t iKmer, uint8_t iCount = 1)
    {

        bool bInserted = false;
        tsx_reprobe_t iReprobe = 0;

        tsx_key_t iKey = m_pHashingFunction->apply(iKmer);


        while (!bInserted)
        {

            // get possible position
            uint64_t iPos = this->getPosition(iKey, iReprobe);

            // does this position match to key?
            bool bMatchesKey = positionMatchesKey(iPos, iKey);

            if (bMatchesKey) {
                // insert new value


                bInserted = true;
            }

            ++iReprobe;
        }

    }


    /*
     *
     * Variables
     *
     */


    uint8_t m_iCounterArrays = -1;
    uint8_t** m_pCounterArray;

    tsx_key_t m_iKeyMask = 0;
    tsx_val_t m_iValMask = 0;
    tsx_key_t m_iReprobeMask = 0;
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

    const UBigInt m_iMapSize;

    IBijectiveFunction* m_pHashingFunction;


};


#endif //TSXCOUNT_TSXHASHMAP_H
