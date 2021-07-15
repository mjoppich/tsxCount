//
// Created by joppich on 6/28/16.
//

#ifndef TSXCOUNT_BIGINT_H
#define TSXCOUNT_BIGINT_H


#include <cstddef>
#include <cstdint>
#include <string>
#include <cstring>

#include <cmath>

#include <algorithm>
#include <exception>
#include <sstream>
#include <functional>
#include <iostream>
#include <bitset>
#include <stdlib.h>

#include "MemoryPool.h"

#ifndef BITSTOFIELDS
#define BITSTOFIELDS(x) (x >> 3)
#endif

#include <src/tsxcount/commons.h>


class UBigInt;

class UBigIntGeneralException: public std::exception
{
public:
    UBigIntGeneralException(UBigInt* pBigInt)
            : std::exception()
    {
        m_pBigInt = pBigInt;
    }


    virtual const char* what() const throw()
    {
        return this->getExceptText();
    }

    virtual const char* getExceptText() const = 0;

protected:

    UBigInt* m_pBigInt = NULL;

};

class UBigIntNotLargeEnough: public UBigIntGeneralException
{
public:
    UBigIntNotLargeEnough(UBigInt* pBigInt, uint64_t iBitsQueried)
            : UBigIntGeneralException(pBigInt)
    {
        m_pBigInt = pBigInt;
        m_iBitsQueried = iBitsQueried;
    }

    virtual const char* getExceptText() const;

protected:

    uint64_t m_iBitsQueried;
};

class UBigIntNotImplemented: public UBigIntGeneralException
{
public:
    UBigIntNotImplemented(UBigInt* pBigInt, std::string sMethod)
            : UBigIntGeneralException(pBigInt)
    {
        m_pBigInt = pBigInt;
        m_sMethod = sMethod;
    }

    virtual const char* getExceptText() const;

protected:

    std::string m_sMethod;
};

/**
 * \class UBigInt
 *
 *
 * \brief Big Integer class for TSX Count
 *
 * This class implements Big Integers (needing more than 64bit)
 *
 * \note Performance check this
 *
 * \author Markus Joppich
 *
 */
class UBigInt {

public:

    uint32_t getBitCount()
    {
        return this->m_iBits;
    }

    static UBigInt createFromBitShift(uint32_t iBits, uint32_t iPosition, MemoryPool<FIELDTYPE>* pPool)
    {

        UBigInt oRet = UBigInt(iBits, true, pPool);

        oRet = 1;
        oRet = oRet << iPosition;

        return oRet;

    }

    void setPool(MemoryPool<FIELDTYPE>* pPool)
    {
        this->m_pPool = pPool;
    }


    /**
     *
     * @return string starting with highest bit to lowest bit (0 = unset, 1 = set)
     */
    std::string to_string() const
    {
        std::stringstream oSS;

        for (uint32_t i = 0; i < m_iBits; ++i)
        {
            uint8_t iValue = this->getBit(i);
            oSS << (char) (48+iValue);
            //std::cerr << (char) (48+iValue) << std::endl;
        }


        std::string sBitSeq = oSS.str();

        std::reverse(sBitSeq.begin(), sBitSeq.end());

        return sBitSeq;
    }

    std::string to_debug() const {
        std::stringstream oSS;

        uint32_t iBits = this->m_iFields * this->m_iFieldSize;

        for (uint32_t i = 0; i < iBits; ++i)
        {
            uint8_t iValue = this->getBit(i);
            oSS << (char) (48+iValue);
            //std::cerr << (char) (48+iValue) << std::endl;
        }


        std::string sBitSeq = oSS.str();

        std::reverse(sBitSeq.begin(), sBitSeq.end());

        return sBitSeq;
    }

    void print_string()
    {
        std::cout << this->to_string() << std::endl;
    }

    /**
     *
     * @param pPosition start position of array
     * @param iIdx Index in Array
     * @param iOffsetBits starting bit in Array
     * @param iLength of memory block in bits
     * @return
     */
    static UBigInt createFromMemory(FIELDTYPE* pPosition, uint64_t iIdx, uint8_t iOffsetBits, uint32_t iLength, MemoryPool<FIELDTYPE>* pPool)
    {

        FIELDTYPE* pPos = pPosition + iIdx;
        UBigInt oReturn = UBigInt(iLength, true, pPool);
        oReturn.copy_content_bits(pPos, iOffsetBits, iLength);

        return oReturn;
    }

    /**
     *
     * @param pPosition array start
     * @param iIdx array offset (bytes)
     * @param iOffsetBits offset in byte
     * @param oValue @UBigInt to store in array
     * @return
     */
    inline
    static void storeIntoMemory(FIELDTYPE* pPosition, uint64_t iIdx, uint8_t iOffsetBits, UBigInt& oValue, size_t iBits = -1)
    {

        FIELDTYPE* pPos = pPosition + iIdx;

        if (iBits == -1)
            iBits = oValue.getBitCount();

        oValue.copy_content_to_array(pPos, iOffsetBits, iBits);

    }

    UBigInt(uint64_t iBits, uint8_t* pMemPos, uint8_t iOffset, MemoryPool<FIELDTYPE>* pPool)
    : UBigInt(iBits, true, pPool)
    {
        uint32_t iFieldsNeeded = std::ceil(iBits / this->m_iFieldSize);

        UBigInt oIntermediate = UBigInt(iFieldsNeeded * this->m_iFieldSize, false, pPool);
        oIntermediate.copy_content_from_le((FIELDTYPE*) pMemPos, iBits);

    }

    UBigInt()
            : UBigInt(0, true, NULL)
    {

    }

    /**
     * \brief default constructor
     */
    UBigInt(MemoryPool<FIELDTYPE>* pPool)
    : UBigInt(64, true, pPool)
    {

    }

    UBigInt(std::string sInput, MemoryPool<FIELDTYPE>* pPool)
    : UBigInt(sInput.size(), true, pPool)
    {
        for (size_t i = 0; i < sInput.size(); ++i)
        {
            size_t iIdx = sInput.size() - 1 -i;

            this->setBit(i, sInput.at(iIdx) == '1' ? 1 : 0);
        }

        std::cerr << sInput << std::endl;
        std::cerr << this->to_string() << std::endl;
    }


    /** \brief UBigInt constructor
      * \param iBits number of bits required.
      * \param bZero if set true, the object will have value 0
      *
      * This implements a new UBigInt constructor. The Array can be directly initialized with 0
      */
    UBigInt(uint64_t iBits, bool bZero, MemoryPool<FIELDTYPE>* pPool)
    {
        initialize(iBits, pPool);
    }

    UBigInt(UBigInt & oOther )
     : UBigInt(oOther.m_iBits, true, oOther.m_pPool)
    {
        copy_content(oOther);
    }

    UBigInt(UBigInt & oOther, uint32_t iBits )
            : UBigInt(iBits, true, oOther.m_pPool)
    {
        copy_content(oOther);
    }

    UBigInt(UBigInt & oOther, MemoryPool<FIELDTYPE>* pPool)
            : UBigInt(oOther.m_iBits, true, pPool)
    {
        copy_content(oOther);
    }

    UBigInt(const UBigInt & oOther )
            : UBigInt(oOther.m_iBits, true, oOther.m_pPool)
    {
        copy_content(oOther);
    }

    ~UBigInt()
    {
        this->reset();
    }

    /**
     *
     * @brief resizes this UBigInt such that it can hold iBits
     * @param iBits number of bits this element can hold
     */
    void resize(uint32_t iBits)
    {

        // Nothing to be done?
        if (iBits == m_iBits)
            return;

        FIELDTYPE* pOldArray = m_pArray;
        MemLoc oOldMemLoc = m_oAlloc;

        uint32_t iOldBits = m_iBits;

        m_pArray = NULL;

        this->initialize(iBits);

        if (pOldArray != NULL)
        {

            uint64_t bits = std::min(iOldBits, m_iBits);
            div_t fields = div( bits, m_iFieldSize );

            // copy old data over
            for (uint32_t i = 0; i < fields.quot; ++i)
            {
                m_pArray[i] = pOldArray[i];
            }

            FIELDTYPE oInter = (pOldArray[fields.quot] << (m_iFieldSize-fields.rem));

            m_pArray[ fields.quot ] = oInter >> (m_iFieldSize-fields.rem);

        }

        this->m_pPool->free(oOldMemLoc);

    }
    /**
     *
     * @param value
     * @return big int with x bits and value copied
     */
    static UBigInt fromUint8(uint8_t value, MemoryPool<FIELDTYPE>* pPool)
    {

        FIELDTYPE oVal = value;

        UBigInt oRet(8, true, pPool);
        oRet.copy_content_to_array(&oVal, 0, 8);

        return oRet;
    }

    /**
     * TODO does this work?
     *
     * @param value
     * @return
     */
    static UBigInt fromUint16(uint16_t value, MemoryPool<FIELDTYPE>* pPool)
    {
        FIELDTYPE oVal = value;

        UBigInt oRet(8, true, pPool);
        oRet.copy_content_to_array( &oVal, 0, 16);

        return oRet;
    }




    UBigInt(uint64_t iValue, MemoryPool<FIELDTYPE>* pPool)
    : UBigInt(0, true, pPool)
    {

        uint32_t bits = UBigInt::log_2(iValue)+1;
        this->initialize(bits);
        this->copy_content_bits((FIELDTYPE*)&iValue, 0, bits);

    }

    static inline uint32_t log_2(const uint64_t x)
    {
        if(x == 0)
            return 0;

        return (63 - __builtin_clzll (x));
    }

    UBigInt& operator = (uint64_t iValue)
    {

        this->reset();
        this->initialize(std::max((uint32_t)64, m_iBits));
        this->copy_content_from_be((FIELDTYPE*)&iValue, 64);

        return *this;

    }

    UBigInt& operator= (const UBigInt & oOther)
    {

        this->reset();
        this->setPool(oOther.m_pPool);
        this->initialize(oOther.m_iBits);
        this->copy_content(oOther);

        return *this;

    }

    UBigInt operator- (const uint64_t iLeft)
    {
        UBigInt oLeft(iLeft, this->m_pPool);

        return operator-(oLeft);
    }

    UBigInt operator- (const UBigInt & oOther)
    {
        UBigInt oRet(this->self() );
        oRet.sub(oOther);
        return oRet;
    }

    UBigInt & operator-= (const UBigInt & oOther)
    {
        this->sub(oOther);
        return *this;
    }

    UBigInt operator+ (const uint64_t iLeft)
    {
        UBigInt oLeft(iLeft, this->m_pPool);

        return operator+(oLeft);
    }

    UBigInt operator++()
    {
        this->add(UBigInt(1, this->m_pPool)); // TODO add constant 1
        return *this;
    }

    UBigInt operator++(int dummy)
    {
        UBigInt ret(this->self());

        this->add(UBigInt(1, this->m_pPool)); // TODO add constant 1
        return ret;
    }

    UBigInt operator+ (const UBigInt & oOther)
    {
        UBigInt oRet(this->self() );
        oRet.add(oOther);
        return oRet;
    }

    UBigInt & operator+= (const UBigInt & oOther)
    {
        this->add(oOther);
        return *this;
    }

    UBigInt operator << (unsigned int iShift)
    {
        UBigInt oRet( this->self() );
        oRet.shiftLeft( iShift );
        return oRet;
    }

    UBigInt operator >> (unsigned int iShift)
    {
        UBigInt oRet(this->self() );
        oRet.shiftRight( iShift );
        return oRet;
    }

    UBigInt operator ^ (const UBigInt& oOther)
    {
        UBigInt oRet(this->self());

        oRet.bitXor(oOther);

        return oRet;


    }

    UBigInt operator & (const UBigInt & oOther)
    {
        UBigInt oRet( this->self() );
        oRet.bitAnd(oOther);

        return oRet;
    }

    UBigInt operator | (const UBigInt & oOther)
    {
        UBigInt oRet(this->self() );
        oRet.bitOr(oOther);

        return oRet;
    }


    UBigInt operator ~ ()
    {
        UBigInt oRet(this->self() );
        oRet.bitCompl();

        return oRet;
    }


    bool isZero()
    {
        div_t procs = div(this->m_iBits, m_iFieldSize);
        size_t checkIdx = 0;

        if (procs.quot > 0)
        {
            for (checkIdx = 0; checkIdx <= procs.quot-1; ++checkIdx)
            {
                if (this->m_pArray[checkIdx] != 0)
                    return false;
            }
        }

        uint32_t iShift = m_iFieldSize - procs.rem;

        FIELDTYPE thisRemain = this->m_pArray[checkIdx] << iShift;

        bool result = (thisRemain == 0);

        return result;
    }

    bool isEqual(const UBigInt* pLeft, const UBigInt* pRight, uint32_t iBits)
    {
        div_t procs = div(iBits, m_iFieldSize);
        size_t checkIdx = 0;

        if (procs.quot > 0)
        {
            for (checkIdx = 0; checkIdx <= procs.quot-1; ++checkIdx)
            {
                if (pLeft->m_pArray[checkIdx] != pRight->m_pArray[checkIdx])
                    return false;
            }
        }

        size_t iShift = m_iFieldSize - procs.rem;

        FIELDTYPE thisRemain = pLeft->m_pArray[checkIdx] << iShift;
        FIELDTYPE otherRemain = pRight->m_pArray[checkIdx] << iShift;

        bool result = (thisRemain == otherRemain);

        return result;
    }

    void transferSelf(UBigInt* pLoc)
    {
        ::memcpy(pLoc, this, sizeof(UBigInt));

        this->m_pArray = NULL;
        //std::cout << "transfer complete" << std::endl;

        this->initialize(m_iBits);

        //std::cout << "transfer complete reinit" << std::endl;
    }

    bool restZero(const UBigInt* pElement, uint32_t iOffsetBits)
    {

        div_t procs = div(iOffsetBits, pElement->m_iFieldSize);

        FIELDTYPE thisRemain = pElement->m_pArray[procs.quot] >> procs.rem;

        if (thisRemain != 0)
            return false;

        for (uint32_t i = procs.quot; i < pElement->m_iFields-1; ++i)
        {

            if (pElement->m_pArray[i] != 0)
                return false;

        }

        if (pElement->m_iFields > 1)
        {
            FIELDTYPE value = pElement->m_pArray[pElement->m_iFields-1] & ~m_iUnusedBitsMask;

            if (value != 0)
                return false;
        }

        return true;

    }

    bool operator == (const UBigInt & oOther)
    {

        // rational: if this and other have same bits: all bits must equal

        if (this->m_iBits == oOther.m_iBits)
        {
            return this->isEqual(this, &oOther, this->m_iBits);

        } else {

            if (this->m_iBits < oOther.m_iBits)
            {

                bool commonEqual = this->isEqual(this, &oOther, this->m_iBits);

                bool restZero = this->restZero(&oOther, this->m_iBits);

                return commonEqual && restZero;

            } else {
                bool commonEqual = this->isEqual(this, &oOther, oOther.m_iBits);

                bool restZero = this->restZero(this, oOther.m_iBits);

                return commonEqual && restZero;
            }

        }

        return false;

    }

    uint64_t toUInt(bool verbose=false)
    {
        if (m_iBits > 64)
        {
            throw ("Cannot represent BigInt in 64bits!");
        }

        if (verbose)
        {
            this->print_string();
        }

        if (m_iBits <= m_iFieldSize)
        {

            uint32_t iDiff = m_iFieldSize-m_iBits;

            FIELDTYPE iField = m_pArray[0];

            iField = (iField << iDiff);
            iField = iField >> iDiff;

            return iField;

        } else {

            // iBits > iFieldSize
            // assemble this uint64_t from multiple fields

            size_t iNeeded = m_iFields;
            uint64_t iReturn = 0;

            for (uint8_t i = 0; i < iNeeded; ++i)
            {
                //std::cerr << std::bitset<8>(m_pArray[i]) << std::endl;
                iReturn |= m_pArray[i] << (i*m_iFieldSize);

            }

            return iReturn;
        }

        // unreachable :)
        return 0;
    }

    uint64_t operator % (const uint64_t iRight)
    {

        throw "not yet implemented";
        return 0;
    }

    /**
     *
     * \brief calculates mod 2^BitPos
     * \param iBitPos the position for which bit the value is set to 1. All other values are 0
     *
     */
    UBigInt mod2(uint64_t iBitPos)
    {

        // return this & ((1 << iBitPos) -1)

        UBigInt oRet(iBitPos, true, this->m_pPool);
        oRet.copy_content_bits( this->m_pArray, 0, iBitPos);
        return oRet;

    }

    bool operator != (const UBigInt & oOther)
    {
        return !( operator==(oOther) );
    }



    struct sBitPosition
    {
        const uint64_t iField;
        const uint8_t iPosInField;

        sBitPosition(uint64_t iField, uint8_t iPos)
        : iField(iField), iPosInField(iPos)
        {

        }
    };

    sBitPosition getBitPosition(uint32_t iBit) const
    {

        uint64_t iField = iBit / m_iFieldSize;
        uint8_t iPos = iBit - m_iFieldSize * iField;


        return sBitPosition(iField, iPos);
    }

    uint8_t getBit(uint32_t iBitPosition) const
    {
        sBitPosition oPos = this->getBitPosition(iBitPosition);

        FIELDTYPE iFieldValue = m_pArray[oPos.iField];

        iFieldValue = (iFieldValue >> oPos.iPosInField) & 1;

        return (uint8_t) iFieldValue;
    }

    void preloadBit(uint32_t iBitPosition) const
    {
        sBitPosition oPos = this->getBitPosition(iBitPosition);

        FIELDTYPE* pFieldValue = m_pArray+ oPos.iField;
        __atomic_compare_exchange(pFieldValue, pFieldValue, pFieldValue, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
    }

    void setBit(uint64_t iBit, uint8_t iValue)
    {

        sBitPosition oPos = getBitPosition(iBit);

        this->setBitDirect(oPos.iField, oPos.iPosInField, iValue);
    }

    void toggleBit(uint64_t iBit)
    {
        sBitPosition oPos = getBitPosition(iBit);

        FIELDTYPE iFieldValue = m_pArray[oPos.iField];
        iFieldValue ^= (1 << oPos.iPosInField);
        m_pArray[oPos.iField] = iFieldValue;

    }


    uint32_t sumBits() {

        uint32_t iCount = 0;

        for (uint32_t i = 0; i < m_iFields; ++i)
        {

            uint32_t iValue = 0;

            if (m_iFieldSize > 32)
            {

                int32_t iBitsRemaining = m_iFieldSize;
                uint32_t iChunk = 0;

                FIELDTYPE oFieldValue = m_pArray[i];
                iChunk = oFieldValue;

                while (iBitsRemaining > 0)
                {
                    iValue += __builtin_popcount( iChunk );

                    oFieldValue = oFieldValue >> 32;
                    iBitsRemaining -= 32;

                    iChunk = oFieldValue;

                }

            } else {
                iValue = __builtin_popcount( m_pArray[i] );
            }


            iCount += iValue;


        }

        return iCount;

    }

    void bitXor(const UBigInt & oOther)
    {

        for (uint32_t i = 0; i < m_iFields; ++i)
        {

            if ( m_pArray[i] == oOther.m_pArray[i])
            {
                m_pArray[i] = 0;
                continue;
            }

            m_pArray[i] = m_pArray[i] ^ oOther.m_pArray[i];


        }
        
    }

    std::vector<uint64_t> onePositions() {

        std::vector<uint64_t> oRet;

        for (uint32_t i = 0; i < m_iFields; ++i)
        {

            if (m_pArray[i] == 0)
            {
                continue;
            }

            for (uint8_t j = 0; j < m_iFieldSize; ++j)
            {
                FIELDTYPE oPosValue = ((m_pArray[i] >> j) & 1);

                if (oPosValue != 0)
                {
                    uint64_t iCurPos = i*m_iFields + j;
                    oRet.push_back(iCurPos);
                }
            }
        }

        return oRet;

    }

    static UBigInt fromString(std::string sInput, MemoryPool<FIELDTYPE>* pPool) {

        UBigInt oBigInt(sInput.size(), true, pPool);

        for (size_t i = 0; i < sInput.size(); ++i)
        {
            size_t iBitPos = sInput.size()-1-i;

            if (sInput.at(i) == '1')
            {
                oBigInt.setBit(iBitPos, 1);
            } else {
                oBigInt.setBit(iBitPos, 0);
            }

        }

        return oBigInt;

    }

    static void shift_array_left(FIELDTYPE* pArray, uint16_t iFieldSize, uint16_t iFields, uint16_t iShift)
    {

        if (iShift >= iFieldSize)
        {

            uint16_t iFields = iShift / iFieldSize;

            uint32_t i;
            for (i = 0; i < iFields-iFields; ++i)
            {
                pArray[iFields-1 -i] = pArray[iFields-1-i-iFields];
            }

            for (i; i < iFields; ++i)
            {
                pArray[iFields-1 -i] = 0;
            }
        }

        iShift = iShift % iFieldSize;
        uint32_t iCounterShift = iFieldSize - iShift;

        FIELDTYPE iCurField, iNextFieldPart;

        for (uint32_t i = 0; i < iFields-1; ++i)
        {
            iCurField = (pArray[iFields-1-i] << iShift);
            iNextFieldPart = (pArray[iFields-1-i-1] >> iCounterShift);
            iCurField = iCurField | iNextFieldPart;

            pArray[iFields-1 -i] = iCurField;
        }

        pArray[0] = (pArray[0] << iShift);

    }


    void shiftLeft( uint16_t iShift)
    {

        if (iShift >= m_iFieldSize)
        {

            uint16_t iFields = iShift / m_iFieldSize;

            uint32_t i;
            for (i = 0; i < m_iFields-iFields; ++i)
            {
                m_pArray[m_iFields-1 -i] = m_pArray[m_iFields-1-i-iFields];
            }

            for (i; i < m_iFields; ++i)
            {
                m_pArray[m_iFields-1 -i] = 0;
            }
        }

        iShift = iShift % m_iFieldSize;
        uint32_t iCounterShift = m_iFieldSize - iShift;

        FIELDTYPE iCurField, iNextFieldPart;

        for (uint32_t i = 0; i < m_iFields-1; ++i)
        {
            iCurField = (m_pArray[m_iFields-1-i] << iShift);
            iNextFieldPart = (m_pArray[m_iFields-1-i-1] >> iCounterShift);
            iCurField = iCurField | iNextFieldPart;

            m_pArray[m_iFields-1 -i] = iCurField;
        }

        m_pArray[0] = (m_pArray[0] << iShift);

    }

    void shiftRight( uint32_t iShift)
    {

        if (iShift >= m_iFieldSize)
        {
            uint64_t iFields = iShift / m_iFieldSize;

            uint32_t i;
            for (i = 0; i < m_iFields - iFields; ++i)
            {
                m_pArray[i] = m_pArray[i + iFields];
            }

            for (i; i < m_iFields; ++i)
            {
                m_pArray[i] = 0;
            }
        }

        iShift = iShift % m_iFieldSize;
        uint32_t iCounterShift = m_iFieldSize - iShift;

        for (uint32_t i = 0; i < m_iFields-1; ++i)
        {

            FIELDTYPE iNewRight = (m_pArray[i] >> iShift);
            FIELDTYPE iNewLeft = (m_pArray[i+1] << iCounterShift);

            FIELDTYPE iNewValue =  iNewLeft | iNewRight;
            m_pArray[i] = iNewValue;

        }

        m_pArray[m_iFields-1] = (m_pArray[m_iFields-1] >> iShift);
    }

    /**
     *
     * @param pSrc destination array
     * @param iBitStart offset in first field of pSrc
     * @param iBitCount how many bits to copy
     */
    void copy_content_to_array(FIELDTYPE* pDest, uint32_t iBitStart, uint32_t iBitCount)
    {

        // how many bits can I copy into first, not fully occupied field?
        uint8_t iOffset = iBitStart % m_iFieldSize;

        // TODO copy first field
        FIELDTYPE iThisVal = m_pArray[0] << iOffset;
        FIELDTYPE iOldVal = pDest[0] << (m_iFieldSize - iOffset);
        iOldVal = iOldVal >> (m_iFieldSize - iOffset);

        pDest[0] = iThisVal | iOldVal;

        // how many bits remain to be copied?
        uint32_t iBitsRemaining = iBitCount - (m_iFieldSize - iOffset);
        // how many full fields are this?
        div_t fields = div(iBitsRemaining, m_iFieldSize);

        // copy full fields over
        FIELDTYPE iPartRight;
        FIELDTYPE iPartLeft;

        uint32_t i = 0;
        for ( i = 0; i < fields.quot; ++i)
        {
            // TODO copy full fields
            iPartRight = m_pArray[i] >> (m_iFieldSize-iOffset);
            iPartLeft = m_pArray[i+1] << iOffset;
            iThisVal = iPartLeft | iPartRight;

            pDest[i+1] = iThisVal;

            iBitsRemaining -= m_iFieldSize;
        }

        if (iBitsRemaining > 0)
        {
            if (iBitsRemaining <= iOffset)
            {

                // we do not need a rest from the next field
                iPartRight = m_pArray[fields.quot] << (m_iFieldSize-iOffset);
                iPartRight = iPartRight << (iOffset-iBitsRemaining);
                iPartRight = iPartRight >> (iOffset-iBitsRemaining);
                iPartRight = iPartRight >> m_iFieldSize-iOffset;

                iOldVal = pDest[fields.quot+1];
                iOldVal = iOldVal >> iBitsRemaining;
                iOldVal = iOldVal << iBitsRemaining;

                pDest[fields.quot+1] = iOldVal | iPartRight;

            } else {

                // we need a rest from the next field
                iPartRight = m_pArray[fields.quot] >> (m_iFieldSize-iOffset);

                uint32_t iRemaining = iBitsRemaining - iOffset;

                // need the rightmost iRemaining bits
                iPartLeft = m_pArray[fields.quot+1] << (m_iFieldSize-iRemaining); // clears all but iRemaining bits
                iPartLeft = iPartLeft >> (m_iFieldSize-iRemaining);
                iPartLeft = iPartLeft << iOffset;

                iThisVal = iPartLeft | iPartRight;

                iOldVal = pDest[fields.quot+1];
                iOldVal = iOldVal >> iBitsRemaining;
                iOldVal = iOldVal << iBitsRemaining;

                pDest[fields.quot+1] = iOldVal | iThisVal;

            }
            /*
            std::cout << "Having i: " << i << " and fields " << fields.quot << std::endl;
            std::cout << "Using the last bits: " << iBitsRemaining << std::endl;
            std::cout<< std::bitset<8>(pDest[fields.quot+1]) << " into field " << fields.quot+1 << " " << static_cast<void*>(pDest + fields.quot+1) << std::endl;
            */
        }



/*
        // how many bits remain?
        if (iBitsRemaining > 0)
        {

            FIELDTYPE iPartRight = m_pArray[fields.quot] >> (m_iFieldSize-iOffset);

            iBitsRemaining -= iOffset;

            FIELDTYPE iThisVal = 0;
            iThisVal |= m_pArray[fields.quot+1] << (m_iFieldSize-iBitsRemaining);
            iThisVal = iThisVal >> (m_iFieldSize-iBitsRemaining);

            if (this->m_iFields < fields.quot+1)
            {
                iThisVal |= m_pArray[fields.quot+1] << (m_iFieldSize-iBitsRemaining); // TODO can this field value exist?
            }

            //FIELDTYPE  iThisVal = m_pArray[fields.quot+1] << (m_iFieldSize-iBitsRemaining);
            //iThisVal = iThisVal >> (m_iFieldSize-iBitsRemaining);

            FIELDTYPE iOldVal = pDest[fields.quot+1] >> iBitsRemaining;
            iOldVal = iOldVal << iBitsRemaining;

            // TODO copy into last field
            pDest[fields.quot+1] = iOldVal | iThisVal;
        }
*/
    }

    /**
 *
 * @brief copies x bits from another instance to this instance
 * @param pSrc the source field
 * @param iBitStart start in pSrc field. Next bit will be first bit in this instance
 * @param iBitCount bitcount
 *
 */
    void copy_content_bits(FIELDTYPE* pSrc, uint32_t iBitStart, uint32_t iBitCount)
    {

        uint8_t iOffset = iBitStart % m_iFieldSize;

        if (iOffset == 0)
        {

            this->copy_content_from_be(pSrc, iBitCount);

        } else {

            uint32_t iRemainingBits = iBitCount;
            div_t fields = div(iRemainingBits, m_iFieldSize);

            //
            FIELDTYPE iPartOne;
            FIELDTYPE iPartTwo;
            FIELDTYPE iFieldValue;
            for (uint32_t i = 0;  i < fields.quot; ++i)
            {

                iPartOne = pSrc[i] >> iOffset;
                iPartTwo = pSrc[i+1] << (m_iFieldSize-iOffset);
                iFieldValue = iPartOne | iPartTwo;

                /*
                std::cout << "i   " << std::bitset<8>(pSrc[i]) << " " << static_cast<void*>(pSrc + i) << std::endl;
                std::cout << "i+1 " << std::bitset<8>(pSrc[i+1]) << " " << static_cast<void*>(pSrc + i+1) << std::endl;
                std::cout<<std::bitset<8>(iFieldValue)<< " from field " << i << " " << static_cast<void*>(pSrc + i) << std::endl;
                 */

                m_pArray[i] = iPartOne | iPartTwo;

                iRemainingBits -= m_iFieldSize;

            }

            if (iRemainingBits > 0)
            {
                uint8_t iUnusedBits = (iBitStart+iBitCount) % m_iFieldSize;
                uint8_t iMissingBits = iBitCount % m_iFieldSize;

                /*
                if (iMissingBits != iRemainingBits)
                    std::cerr << "offset " << std::to_string(iOffset) << " and field size-remain inqueal " << std::to_string(m_iFieldSize-iRemainingBits) << std::endl;
                */

                FIELDTYPE iPart = pSrc[fields.quot];

                iPart = iPart << (m_iFieldSize-iOffset-iRemainingBits);
                iPart = iPart >> (m_iFieldSize-iOffset-iRemainingBits);

                iPart = iPart >> iOffset;

                //FIELDTYPE iPart = pSrc[fields.quot] >> (m_iFieldSize-iRemainingBits);
                //iPart = iPart << (m_iFieldSize-iRemainingBits);
                //iPart = iPart >> iOffset;

                m_pArray[fields.quot] = iPart;

            }

        }
    }

    /**
     *
     * \brief copies content from one bigint to other (starts at least)
     * \note copies least significant bytes only!
     * \param oOther Other big int to copy stuff from
     */
    void copy_content(const UBigInt & oOther)
    {

        uint32_t iMinBits = std::min(oOther.m_iBits, this->m_iBits);
        this->copy_content_bits(oOther.m_pArray, 0, iMinBits);

    }

    void bitAnd(const UBigInt & oOther)
    {
        const FIELDTYPE* pOtherField = oOther.m_pArray;
        const uint32_t iMaxField = oOther.m_iFields;

        this->applyToAllFields( [pOtherField, iMaxField] (FIELDTYPE iValue, uint32_t iField) {

            if (iField < iMaxField)
                return (FIELDTYPE) (iValue & pOtherField[iField]);

            return (FIELDTYPE) 0;
        } );

    }

    void bitOr(const UBigInt & oOther)
    {

        const FIELDTYPE* pOtherField = oOther.m_pArray;
        const uint32_t iMaxField = oOther.m_iFields;

        this->applyToAllFields( [pOtherField, iMaxField] (FIELDTYPE iValue, uint32_t iField) {

            return (iField < iMaxField) ? (iValue | pOtherField[iField]) : iValue;

        } );

    }



    uint32_t m_iBits = 0;
    uint32_t m_iFields = 0;

    FIELDTYPE* m_pArray = NULL;
    static const uint8_t m_iFieldSize = sizeof(FIELDTYPE) * 8;


protected:

    /**
     * \return const reference to this
     */
    const UBigInt & self() { return *this; }


    /**
     *
     * \brief initializes big int variable
     * \param bZero if set true, the big int value is set to 0
     *
     */
    void initialize(uint32_t iBits, MemoryPool<FIELDTYPE>* pPool=NULL)
    {

        if (pPool != NULL)
        {
            m_pPool = pPool;
        }

        if (iBits == 0)
            return;

        m_iBits = iBits;
        m_iFields = std::ceil( iBits / (8.0f*sizeof(FIELDTYPE)) );

        m_oAlloc = m_pPool->malloc(m_iFields);
        m_pArray = (FIELDTYPE*) m_oAlloc.addr;

        // NULLING the array
        for (uint32_t i = 0; i < m_iFields; ++i)
            m_pArray[i] = 0;


        m_iUnusedBitsMask = this->createUnusedBitsMask();

    }

    /**
     *
     * @return mask of bits in fieldtype that are unused (e.g. 11100000 if last 3 bits are unused)
     */
    FIELDTYPE createUnusedBitsMask()
    {
        FIELDTYPE iUnusedBitsMask = 0;

        uint32_t iRemainingBits = m_iFields * sizeof(FIELDTYPE) * 8 - m_iBits;

        const uint32_t stop = m_iFieldSize - iRemainingBits;

        for (uint32_t i = 0; i < stop; ++i)
            iUnusedBitsMask = (iUnusedBitsMask << 1) | 1;

        iUnusedBitsMask = ~iUnusedBitsMask;

        return iUnusedBitsMask;
    }

    /**
     * sets the posinfield-th bit in the iFields-th field to value iValue [0,1]
     * \param iField
     * \param iPosInField
     * \param iValue
     */
    void setBitDirect( uint64_t iField, uint8_t iPosInField, uint8_t iValue )
    {
        if (iValue == 0)
        {

            FIELDTYPE iFieldValue = m_pArray[iField];
            iFieldValue &= ~(1 << iPosInField);
            m_pArray[iField] = iFieldValue;

        } else {

            FIELDTYPE iFieldValue = m_pArray[iField];
            iFieldValue |= (1 << iPosInField);
            m_pArray[iField] = iFieldValue;

        }
    }








    void applyToAllFields(std::function<FIELDTYPE(FIELDTYPE, uint32_t)> oOp)
    {
        for (uint32_t i = 0; i < m_iFields; ++i)
        {
            m_pArray[i] = oOp(m_pArray[i], i);
        }
    }

    void sub(const UBigInt & oOther)
    {


        uint32_t iBitsUsed = this->m_iBits;

        this->bitCompl();
        this->add(oOther);
        this->bitCompl();


    }



    void add(const UBigInt & oOther)
    {

        FIELDTYPE iCarriage = 0;
        FIELDTYPE iMax = 0;

        if (oOther.m_iFields > this->m_iFields)
        {
            this->resize( oOther.m_iBits );
        }

        uint32_t iFields = std::min(oOther.m_iFields, this->m_iFields);


        uint64_t i = 0;
        for ( ; i < iFields; ++i)
        {

            iMax = std::max(m_pArray[i], oOther.m_pArray[i]);
            m_pArray[i] = m_pArray[i] + oOther.m_pArray[i] + iCarriage;

            iCarriage = 0;

            if (m_pArray[i] < iMax)
                iCarriage = 1;

        }

        if (iFields < this->m_iFields)
        {

            for (; i < m_iFields; ++i)
            {
                iMax = m_pArray[i];
                m_pArray[i] = m_pArray[i] + iCarriage;

                iCarriage = 0;

                if (m_pArray[i] < iMax)
                    iCarriage = 1;
            }
        }

        // cases:
        // carriage == 0 and no more bits used ==> do nothing
        // carriage == 0 and more bits used (but fits into fields) ==> add bit
        // carriage == 1  ==> add field + bit

        bool bMoreBitsUsed = (m_pArray[i-1] & m_iUnusedBitsMask) > 0;

        if ((iCarriage == 0) and (bMoreBitsUsed))
        {
            // TODO could this be more than one additional bit?
            this->resize(m_iBits+1);

        } else if (iCarriage == 1)
        {

            this->resize(m_iBits+1);

            this->setBit(m_iBits, 1);

            uint8_t iBit = m_iBits % m_iFieldSize;

        }

    }

    void bitCompl()
    {

        for (uint64_t i = 0; i < m_iFields-1; ++i)
        {
            m_pArray[i] = ~m_pArray[i];
        }

        FIELDTYPE iLastField = m_pArray[m_iFields-1];

        // 0b00vvvvvv & 0b00111111 --> 0b00vvvvvv -compl-> 0b11cccccc -&-> 0b00cccccc
        iLastField = (~(iLastField & ~m_iUnusedBitsMask)) & ~m_iUnusedBitsMask;

        m_pArray[m_iFields-1] = iLastField;

    }

    void reset()
    {

        if (m_pArray != NULL)
        {
            m_pPool->free(m_oAlloc);
            m_oAlloc.addr = NULL;
            m_pArray = NULL;
        }

    }

    /**
     *
     * \brief copies iLength bits from pSrc+iStart to pDest. It is assumed that the highest bit is at the lowest address.
     * \note this only applies if src is a contiguous array of values, e.g. little endian
     *
     */
    void copy_content_from_be(FIELDTYPE* pSrc, uint64_t iLength)
    {

        /*
        if (this->m_iBits < iLength)
            throw new UBigIntNotLargeEnough(this, iLength);
        */

        //const uint8_t iFields = std::ceil( iLength / m_iFieldSize );
        const div_t fields = div(iLength, m_iFieldSize);

        uint8_t i = 0;
        FIELDTYPE iInsert;
        for (i = 0; i < fields.quot; ++i)
        {
            iInsert = pSrc[i];
            m_pArray[i] = iInsert;
        }

        if (fields.rem > 0)
        {

            iInsert = pSrc[fields.quot] << (m_iFieldSize-fields.rem);
            m_pArray[ fields.quot ] = iInsert >> (m_iFieldSize-fields.rem);

        }

    }

    void reverse(FIELDTYPE* pField)
    {

        uint8_t iBits = sizeof(FIELDTYPE) * 8;

        FIELDTYPE oValue = 0;

        for (int i = 0; i < iBits; ++i)
        {

            oValue = oValue << 1;

            if ((pField[0] >> i) & 1 == 1)
            {
                oValue |= 1;
            }

        }

        *pField = oValue;
    }


    /**
 *
 * \brief copies iLength bits from pSrc+iStart to pDest. It is assumed that the highest bit is at the highest address.
 * \note this only applies if src is a contiguous array of values, e.g. little endian
 *
 */
    void copy_content_from_le(FIELDTYPE* pSrc, uint64_t iLength)
    {

        if (this->m_iBits < iLength)
            throw UBigIntNotLargeEnough(this, iLength);

        const uint8_t iFields = std::ceil( iLength / m_iFieldSize );

        for (uint8_t i = 0; i < iFields; ++i)
        {
            m_pArray[i] = pSrc[i];
        }

        // reverse fields
        for (uint32_t i = 0; i < m_iFields; ++i)
        {
            this->reverse(m_pArray + i);
        }


        uint64_t iOvercopied = iFields * m_iFieldSize - iLength;
        if (iOvercopied > 0)
        {
            this->shiftRight( iOvercopied );
        }

    }

    void reverseFields()
    {
        uint64_t iHalf = m_iFields/2;

        for (uint8_t i = 0; i < iHalf; ++i)
        {

            FIELDTYPE pTmp = m_pArray[i];
            m_pArray[i] = m_pArray[m_iFields-1-i];
            m_pArray[m_iFields-1-i] = pTmp;

        }

    }

    /**
     *
     * \brief copies the fields content from another instance to this instance
     *
     */
    void copy_content_fields(FIELDTYPE* pDest, FIELDTYPE* pSrc, uint64_t iMaxField, uint8_t iStartField)
    {

        for (uint64_t i = iStartField; i < iMaxField; ++i)
        {
            pDest[iStartField] = pSrc[iStartField];
        }


    }



    /*void copy_content_to_array_old(FIELDTYPE* pSrc, uint32_t iBitStart, uint32_t iBitCount)
    {

        uint8_t iOffset = iBitStart % m_iFieldSize;
        uint32_t iStartField = iBitStart / m_iFieldSize;
        //uint32_t iBitsRead = 0;
        uint32_t iMaxField = (iBitStart + iBitCount) / m_iFieldSize;


        if (iOffset == 0)
        {

            uint8_t iUnusedBits = this->m_iBits % m_iFieldSize;
            if (iUnusedBits != 0)
                iMaxField -= 1;

            for (size_t i = iStartField; i < iMaxField; ++i)
            {
                pSrc[i] = this->m_pArray[i];
            }

            if (iUnusedBits != 0)
            {
                iMaxField -= 1;

                FIELDTYPE oSrc = pSrc[iMaxField];
                oSrc = (oSrc >> (m_iFieldSize-iUnusedBits)) << (m_iFieldSize-iUnusedBits);

                FIELDTYPE  oOrigin = m_pArray[iMaxField];
                oOrigin = (oOrigin << iUnusedBits) >> iUnusedBits;

                pSrc[iMaxField] = oSrc | oOrigin;
            }

        } else {

            // first field: length-offset high bits
            FIELDTYPE newval = m_pArray[iStartField] << iOffset;

            // old field: offset bits
            FIELDTYPE oldval = pSrc[iStartField] << iOffset;
            oldval = oldval >> iOffset;

            pSrc[iStartField] = newval | oldval;

            size_t iBitsDone = m_iFieldSize - iOffset;
            for (size_t i = iStartField; i < iMaxField-1; ++i)
            {

                FIELDTYPE rem = m_pArray[i] >> (m_iFieldSize-iOffset);
                rem = rem | (m_pArray[i+1] << iOffset);

                pSrc[i] = rem;

                iBitsDone += 8;
            }

            uint8_t iBitsRemaining = m_iBits - iBitsDone;

            if (iBitsRemaining > 0)
            {

                FIELDTYPE rem = m_pArray[iMaxField-1] >> (m_iFieldSize-iOffset);
                rem = rem | (m_pArray[iMaxField] << iOffset);

                rem = rem << m_iFieldSize-iBitsRemaining;
                rem = rem >> m_iFieldSize-iBitsRemaining;

                pSrc[iMaxField] = pSrc[iMaxField] << iBitsRemaining;
                pSrc[iMaxField] = pSrc[iMaxField] >> iBitsRemaining;

                pSrc[iMaxField] = pSrc[iMaxField] | rem;
            }

        }

    }*/





    FIELDTYPE m_iUnusedBitsMask;


    MemoryPool<FIELDTYPE>* m_pPool=NULL;
    MemLoc m_oAlloc;
};


#endif //TSXCOUNT_BIGINT_H
