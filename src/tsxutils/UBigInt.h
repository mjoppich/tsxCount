//
// Created by joppich on 6/28/16.
//

#ifndef TSXCOUNT_BIGINT_H
#define TSXCOUNT_BIGINT_H


#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <exception>
#include <sstream>
#include <functional>
#include <iostream>
#include <bitset>


#ifndef BITSTOFIELDS
#define BITSTOFIELDS(x) (x >> 3)
#endif

typedef uint8_t FIELDTYPE;

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

    static UBigInt createFromBitShift(uint32_t iBits, uint32_t iPosition)
    {

        UBigInt oRet = UBigInt(iBits, true);

        oRet = 1;

        oRet = oRet << iPosition;

        return oRet;

    }


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

    /**
     *
     * @param pPosition start position of array
     * @param iIdx Index in Array
     * @param iOffsetBits starting bit in Array
     * @param iLength of memory block in bits
     * @return
     */
    static UBigInt createFromMemory(uint8_t* pPosition, uint64_t iIdx, uint8_t iOffsetBits, uint32_t iLength)
    {

        uint8_t* pPos = pPosition + iIdx;

        UBigInt oReturn = UBigInt(iLength, true);

        oReturn.copy_content_bits(pPos, iOffsetBits, iLength);

    }

    UBigInt(uint64_t iBits, uint8_t* pMemPos, uint8_t iOffset)
    : UBigInt(iBits, true)
    {
        uint32_t iFieldsNeeded = std::ceil(iBits / this->m_iFieldSize);

        UBigInt oIntermediate = UBigInt(iFieldsNeeded * this->m_iFieldSize, false);
        oIntermediate.copy_content_from_le((FIELDTYPE*) pMemPos, iBits);

    }

    /**
     * \brief default constructor
     */
    UBigInt()
    : UBigInt(64, true)
    {

    }

    /** \brief UBigInt constructor
      * \param iBits number of bits required.
      * \param bZero if set true, the object will have value 0
      *
      * This implements a new UBigInt constructor. The Array can be directly initialized with 0
      */
    UBigInt(uint64_t iBits, bool bZero)
    {

        initialize(iBits);

    }

    UBigInt(UBigInt & oOther )
     : UBigInt(oOther.m_iBits, false)
    {
        copy_content(oOther);
    }

    UBigInt(const UBigInt & oOther )
            : UBigInt(oOther.m_iBits, false)
    {
        copy_content(oOther);
    }

    ~UBigInt()
    {
        this->reset();
    }

    UBigInt(uint64_t iValue)
    : UBigInt(64, true)
    {

        this->copy_content_bits((FIELDTYPE*)&iValue, 0, 64);

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
        this->initialize(oOther.m_iBits);
        this->copy_content(oOther);

        return *this;

    }

    UBigInt operator- (const uint64_t iLeft)
    {
        UBigInt oLeft = iLeft;

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
        UBigInt oLeft = iLeft;

        return operator+(oLeft);
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
        UBigInt oRet(this->self() );
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

    bool operator == (const UBigInt & oOther)
    {


        throw "not yet implemented";




    }

    uint64_t toUInt()
    {
        if (m_iBits <= m_iFieldSize)
        {

            uint32_t iDiff = m_iFieldSize-m_iBits;

            uint64_t iField = m_pArray[0];

            iField = (iField << iDiff) >> iDiff;

            return iField;

        } else {

            // iBits > iFieldSize
            // assemble this uint64_t from multiple fields

            uint32_t iNeeded = std::ceil( 64.0 / m_iFieldSize );
            uint64_t iReturn = 0;

            for (uint8_t i = iNeeded; i > 0; --i)
            {

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

        UBigInt oRet(iBitPos-1, true);
        oRet.copy_content_bits( this->m_pArray, 0, iBitPos-1);
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
    void initialize(uint32_t iBits)
    {

        m_iBits = iBits;
        m_iFields = std::ceil( iBits / (8.0f*sizeof(FIELDTYPE)) );

        m_pArray = new FIELDTYPE[m_iFields];

        // NULLING the array
        for (uint32_t i = 0; i < m_iFields; ++i)
            m_pArray[i] = 0;


        m_iUnusedBitsMask = this->createUnusedBitsMask();


    }

    FIELDTYPE createUnusedBitsMask()
    {
        FIELDTYPE iUnusedBitsMask = -1;

        uint64_t iRemainingBits = m_iFields * sizeof(FIELDTYPE) - m_iBits;

        if (iRemainingBits > 0)
        {
            m_iUnusedBitsMask -= ( (1 << (m_iFieldSize-iRemainingBits)) -1);
        }

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


    void shiftLeft( uint32_t iShift)
    {

        if (iShift > m_iFieldSize)
        {

            uint64_t iFields = iShift / m_iFieldSize;

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

        for (uint32_t i = 0; i < m_iFields-1; ++i)
        {
            FIELDTYPE iCurField = (m_pArray[m_iFields-1-i] << iShift);
            FIELDTYPE iNextFieldPart = (m_pArray[m_iFields-1-i-1] >> iCounterShift);

            m_pArray[m_iFields-1 -i] = iCurField | iNextFieldPart;

        }

        m_pArray[0] = (m_pArray[0] << iShift);


    }

    void shiftRight( uint32_t iShift)
    {

        if (iShift > m_iFieldSize)
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

    void bitAnd(const UBigInt & oOther)
    {
        const FIELDTYPE* pOtherField = oOther.m_pArray;
        this->applyToAllFields( [pOtherField] (FIELDTYPE iValue, uint64_t iField) {return (iValue & pOtherField[iField]);} );
    }

    void bitOr(const UBigInt & oOther)
    {

        const FIELDTYPE* pOtherField = oOther.m_pArray;
        this->applyToAllFields( [pOtherField] (FIELDTYPE iValue, uint64_t iField) {return (iValue | pOtherField[iField]);} );

    }

    void bitXor(const UBigInt & oOther)
    {

        const FIELDTYPE* pOtherField = oOther.m_pArray;
        this->applyToAllFields( [pOtherField] (FIELDTYPE iValue, uint64_t iField) {return (iValue ^ pOtherField[iField]);} );

    }

    void applyToAllFields(std::function<FIELDTYPE(FIELDTYPE, uint64_t)> oOp)
    {
        for (uint64_t i = 0; i < m_iFields; ++i)
        {
            m_pArray[i] = oOp(m_pArray[i], i);
        }
    }

    void sub(const UBigInt & oOther)
    {

        this->bitCompl();
        this->add(oOther);
        this->bitCompl();

    }

    /**
     *
     * @brief resizes this UBigInt such that it can hold iBits
     *
     * @param iBits number of bits this element can hold
     */
    void resize(uint32_t iBits)
    {
        FIELDTYPE* pOldArray = m_pArray;
        uint32_t iOldFields = m_iFields;

        m_pArray = NULL;

        this->initialize(iBits);

        m_iUnusedBitsMask = this->createUnusedBitsMask();

        if (pOldArray != NULL)
        {
            // copy old data over
            for (uint32_t i = 0; i < iOldFields; ++i)
            {
                m_pArray[i] = pOldArray[i];
            }

            delete pOldArray;

        }

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
        iLastField = (~(iLastField & m_iUnusedBitsMask)) & m_iUnusedBitsMask;

        m_pArray[m_iFields-1] = iLastField;

    }

    void reset()
    {
        if (m_pArray != NULL)
        {
            delete m_pArray;
            m_pArray == NULL;
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

        if (this->m_iBits < iLength)
            throw new UBigIntNotLargeEnough(this, iLength);

        const uint8_t iFields = std::ceil( iLength / m_iFieldSize );

        uint8_t i = 0;
        for (i = 0; i < iFields; ++i)
        {
            FIELDTYPE iInsert = pSrc[i];
            m_pArray[i] = iInsert;
        }

        uint64_t iOvercopied = iFields * m_iFieldSize - iLength;
        if (iOvercopied > 0)
        {

            m_pArray[i-1] = (m_pArray[i-1] << iOvercopied) >> iOvercopied;

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


    /**
     *
     * @brief copies x bits from another instance to this instance
     * @param pSrc the source field
     * @param iBitStart bit start
     * @param iBitCount bitcount
     *
     */
    void copy_content_bits(FIELDTYPE* pSrc, uint64_t iBitStart, uint64_t iBitCount)
    {

        uint8_t iOffset = iBitStart % m_iFieldSize;
        uint64_t iStartField = iBitStart / m_iFieldSize;
        uint64_t iBitsRead = 0;
        uint64_t iMaxField = (iBitStart + iBitCount) / m_iFieldSize;

        if (iOffset == 0)
        {

            this->copy_content_from_be(pSrc, iBitCount);

        } else {


            for (uint64_t i = iStartField; i < iMaxField; ++i)
            {
                FIELDTYPE iValue = 0;

                iValue |= (pSrc[i+1] << (m_iFieldSize - iOffset) | (pSrc[ i ] >> iOffset));

                m_pArray[i] = iValue;

                iBitsRead += m_iFieldSize;

                if ( iBitCount - iBitsRead < m_iFieldSize)
                {

                    break;

                }


            }


            uint8_t iBitsRemaining = iBitCount - iBitsRead;

            if (iBitsRemaining < iOffset)
            {
                FIELDTYPE iValue = 0;

                iValue |= (pSrc[iMaxField] >> iBitsRemaining);

                m_pArray[iMaxField] = iValue;
            } else {
                FIELDTYPE iValue = 0;

                FIELDTYPE iMask = ~0 >> iBitsRemaining;
                iValue |=  ((pSrc[iMaxField+1] << (m_iFieldSize - iOffset) | (pSrc[ iMaxField ] >> iOffset)) & iMask);
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

        this->copy_content_bits(oOther.m_pArray, 0, this->m_iBits);

    }

    uint32_t m_iBits = 0;
    uint32_t m_iFields = 0;

    FIELDTYPE* m_pArray = NULL;
    static const uint8_t m_iFieldSize = sizeof(FIELDTYPE) * 8;

    FIELDTYPE m_iUnusedBitsMask;

};


#endif //TSXCOUNT_BIGINT_H
