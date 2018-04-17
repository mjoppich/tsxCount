//
// Created by mjoppich on 11/10/17.
//

#ifndef TSXCOUNT_SEQUENCEUTILS_H
#define TSXCOUNT_SEQUENCEUTILS_H

#include <string>
#include <tsxutils/UBigInt.h>
#include <tsxcount/TSXTypes.h>

namespace TSXSeqUtils {


    class TSXSeqException: public std::exception
    {
    public:
        TSXSeqException(std::string sText)
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


    static std::string toSequence(TSX::tsx_kmer_t& oSeq)
    {

        if (oSeq.getBitCount() % 2 != 0)
        {
            throw new TSXSeqException("Sequence length not divisible by 2");
        }

        std::string sRetSeq;
        sRetSeq.resize(oSeq.getBitCount()/2);

        UBigInt oTest(oSeq);

        for (size_t i = 0; i < sRetSeq.size(); ++i)
        {

            UBigInt oVal = oTest;
            oVal.resize(2);

            oTest = oTest >> 2;

            uint8_t iVal = oVal.toUInt();

            switch (iVal)
            {
                case 0: sRetSeq[i] = 'A'; break;
                case 1: sRetSeq[i] = 'C'; break;
                case 2: sRetSeq[i] = 'G'; break;
                case 3: sRetSeq[i] = 'T'; break;

                default: sRetSeq[i] = 'N'; break;

            }

        }

        return sRetSeq;
    }

    static TSX::tsx_kmer_t fromSequence(std::string& seq)
    {

        UBigInt oRet(seq.length()*2, false);
        bool hadN = false;

        for (size_t i = 0; i < seq.length(); ++i)
        {

            //size_t iBitPos = 2*(seq.length() -1-i);
            size_t iBitPos = 2*i;

            switch (seq.at(i))
            {
                case 'A': // 0 00

                    oRet.setBit(iBitPos, 0);
                    oRet.setBit(iBitPos+1, 0);

                    break;

                case 'C': // 1 01

                    oRet.setBit(iBitPos, 1);
                    oRet.setBit(iBitPos+1, 0);

                    break;

                case 'G': // 2 10

                    oRet.setBit(iBitPos, 0);
                    oRet.setBit(iBitPos+1, 1);

                    break;

                case 'T': // 3 11
                    oRet.setBit(iBitPos, 1);
                    oRet.setBit(iBitPos+1, 1);
                    break;

                default: // random e.g. N

                    uint8_t iBit1 = rand() % 2;
                    uint8_t iBit2 = rand() % 2;

                    oRet.setBit(iBitPos, iBit1);
                    oRet.setBit(iBitPos+1, iBit2);

                    hadN = true;

                    break;
            }
        }

        /*
        UBigInt oTestStr = UBigInt::fromString("11010101100000100000");
        if (oRet == oTestStr)
        {
            std::cout << oTestStr.to_string() << " 11010101100000100000" << std::endl;
            std::cout << seq << " " << TSXSeqUtils::toSequence(oRet) << std::endl;

            std::cout << std::endl;
        }
         */

        if (hadN)
        {
            std::string sTransSeq = TSXSeqUtils::toSequence(oRet);

            std::cout << "HADN: " << seq << "\t" << sTransSeq << std::endl;
        }

        return oRet;

    }


}



#endif //TSXCOUNT_SEQUENCEUTILS_H
