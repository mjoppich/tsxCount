
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>


void evaluate(TSXHashMap* pMap, UBigInt& kmer, const size_t iRefCount)
{

    UBigInt oRes1 = pMap->getKmerCount(kmer);
    uint32_t iKmer1Count = oRes1.toUInt();

    if (iKmer1Count != iRefCount)
    {
        oRes1 = pMap->getKmerCount(kmer);
        iKmer1Count = oRes1.toUInt();
    }

    std::cerr << "kmer: " << kmer.to_string() << ": " << oRes1.to_string() << " " << std::to_string(iKmer1Count) << std::endl;

}

void testHashMap()
{
    UBigInt m_iKeyMask(128, true);
    uint64_t iBla = 0x8000000800000ULL;
    m_iKeyMask = UBigInt::createFromBitShift(128, 64);

    std::cout << "SHIFT1 " << m_iKeyMask.to_string() << std::endl;


    m_iKeyMask = m_iKeyMask + 1;

    std::cout << "BIGINT " << m_iKeyMask.to_string() << std::endl;
    std::cout << "BIGINT " << m_iKeyMask.to_string() << std::endl;

    m_iKeyMask = m_iKeyMask >> 4;

    std::cout << "SHIFT4 " << m_iKeyMask.to_string() << std::endl;
    std::cout << "SHIFT4 " << m_iKeyMask.to_string() << std::endl;

    m_iKeyMask = m_iKeyMask << 4;

    std::cout << "BSHFT4 " << m_iKeyMask.to_string() << std::endl;
    std::cout << "BSHFT4 " << m_iKeyMask.to_string() << std::endl;

    m_iKeyMask = m_iKeyMask >> 20;

    std::cout << "SHIF20 " << m_iKeyMask.to_string() << std::endl;
    std::cout << "SHIF20 " << m_iKeyMask.to_string() << std::endl;

    m_iKeyMask = m_iKeyMask << 20;

    std::cout << "BSHF20 " << m_iKeyMask.to_string() << std::endl;
    std::cout << "BSHF20 " << m_iKeyMask.to_string() << std::endl;

    uint32_t iK = 8;

    TSXHashMap* pMap = new TSXHashMap(8, 4, iK);


    UBigInt oTest ("110110010001");
    std::cout << (oTest >> 7).to_string() << std::endl;


    UBigInt oKmer1(0);
    UBigInt oKmer2(1);
    UBigInt oKmer3(178);

    oKmer1.resize( 2*iK );
    oKmer2.resize( 2*iK );
    oKmer3.resize( 2*iK );



    /*
    oTest = ~oTest;
    std::cerr << oTest.to_string() << std::endl;
    std::cerr << (oTest == 0) << std::endl;
*/
/*
    pMap->addKmer( oKmer1 );
    pMap->addKmer( oKmer2 );
    pMap->addKmer( oKmer1 );
*/

    const size_t iMaxCount = 2048*4;

    size_t i;
    for ( i = 0; i < iMaxCount; ++i)
    {
        std::cerr << "adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
        pMap->addKmer( oKmer1 );


        if ((i % 2) == 0)
        {
            std::cerr << "adding kmer: " << oKmer2.to_string() << " " << std::to_string(i) << std::endl;
            pMap->addKmer( oKmer2 );
        }


        if ((i % 4) == 0)
        {
            std::cerr << "adding kmer: " << oKmer3.to_string() << " " << std::to_string(i) << std::endl;
            pMap->addKmer(oKmer3);
        }
    }


    evaluate(pMap, oKmer1, iMaxCount);
    evaluate(pMap, oKmer2, iMaxCount);
    evaluate(pMap, oKmer3, iMaxCount);
}

std::vector<TSX::tsx_kmer_t> createKMers(std::string& sSequence, size_t iK)
{
    std::vector<TSX::tsx_kmer_t> oRetVec;

    for (size_t i = 0; i < sSequence.length()-iK+1; ++i)
    {
        std::string sSubSeq = sSequence.substr(i, iK);
        TSX::tsx_kmer_t oKmer = TSXSeqUtils::fromSequence(sSubSeq);
        oRetVec.push_back(oKmer);
    }

    return oRetVec;

}

int main(int argc, char *argv[])
{

    CLParser oParser(argc, argv);

    oParser.setArgument("k", std::to_string(10));
    oParser.setArgument("fastq", "/mnt/c/ownCloud/data/tsx/small.fq");

    FASTXreader<FASTQEntry>* pReader = FASTXreader<FASTQEntry>::createFQReader(&oParser);
    size_t iK = oParser.isSet("k") ? (size_t) oParser.getIntArgument("k") : 15;

    TSXHashMap* pMap = new TSXHashMap(16, 6, iK);

    pMap->testHashFunction();

    bool exitNow = false;

    while (pReader->hasNext() and !exitNow)
    {

        std::vector<FASTQEntry> oEntries = pReader->readEntries(100);

        for (size_t i = 0; i < oEntries.size(); ++i)
        {
            FASTQEntry* pEntry = &(oEntries.at(i));

            std::string sSeq = pEntry->getSequence();
            std::vector<TSX::tsx_kmer_t> allKmers = createKMers(sSeq, iK);

            for (auto kmer : allKmers)
            {
                pMap->addKmer(kmer);

                exitNow=false;
            }

        }


    }

    std::cout << "Added a total of " << pMap->getKmerCount() << " different kmers" << std::endl;

    std::vector<TSX::tsx_kmer_t> allKmers = pMap->getAllKmers();

    std::ofstream("/mnt/c/ownCloud/data/tsx/small.fq")

#pragma omp parallel for
    for (size_t i=0; i < allKmers.size(); ++i)
    {

        TSX::tsx_kmer_t& kmer = allKmers.at(i);

        UBigInt oCount = pMap->getKmerCount(kmer);

#pragma omp critical
        {
            //std::cout << kmer.to_string() << "\t" << TSXSeqUtils::toSequence(kmer) << "\t" << oCount.toUInt() << std::endl;
        }
    }

    std::cout << "Printed " << allKmers.size() << " kmers" << std::endl;



    return 0;
}
