//
// Created by mjopp on 24/10/2018.
//

#ifndef TSXCOUNT_TESTEXECUTION_H
#define TSXCOUNT_TESTEXECUTION_H

#include <vector>
#include <tsxcount/TSXTypes.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMap.h>

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

void testHashMap(TSXHashMap* pMap, bool parallel=false)
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




    UBigInt oTest ("110110010001");
    std::cout << (oTest >> 7).to_string() << std::endl;


    UBigInt oKmer1(0);
    UBigInt oKmer2(1);
    UBigInt oKmer3(178);

    oKmer1.resize( 2*pMap->getK() );
    oKmer2.resize( 2*pMap->getK() );
    oKmer3.resize( 2*pMap->getK() );


    const size_t iMaxCount = 2048*4;

    if (!parallel)
    {
        omp_set_num_threads(1);
    } else {
        omp_set_num_threads(pMap->getThreads());
    }

    std::cerr << "Parellel=" << parallel << " Running on " << omp_get_num_threads() << " threads on map configured for " << pMap->getThreads() << std::endl;

    size_t i;

#pragma omp parallel for
    for ( i = 0; i < iMaxCount; ++i)
    {

        std::cerr << "Thread "<< omp_get_thread_num() << " adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
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

#endif //TSXCOUNT_TESTEXECUTION_H
