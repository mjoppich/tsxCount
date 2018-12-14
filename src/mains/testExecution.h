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

    std::cerr << "kmer: " << kmer.to_string() << ": " << oRes1.to_string() << " " << std::to_string(iKmer1Count) << " Should be " << iRefCount <<  std::endl;

}

void testHashMap(TSXHashMap* pMap, bool parallel=false)
{

    /*

    UBigInt m_iKeyMask(128, true, pMap->getMemoryPool());
    uint64_t iBla = 0x8000000800000ULL;
    m_iKeyMask = UBigInt::createFromBitShift(128, 64, pMap->getMemoryPool());

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



    UBigInt oTest ("110110010001", pMap->getMemoryPool());
    std::cout << (oTest >> 7).to_string() << std::endl;
    */

    UBigInt oKmer1(0, pMap->getMemoryPool());
    UBigInt oKmer2(788, pMap->getMemoryPool());
    UBigInt oKmer3(1, pMap->getMemoryPool());
    UBigInt oKmer4(178, pMap->getMemoryPool());

    oKmer1.resize( 2*pMap->getK() );
    oKmer2.resize( 2*pMap->getK() );
    oKmer3.resize( 2*pMap->getK() );
    oKmer4.resize( 2*pMap->getK() );


    std::cout << "Kmer1 " << oKmer1.to_string() << std::endl;
    std::cout << "Kmer2 " << oKmer2.to_string() << std::endl;
    std::cout << "Kmer3 " << oKmer3.to_string() << std::endl;
    std::cout << "Kmer4 " << oKmer4.to_string() << std::endl;


    const size_t iMaxCount = 2048*4*16;

    uint8_t threads = 4;
    pMap->setThreads(threads);

    std::cout << "Running on " << (int) threads << " threads" << std::endl;
omp_set_dynamic(0);

    if (!parallel)
    {
        omp_set_num_threads(1);
        threads = 1;
    } else {
        threads = pMap->getThreads();
        omp_set_num_threads(pMap->getThreads());
    }


    size_t i;

	std::cout << "Threads specified: " << (int) threads << "/" << pMap->getThreads()<<"/"<< omp_get_max_threads() << " OMP GET NUM THREADS " << omp_get_num_threads() << std::endl;

#pragma omp parallel num_threads(threads)
    {
        int iThreadCount = omp_get_num_threads();
#pragma omp critical
        {
            std::cerr << "Parellel=" << parallel << " Running on thread " << omp_get_thread_num() << " on map configured for " << threads << " (omp threads) " << iThreadCount << std::endl;
		
        }
}



        #pragma omp parallel for private(i) num_threads(threads)
        for (i = 0; i < iMaxCount; ++i) {


            //std::cerr << "Thread " << omp_get_thread_num() << " adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
            pMap->addKmer(oKmer1);

            if ((!parallel) && (i % 1000 == 0))
            {
                UBigInt oRes1 = pMap->getKmerCount(oKmer1);
                uint32_t iKmer1Count = oRes1.toUInt();
                std::cout << "Kmer1 " << iKmer1Count << " " << i << " " << oKmer1.to_string() << " " << omp_get_thread_num() << std::endl;
            }

            sleep(0.1);

            //
            //
            //

            if ((i % 2) == 0) {
                //std::cerr << "adding kmer: " << oKmer2.to_string() << " " << std::to_string(i) << std::endl;
                pMap->addKmer(oKmer2);
            } else {
                pMap->addKmer(oKmer3);

            }


            if ((i % 4) == 0) {
                //std::cerr << "adding kmer: " << oKmer3.to_string() << " " << std::to_string(i) << std::endl;
                pMap->addKmer(oKmer4);
            }
        }




    evaluate(pMap, oKmer1, iMaxCount);
    evaluate(pMap, oKmer2, iMaxCount/2);
    evaluate(pMap, oKmer3, iMaxCount/2);
    evaluate(pMap, oKmer4, iMaxCount/4);
}

void testHashMapOld(TSXHashMap* pMap, bool parallel=false)
{

    /*

    UBigInt m_iKeyMask(128, true, pMap->getMemoryPool());
    uint64_t iBla = 0x8000000800000ULL;
    m_iKeyMask = UBigInt::createFromBitShift(128, 64, pMap->getMemoryPool());

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



    UBigInt oTest ("110110010001", pMap->getMemoryPool());
    std::cout << (oTest >> 7).to_string() << std::endl;
    */

    UBigInt oKmer1(0, pMap->getMemoryPool());
    UBigInt oKmer2(788, pMap->getMemoryPool());
    UBigInt oKmer3(1, pMap->getMemoryPool());
    UBigInt oKmer4(178, pMap->getMemoryPool());

    oKmer1.resize( 2*pMap->getK() );
    oKmer2.resize( 2*pMap->getK() );
    oKmer3.resize( 2*pMap->getK() );
    oKmer4.resize( 2*pMap->getK() );


    std::cout << "Kmer1 " << oKmer1.to_string() << std::endl;
    std::cout << "Kmer2 " << oKmer2.to_string() << std::endl;
    std::cout << "Kmer3 " << oKmer3.to_string() << std::endl;
    std::cout << "Kmer4 " << oKmer4.to_string() << std::endl;


    const size_t iMaxCount = 2048*4*16;

    uint8_t threads = 4;
    pMap->setThreads(threads);

    std::cout << "Running on " << (int) threads << " threads" << std::endl;
    omp_set_dynamic(0);

    if (!parallel)
    {
        omp_set_num_threads(1);
        threads = 1;
    } else {
        threads = pMap->getThreads();
        omp_set_num_threads(pMap->getThreads());
    }


    size_t i;

    std::cout << "Threads specified: " << (int) threads << "/" << pMap->getThreads()<<"/"<< omp_get_max_threads() << " OMP GET NUM THREADS " << omp_get_num_threads() << std::endl;

#pragma omp parallel num_threads(threads)
    {
        int iThreadCount = omp_get_num_threads();
#pragma omp critical
        {
            std::cerr << "Parellel=" << parallel << " Running on thread " << omp_get_thread_num() << " on map configured for " << threads << " (omp threads) " << iThreadCount << std::endl;

        }
        sleep(3);
    }



#pragma omp parallel for private(i) num_threads(threads)
    for (i = 0; i < iMaxCount; ++i) {


        //std::cerr << "Thread " << omp_get_thread_num() << " adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
        pMap->addKmer(oKmer1);

        if ((!parallel) && (i % 1000 == 0))
        {
            UBigInt oRes1 = pMap->getKmerCount(oKmer1);
            uint32_t iKmer1Count = oRes1.toUInt();
            std::cout << "Kmer1 " << iKmer1Count << " " << i << " " << oKmer1.to_string() << " " << omp_get_thread_num() << std::endl;
        }

        //
        //
        //

        if ((i % 2) == 0) {
            //std::cerr << "adding kmer: " << oKmer2.to_string() << " " << std::to_string(i) << std::endl;
            pMap->addKmer(oKmer2);
        } else {
            pMap->addKmer(oKmer3);

        }


        if ((i % 4) == 0) {
            //std::cerr << "adding kmer: " << oKmer3.to_string() << " " << std::to_string(i) << std::endl;
            pMap->addKmer(oKmer4);
        }
    }




    evaluate(pMap, oKmer1, iMaxCount);
    evaluate(pMap, oKmer2, iMaxCount/2);
    evaluate(pMap, oKmer3, iMaxCount/2);
    evaluate(pMap, oKmer4, iMaxCount/4);
}

std::vector<TSX::tsx_kmer_t> createKMers(std::string& sSequence, size_t iK, MemoryPool<FIELDTYPE>* pPool)
{
    std::vector<TSX::tsx_kmer_t> oRetVec;

    for (size_t i = 0; i < sSequence.length()-iK+1; ++i)
    {
        std::string sSubSeq = sSequence.substr(i, iK);
        TSX::tsx_kmer_t oKmer = TSXSeqUtils::fromSequence(sSubSeq, pPool);
        oRetVec.push_back(oKmer);
    }

    return oRetVec;

}

#endif //TSXCOUNT_TESTEXECUTION_H
