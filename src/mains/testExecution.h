//
// Created by mjopp on 24/10/2018.
//

#ifndef TSXCOUNT_TESTEXECUTION_H
#define TSXCOUNT_TESTEXECUTION_H

#include <vector>
#include <tsxcount/TSXTypes.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMap.h>
#include <fastxutils/FastXReader.h>
#include <fstream>

std::vector<std::string> createKMers(std::string& sSequence, size_t iK, MemoryPool<FIELDTYPE>* pPool)
{
    std::vector<std::string> oRetVec;

    if (sSequence.size() < iK)
        return oRetVec;

    for (size_t i = 0; i < sSequence.length()-iK+1; ++i)
    {
        std::string sSubSeq = sSequence.substr(i, iK);

        if (sSubSeq.size() != iK)
        {
            continue;
        }

        oRetVec.push_back(sSubSeq);
    }

    return oRetVec;

}

void evaluate(TSXHashMap* pMap, UBigInt& kmer, const size_t iRefCount, const std::string* pKmer=NULL)
{

    UBigInt oRes1 = pMap->getKmerCount(kmer);
    uint32_t iKmer1Count = oRes1.toUInt();

    //if (iKmer1Count != iRefCount)
    //{
    //    oRes1 = pMap->getKmerCount(kmer);
    //    iKmer1Count = oRes1.toUInt();
    //}

    if (iRefCount != iKmer1Count)
    {
#pragma omp critical
        {
            std::cout << "kmer: " << kmer.to_string() << "( " << TSXSeqUtils::toSequence(kmer) << " )" << ": " << oRes1.to_string() << " " << std::to_string(iKmer1Count) << " Should be " << iRefCount << std::endl;


            UBigInt oRes2 = pMap->getKmerCount(kmer, true);
            if (pKmer != NULL)
            {
                std::cout << " " << *pKmer;
            }

            std::cout <<  std::endl;
        };

    }



}

std::vector<std::string> split(const std::string& str, char delim = '\t')
{

    std::vector<std::string> cont;
    std::size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos) {
        cont.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    cont.push_back(str.substr(previous, current - previous));

    return cont;
}

std::map<std::string, uint32_t> loadReferences(std::string sFilename)
{

    std::cerr << "Loading reference file: " << sFilename << std::endl;
    std::ifstream file(sFilename);

    std::map<std::string, uint32_t> kmer2count;

    if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {

            std::vector<std::string> vElems = split(line);

            uint32_t iCount = std::atoi(vElems[1].c_str());
            kmer2count[vElems[0]] = iCount;


        }
        file.close();
    }

    return kmer2count;

}

void testHashMap(TSXHashMap* pMap, bool parallel=false)
{

    const size_t iMaxCount = 2048*4*16;

    uint8_t threads = pMap->getThreads();
    omp_set_num_threads(pMap->getThreads());

    std::cout << "Running on " << (int) threads << " threads" << std::endl;
    omp_set_dynamic(0);




    std::string sFileName = "/mnt/d/owncloud/data/tsx/usmall_t7.fastq";
    //sFileName = "/mnt/d/owncloud/data/tsx/small_t7.1000.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small_t7.3000.fastq";
    //sFileName = "/mnt/d/owncloud/data/tsx/small_t7.5000.fastq";
    //sFileName = "/mnt/d/owncloud/data/tsx/small_t7.8000.fastq";
    //sFileName = "/mnt/d/owncloud/data/tsx/small_t7.fastq";
    //sFileName = "/mnt/d/owncloud/data/tsx/small2_t7.fastq";

    FASTXreader<FASTQEntry>* pReader = new FASTXreader<FASTQEntry>(&sFileName);

    uint32_t iK = pMap->getK();


#pragma omp parallel num_threads(threads)
    {
#pragma omp master
        {
            bool exitNow = false;
            bool verbose = true;

            while (pReader->hasNext())
            {

                std::vector<FASTQEntry>* pEntries = pReader->getEntries(10);
                uint64_t iPosOfInterest = 0;

#pragma omp task firstprivate(pEntries) shared(iPosOfInterest)
                {

#pragma omp critical
                    {
                        std::cout << "In thread " << omp_get_thread_num() << " reading entries " << pEntries->size() << std::endl;
                    }

                    MemoryPool<FIELDTYPE>* pPool = pMap->getMemoryPool();

                    for (size_t i = 0; i < pEntries->size(); ++i)
                    {
                        FASTQEntry* pEntry = &(pEntries->at(i));

                        std::string sSeq = pEntry->getSequence();
                        std::vector<std::string> allKmers = createKMers(sSeq, iK, pPool);

                        uint32_t iAddedKmers = 0;
                        uint32_t iTotalKmers = allKmers.size();

                        std::string targetKmer = "TCCATCCATTCCAT";

                        for (auto kmerStr : allKmers)
                        {

                            TSX::tsx_kmer_t oKmer = TSXSeqUtils::fromSequence(kmerStr, pPool);
                            bool vadd = false;
                            std::vector<uint64_t> vAllPos;

                            //if (iAddedKmers % 10 == 0)
                            //{
                            //   std::cout << "processed " << iAddedKmers << " of " << iTotalKmers << std::endl;
                            //}


                            if (verbose)
                            {


                                if (kmerStr == targetKmer)
                                {
                                    UBigInt oRes1 = pMap->getKmerCount(oKmer);

                                    uint32_t iKmer1Count = oRes1.toUInt();
                                    std::cerr << kmerStr << " " << iKmer1Count << std::endl;
                                }

                                if (iPosOfInterest != 0)
                                {
                                    vAllPos = pMap->getKmerPositions(oKmer);
                                }

                                if (iPosOfInterest != 0)
                                {

                                    for (auto iPos : vAllPos)
                                    {
                                        if ((iPosOfInterest-3 <= iPos) && (iPos <= iPosOfInterest +3))
                                        {
                                            std::cerr << "POI usage BEFORE " << kmerStr << " " << iPos << " " << iPosOfInterest << std::endl;

                                            std::string stest = targetKmer;

                                            TSX::tsx_kmer_t oTestKmer = TSXSeqUtils::fromSequence(stest, pPool);
                                            UBigInt oRes1 = pMap->getKmerCount(oTestKmer);
                                            uint32_t iKmer1Count = oRes1.toUInt(true);
                                            std::cerr << stest << " " << iKmer1Count << " OB" << std::endl;
                                            std::cerr << std::endl;
                                            vadd=true;
                                        }
                                    }
                                }

                            }

                            pMap->addKmer(oKmer, verbose); //vadd
                            //iAddedKmers += 1;


                            if (verbose)
                            {
                                if (iPosOfInterest != 0)
                                {

                                    for (auto iPos : vAllPos)
                                    {
                                        if ((iPosOfInterest-3 <= iPos) && (iPos <= iPosOfInterest +3))
                                        {
                                            std::cerr << "POI usage AFTER " << kmerStr << " " << iPos << " " << iPosOfInterest << std::endl;

                                            std::string stest = targetKmer;

                                            TSX::tsx_kmer_t oTestKmer = TSXSeqUtils::fromSequence(stest, pPool);
                                            UBigInt oRes1 = pMap->getKmerCount(oTestKmer);
                                            uint32_t iKmer1Count = oRes1.toUInt(true);
                                            std::cerr << stest << " " << iKmer1Count << " OA" << std::endl;
                                        }
                                    }
                                }

                                if (kmerStr == targetKmer)
                                {
                                    std::vector<uint64_t> vPos = pMap->getKmerPositions(oKmer);
                                    for (auto iPos : vPos)
                                    {
                                        std::cerr << kmerStr << " vPos " << iPos << std::endl;
                                        iPosOfInterest = iPos;
                                    }

                                    UBigInt oRes1 = pMap->getKmerCount(oKmer, vadd);
                                    uint32_t iKmer1Count = oRes1.toUInt();
                                    std::cerr << kmerStr << " " << iKmer1Count << "E" << std::endl;
                                }
                            }


                        }


                    }

#pragma omp critical
                    {
                      pMap->print_stats();
                    }

                    delete pEntries;
                }



            }


        }

    }



    std::cout << "Added a total of " << pMap->getKmerCount() << " different kmers" << std::endl;

    std::map<std::string, uint32_t> kmer2c = loadReferences(sFileName + "." + std::to_string(iK) + ".count");


    // checking counts
    std::cout << "Checking kmer counts against manual hashmap ..." << std::endl;
    std::cout << "Manual counts: " << kmer2c.size() << std::endl;



#pragma omp parallel
    {
#pragma omp single
        {
            for (auto& elem: kmer2c)
            {
#pragma omp task
                {
                    // Do something with x, e.g.
                    UBigInt tkmer = TSXSeqUtils::fromSequenceD(elem.first, pMap->getMemoryPool());
                    evaluate(pMap, tkmer, elem.second, &(elem.first));
                }
            }
        }
    }


    std::cout << "Kmer count check completed." << std::endl;


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

    UBigInt oKmer1(10067, pMap->getMemoryPool());
    UBigInt oKmer2(2786, pMap->getMemoryPool());
    UBigInt oKmer3(9816, pMap->getMemoryPool());
    UBigInt oKmer4(156, pMap->getMemoryPool());

    oKmer1.resize( 2*pMap->getK() );
    oKmer2.resize( 2*pMap->getK() );
    oKmer3.resize( 2*pMap->getK() );
    oKmer4.resize( 2*pMap->getK() );


    std::cout << "Kmer1 " << oKmer1.to_string() << std::endl;
    std::cout << "Kmer2 " << oKmer2.to_string() << std::endl;
    std::cout << "Kmer3 " << oKmer3.to_string() << std::endl;
    std::cout << "Kmer4 " << oKmer4.to_string() << std::endl;


    const size_t iMaxCount = 2048*4*24;

    uint8_t threads = pMap->getThreads();

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



#pragma omp parallel for schedule(dynamic) private(i) num_threads(threads)
    for (i = 0; i < iMaxCount; ++i) {


        //std::cerr << "Thread " << omp_get_thread_num() << " adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
        pMap->addKmer(oKmer1);

        if ((threads <=4) && (i % 1000 == 0))
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



#endif //TSXCOUNT_TESTEXECUTION_H
