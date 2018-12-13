
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMapPThread.h>
#include "testExecution.h"


int main(int argc, char *argv[])
{

    uint32_t itK = 8;
    TSXHashMapPThread* ptMap = new TSXHashMapPThread(8, 4, itK, 2);
    
    testHashMap(ptMap, true);
    
    return 1;


    CLParser oParser(argc, argv);

    oParser.setArgument("k", std::to_string(10));
    oParser.setArgument("fastq", "/mnt/c/ownCloud/data/tsx/small.fq");

    FASTXreader<FASTQEntry>* pReader = FASTXreader<FASTQEntry>::createFQReader(&oParser);
    size_t iK = oParser.isSet("k") ? (size_t) oParser.getIntArgument("k") : 15;

    TSXHashMapPThread* pMap = new TSXHashMapPThread(16, 6, iK);

    pMap->testHashFunction();

    bool exitNow = false;

    while (pReader->hasNext() and !exitNow)
    {

        std::vector<FASTQEntry> oEntries = pReader->readEntries(100);

        for (size_t i = 0; i < oEntries.size(); ++i)
        {
            FASTQEntry* pEntry = &(oEntries.at(i));

            std::string sSeq = pEntry->getSequence();
            std::vector<TSX::tsx_kmer_t> allKmers = createKMers(sSeq, iK, pMap->getMemoryPool());

            for (auto kmer : allKmers)
            {
                pMap->countKmer(kmer);

                exitNow=false;
            }

        }


    }

    std::cout << "Added a total of " << pMap->getKmerCount() << " different kmers" << std::endl;

    std::vector<TSX::tsx_kmer_t> allKmers = pMap->getAllKmers();

    std::ofstream("/mnt/c/ownCloud/data/tsx/small.fq");

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
