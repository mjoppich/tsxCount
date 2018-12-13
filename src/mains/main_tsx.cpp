
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMapTSX.h>
#include "testExecution.h"
#include <stddef.h>
#include <immintrin.h>




int main(int argc, char *argv[])
{

    /*
    int x = 0;
    int itries = 0;

    while(true) {

        uint status = _xbegin();
        if(status != _XBEGIN_STARTED) {
            //printf("Transaction failed, retrying\n");
            itries += 1;
            continue;
        }

        x++;
        //printf("Transaction done\n");

        _xend();

        std::cout << "status " << status << " retries " << itries << std::endl;
        break;
    }


    exit(0);
     */

    uint32_t itK = 8;
    TSXHashMapTSX* ptMap = new TSXHashMapTSX(8, 4, itK, 2);

    testHashMap(ptMap, false);


    std::cerr << "Used fields: " << ptMap->getUsedPositions() << std::endl;
    std::cerr << "adds: " << ptMap->iAddCount << std::endl;
    std::cerr << "add calls: " << ptMap->iAddKmerCount << std::endl;
    std::cerr << "aborts calls: " << ptMap->iAborts << std::endl;

    return 1;


    CLParser oParser(argc, argv);

    oParser.setArgument("k", std::to_string(10));
    oParser.setArgument("fastq", "/mnt/c/ownCloud/data/tsx/small.fq");

    FASTXreader<FASTQEntry>* pReader = FASTXreader<FASTQEntry>::createFQReader(&oParser);
    size_t iK = oParser.isSet("k") ? (size_t) oParser.getIntArgument("k") : 15;

    TSXHashMapTSX* pMap = new TSXHashMapTSX(16, 6, iK);

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
