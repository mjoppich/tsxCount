
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMapTSXSmall.h>
#include "testExecution.h"
#include <stddef.h>
#include <immintrin.h>
#include <tsxcount/TSXHashMapTSXPerf.h>


int main(int argc, char *argv[])
{


    uint32_t itK = 14;
    uint32_t iThreads = 4;

    TSXHashMapTSXPerf* ptMap = new TSXHashMapTSXPerf(26, 4, itK, iThreads);


    UBigInt testkmer = UBigInt::fromString("1010010110000010111000001110", ptMap->getMemoryPool());
    TSX::tsx_key_t kmerkey = ptMap->getHashingFunc()->apply(testkmer);

    uint64_t iPos = ptMap->getPosition(kmerkey, 1);
    std::cout << "position " << iPos << std::endl;

    iPos = ptMap->getPosition(kmerkey, 2);
    std::cout << "position " << iPos << std::endl;

    iPos = ptMap->getPosition(kmerkey, 3);
    std::cout << "position " << iPos << std::endl;

    iPos = ptMap->getPosition(kmerkey, 4);
    std::cout << "position " << iPos << std::endl;


    UBigInt oReprobe = ptMap->makeOverflowReprobe(2,3);
    std::cout << oReprobe.to_string() << std::endl;


    //testHashMapOld(ptMap, true);
    testHashMap(ptMap, iThreads != 1);


    std::cerr << "Used fields: " << ptMap->getUsedPositions() << std::endl;
    std::cerr << "adds: " << ptMap->iAddCount << std::endl;
    std::cerr << "add calls: " << ptMap->iAddKmerCount << std::endl;
    std::cerr << "aborts calls: " << ptMap->iAborts << std::endl;
    std::cerr << "general abort calls: " << ptMap->iTotalAborts << std::endl;

    return 1;


}
