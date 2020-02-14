
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
    uint32_t iThreads = 1;

    TSXHashMapTSXPerf* ptMap = new TSXHashMapTSXPerf(26, 4, itK, iThreads);

    //testHashMapOld(ptMap, true);
    testHashMap(ptMap, iThreads != 1);


    std::cerr << "Used fields: " << ptMap->getUsedPositions() << std::endl;
    std::cerr << "adds: " << ptMap->iAddCount << std::endl;
    std::cerr << "add calls: " << ptMap->iAddKmerCount << std::endl;
    std::cerr << "aborts calls: " << ptMap->iAborts << std::endl;
    std::cerr << "general abort calls: " << ptMap->iTotalAborts << std::endl;

    return 1;


}
