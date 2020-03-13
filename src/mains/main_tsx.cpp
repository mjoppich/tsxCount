
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


    uint32_t itK = 14;
    uint32_t iThreads = 4;

    TSXHashMapTSX* ptMap = new TSXHashMapTSX(26, 4, itK, iThreads);

    testHashMap(ptMap, true);


    std::cerr << "Used fields: " << ptMap->getUsedPositions() << std::endl;
    std::cerr << "adds: " << ptMap->iAddCount << std::endl;
    std::cerr << "add calls: " << ptMap->iAddKmerCount << std::endl;
    std::cerr << "aborts calls: " << ptMap->iAborts << std::endl;

    return 1;


}
