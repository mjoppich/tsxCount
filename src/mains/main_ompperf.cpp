
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMapOMPPerf.h>
#include "testExecution.h"

int main(int argc, char *argv[])
{

    uint32_t itK = 14;
    TSXHashMapOMPPerf* pTMap = new TSXHashMapOMPPerf(26, 4, itK);
    testHashMap(pTMap, true);

    return 0;

}
