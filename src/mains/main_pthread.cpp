
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

    uint32_t itK = 14;
    TSXHashMapPThread* ptMap = new TSXHashMapPThread(26, 4, itK, 2);
    
    testHashMap(ptMap, true);
    
    return 1;


}
