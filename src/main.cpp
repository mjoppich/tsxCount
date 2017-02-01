
#include <iostream>


#include <tsxcount/TSXHashMap.h>

int main(int argc, char *argv[])
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

    TSXHashMap* pMap = new TSXHashMap(12, 12, 8);

    pMap->addKmer(0);
    pMap->addKmer(1);


    return 0;
}
