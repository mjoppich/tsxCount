
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>

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

    TSXHashMap* pMap = new TSXHashMap(3, 4, 4);

    UBigInt oKmer1(0);
    UBigInt oKmer2(1);

    oKmer1.resize(8);
    oKmer2.resize(8);


    /*
    oTest = ~oTest;
    std::cerr << oTest.to_string() << std::endl;
    std::cerr << (oTest == 0) << std::endl;
*/
/*
    pMap->addKmer( oKmer1 );
    pMap->addKmer( oKmer2 );
    pMap->addKmer( oKmer1 );
*/
    size_t i;
    for ( i = 0; i < 8*256; ++i)
    {
        std::cerr << "adding kmer: " << oKmer1.to_string() << " " << std::to_string(i) << std::endl;
        pMap->addKmer( oKmer1 );

        std::cerr << "adding kmer: " << oKmer2.to_string() << " " << std::to_string(i) << std::endl;
        pMap->addKmer( oKmer2 );
    }

    UBigInt oRes1 = pMap->getKmerCount(oKmer1);
    UBigInt oRes2 = pMap->getKmerCount(oKmer2);

    std::cerr << "Added each kmer: " << std::to_string(i) << std::endl;

    std::cerr << "kmer: " << oKmer1.to_string() << ": " << oRes1.to_string() << " " << std::to_string(oRes1.toUInt()) << std::endl;
    std::cerr << "kmer: " << oKmer2.to_string() << ": " << oRes2.to_string() << " " << std::to_string(oRes2.toUInt()) << std::endl;

    return 0;
}
