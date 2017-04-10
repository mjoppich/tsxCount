//
// Created by mjopp on 07.02.2017.
//

#include "TSXTypes.h"
#include <sstream>

std::string TSX::print_kmer_t(TSX::tsx_kmer_t &elem) {
    std::stringstream oSS;

    return oSS.str();
}

std::string TSX::print_key_t(TSX::tsx_key_t& elem)
{
    std::stringstream oSS;

    oSS << "key: " << elem.to_string() << " value: ";

    return oSS.str();
}

std::string TSX::print_val_t(TSX::tsx_val_t& elem, uint32_t ikeybits)
{
    std::stringstream oSS;

    std::string keystring(ikeybits, '-');

    oSS << "key: " << keystring << " value: " << elem.to_string();

    return oSS.str();
}

std::string TSX::print_keyval_t(TSX::tsx_keyval_t& elem, uint32_t iStorageBits)
{
    std::stringstream oSS;

    TSX::tsx_key_t key = elem >> iStorageBits;
    TSX::tsx_val_t val = elem.mod2( iStorageBits+1 );

    oSS << "key: " << key.to_string() << " value: " << val.to_string();

    return oSS.str();
}
