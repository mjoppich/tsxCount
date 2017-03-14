//
// Created by mjopp on 14.06.2016.
//

#ifndef TSXCOUNT_TSXTYPES_H
#define TSXCOUNT_TSXTYPES_H

#include <stdint.h>
#include <src/tsxutils/UBigInt.h>


/**
 * \namespace TSX
 *
 *
 * \brief contains the typedefs for tsx variables
 *
 * This namespace contains variable type definitions
 *
 * \note Performance check this
 *
 * \author Markus Joppich
 *
 */
namespace TSX
{
    typedef UBigInt tsx_kmer_t;

    // 2k bits
    typedef UBigInt tsx_key_t;

    // storage bits
    typedef UBigInt tsx_val_t;

    // key + val
    typedef UBigInt tsx_keyval_t;

    // lower l bits of key
    typedef UBigInt tsx_reprobe_t;

    // higher 2k-l bits of key
    typedef UBigInt tsx_func_t;

    std::string print_kmer_t(TSX::tsx_kmer_t& elem);
    std::string print_key_t(TSX::tsx_key_t& elem);
    std::string print_val_t(TSX::tsx_val_t& elem, uint32_t ikeybits = 32);
    std::string print_keyval_t(TSX::tsx_keyval_t& elem, uint32_t iStorageBits = 3);

}

#endif //TSXCOUNT_TSXTYPES_H
