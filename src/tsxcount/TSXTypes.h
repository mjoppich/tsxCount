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
}

#endif //TSXCOUNT_TSXTYPES_H
