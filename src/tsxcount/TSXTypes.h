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
    typedef UBigInt tsx_key_t;
    typedef UBigInt tsx_val_t;
    typedef UBigInt tsx_keyval_t;
    typedef uint32_t tsx_reprobe_t;
}

#endif //TSXCOUNT_TSXTYPES_H
