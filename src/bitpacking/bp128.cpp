#include "bp128.h"

// #include <iostream>
// void print_ints(BPCells::vec i) {
//     uint32_t array[4];
//     BPCells::store((BPCells::vec *) &array, i);
//     std::cout << array[0] << '\t';
//     std::cout << array[1] << '\t';
//     std::cout << array[2] << '\t';
//     std::cout << array[3] << '\n';
// }

// template <unsigned B>
// void print_vec(vec i) {
//     uint32_t array[4];
//     store((vec *) &array, i);
//     std::cout << std::bitset<B>(array[0]) << '\n';
//     std::cout << std::bitset<B>(array[1]) << '\n';
//     std::cout << std::bitset<B>(array[2]) << '\n';
//     std::cout << std::bitset<B>(array[3]) << '\n';
// }

namespace BPCells {






#define BP128_UNROLL_LOOP1(start, counter, body) \
    counter = start;                             \
    body

#define BP128_UNROLL_LOOP4(start, counter, body)   \
    BP128_UNROLL_LOOP1((start + 0), counter, body) \
    BP128_UNROLL_LOOP1((start + 1), counter, body) \
    BP128_UNROLL_LOOP1((start + 2), counter, body) \
    BP128_UNROLL_LOOP1((start + 3), counter, body)

#define BP128_UNROLL_LOOP32(counter, body) \
    BP128_UNROLL_LOOP4(0, counter, body)   \
    BP128_UNROLL_LOOP4(4, counter, body)   \
    BP128_UNROLL_LOOP4(8, counter, body)   \
    BP128_UNROLL_LOOP4(12, counter, body)  \
    BP128_UNROLL_LOOP4(16, counter, body)  \
    BP128_UNROLL_LOOP4(20, counter, body)  \
    BP128_UNROLL_LOOP4(24, counter, body)  \
    BP128_UNROLL_LOOP4(28, counter, body)

//##############################################################################
//#                         BP128 CODEC                                        #
//##############################################################################

template <unsigned B>
void unpack(const vec *in, vec *out) {
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        unsigned int shift = (i * B) & 31;
        if (shift == 0) {
            InReg = load(in++);
            OutReg = InReg;
        } else {
            OutReg = shift_r(InReg, shift);
        }

        if (shift > 32 - B) {
            InReg = load(in++);
            OutReg = bitwise_or(OutReg, shift_l(InReg, 32 - shift));
        }
        OutReg = bitwise_and(OutReg, mask);
        store(out++, OutReg);
    })
}

template<>
void unpack<0>(const vec *in, vec *out) {
    vec reg = splat(0);
    int i;
    BP128_UNROLL_LOOP32(i, {
        store(out++, reg);
    })
}

template <>
void unpack<32>(const vec *in, vec *out) {
    vec Reg;
    int i;
    BP128_UNROLL_LOOP32(i, {
        Reg = load(in++);
        store(out++, Reg);
    })
}


template <unsigned B, bool MASK>
void pack(const vec *in, vec *out) {
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        unsigned int shift = (i * B) & 31;
        InReg = load(in++);
        if (MASK) {
            InReg = bitwise_and(InReg, mask);
        }

        if (shift == 0) {
            OutReg = InReg;
        } else {
            OutReg = bitwise_or(OutReg, shift_l(InReg, shift));
        }

        if (shift >= 32 - B) {
            store(out++, OutReg);
            if (shift > 32 - B) {
                OutReg = shift_r(InReg, 32 - shift);
            }
        }
    })
}

template <>
void pack<32, true>(const vec *in, vec *out) {
    unpack<32> (in, out);
}

template <>
void pack<32, false>(const vec *in, vec *out) {
    unpack<32> (in, out);
}

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// template <unsigned B>
// void pack_nomask(const vec *in, vec *out) {
//     pack<B, false>(in, out);
// }

template <unsigned B>
void pack_mask(const vec *in, vec *out) {
    pack<B, true>(in, out);
}

//##############################################################################
//#                         BP128D1 CODEC                                      #
//##############################################################################

template <unsigned B, bool ZIGZAG>
vec unpackd1(vec initOffset, const vec *in, vec *out) {
    unsigned int shift;
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    vec zero = splat(0);
    vec one = splat(1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        shift = (i * B) & 31;
        if (shift == 0) {
            InReg = load(in++);
            OutReg = InReg;
        } else {
            OutReg = shift_r(InReg, shift);
        }

        if (shift > 32 - B) {
            InReg = load(in++);
            OutReg = bitwise_or(OutReg, shift_l(InReg, 32 - shift));
        }

        OutReg = bitwise_and(OutReg, mask);

        if (ZIGZAG) {
            // (i >>> 1) ^ -(i & 1)
            OutReg = bitwise_xor(
                shift_r(OutReg, 1), sub(zero, bitwise_and(OutReg, one))
            );
        }

        OutReg = prefixSum(OutReg, initOffset);

        initOffset = OutReg;
        store(out++, OutReg);
    })
    return initOffset;
}

template <>
vec unpackd1<0, false>(vec initOffset, const vec *in, vec *out) {
    int i;
    BP128_UNROLL_LOOP32(i, {
        store(out++, initOffset);
    })
    return initOffset;
}

template <>
vec unpackd1<0, true>(vec initOffset, const vec *in, vec *out) {
    return unpackd1<0, false>(initOffset, in, out);
}

template<>
vec unpackd1<32, false>(vec initOffset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    int i;
    vec Reg;
    BP128_UNROLL_LOOP32(i, {
        Reg = load(in++);
        store(out++, Reg);
    })
    return Reg;
}

template<>
vec unpackd1<32, true>(vec initOffset, const vec *in, vec *out) {
    return unpackd1<32, false>(initOffset, in, out);
}

template<unsigned B>
vec unpackd1_nozigzag(vec initOffset, const vec *in, vec *out) {
    return unpackd1<B, false>(initOffset, in, out);
}

template<unsigned B>
vec unpackd1_zigzag(vec initOffset, const vec *in, vec *out) {
    return unpackd1<B, true>(initOffset, in, out);
}

template <unsigned B, bool MASK, bool ZIGZAG>
void packd1(vec initOffset, const vec *in, vec *out) {
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        unsigned int shift = (i * B) & 31;
        InReg = load(in++);
        
        vec _tmp1 = delta(InReg, initOffset);
        initOffset = InReg;
        InReg = _tmp1;
        
        if (ZIGZAG) {
            // (i >> 31) ^ (i << 1)
            InReg = bitwise_xor(shift_r_arith(InReg, 31), shift_l(InReg, 1));
        }

        if (MASK) {
            InReg = bitwise_and(InReg, mask);
        }

        if (shift == 0) {
            OutReg = InReg;
        } else {
            OutReg = bitwise_or(OutReg, shift_l(InReg, shift));
        }

        if (shift >= 32 - B) {
            store(out++, OutReg);
            if (shift > 32 - B) {
                OutReg = shift_r(InReg, 32 - shift);
            }
        }
    })
}

template <>
void packd1<32, true, true>(vec initOffset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    unpack<32> (in, out);
}
template <> void packd1<32, true, false>(vec initOffset, const vec *in, vec *out) {unpack<32> (in, out);}
template <> void packd1<32, false, true>(vec initOffset, const vec *in, vec *out) {unpack<32> (in, out);}
template <> void packd1<32, false, false>(vec initOffset, const vec *in, vec *out) {unpack<32> (in, out);}

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// template <unsigned B>
// void packd1_nomask(vec initOffset, const vec *in, vec *out) {
//     packd1<B, false>(initOffset, in, out);
// }

template <unsigned B>
void packd1_mask(vec initOffset, const vec *in, vec *out) {
    packd1<B, true, false>(initOffset, in, out);
}

template <unsigned B>
void packd1z_mask(vec initOffset, const vec *in, vec *out) {
    packd1<B, true, true>(initOffset, in, out);
}

//##############################################################################
//#                         BP128FOR CODEC                                     #
//##############################################################################

template <unsigned B>
void unpackFOR(vec initOffset, const vec *in, vec *out) {
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        unsigned int shift = (i * B) & 31;
        if (shift == 0) {
            InReg = load(in++);
            OutReg = InReg;
        } else {
            OutReg = shift_r(InReg, shift);
        }

        if (shift > 32 - B) {
            InReg = load(in++);
            OutReg = bitwise_or(OutReg, shift_l(InReg, 32 - shift));
        }
        OutReg = bitwise_and(OutReg, mask);
        OutReg = add(OutReg, initOffset);
        store(out++, OutReg);
    })
}

template<>
void unpackFOR<0>(vec initOffset, const vec *in, vec *out) {
    int i;
    BP128_UNROLL_LOOP32(i, {
        store(out++, initOffset);
    })
}

template<>
void unpackFOR<32>(vec initOffset, const vec *in, vec *out) {
    int i;
    vec Reg;
    BP128_UNROLL_LOOP32(i, {
        Reg = load(in++);
        store(out++, Reg);
    })
}

template <unsigned B, bool MASK>
void packFOR(uint32_t offset, const vec *in, vec *out) {
    vec InReg, OutReg;
    vec initOffset = splat(offset);
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
    int i;
    BP128_UNROLL_LOOP32(i, {
        unsigned int shift = (i * B) & 31;
        InReg = load(in++);
        InReg = sub(InReg, initOffset);

        if (MASK) {
            InReg = bitwise_and(InReg, mask);
        }

        if (shift == 0) {
            OutReg = InReg;
        } else {
            OutReg = bitwise_or(OutReg, shift_l(InReg, shift));
        }

        if (shift >= 32 - B) {
            store(out++, OutReg);
            if (shift > 32 - B) {
                OutReg = shift_r(InReg, 32 - shift);
            }
        }
    })
}

template <>
void packFOR<32, true>(uint32_t offset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    unpack<32> (in, out);
}

template <>
void packFOR<32, false>(uint32_t offset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    unpack<32> (in, out);
}

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// template <unsigned B>
// void packFOR_nomask(uint32_t initOffset, const vec *in, vec *out) {
//     packFOR<B, false>(initOffset, in, out);
// }

template <unsigned B>
void packFOR_mask(uint32_t initOffset, const vec *in, vec *out) {
    packFOR<B, true>(initOffset, in, out);
}



// Find maximum number of bits required to represent input in d1 encoding
uint32_t simdmaxbits(const uint32_t *in) {
    const vec *_in = (vec *)in;
    vec InReg;
    vec accumulator = splat(0);
    int i;

    BP128_UNROLL_LOOP32(i, {
        InReg = load(_in++);
        accumulator = bitwise_or(accumulator, InReg);
    })

    return maxbitas32int(accumulator);
}

// Find maximum number of bits required to represent input in d1 encoding
uint32_t simdmaxbitsd1(uint32_t initvalue, const uint32_t *in) {
    const vec *_in = (vec *)in;
    vec InReg;
    vec initOffset = splat(initvalue);
    vec accumulator = splat(0);
    int i;

    BP128_UNROLL_LOOP32(i, {
        InReg = load(_in++);
        accumulator = bitwise_or(accumulator, delta(InReg, initOffset));
        initOffset = InReg;
    })

    return maxbitas32int(accumulator);
}

// Find maximum number of bits required to represent input in d1z encoding
uint32_t simdmaxbitsd1z(uint32_t initvalue, const uint32_t *in) {
    const vec *_in = (vec *)in;
    vec InReg, tmp;
    vec initOffset = splat(initvalue);
    vec accumulator = splat(0);
    int i;

    BP128_UNROLL_LOOP32(i, {
        InReg = load(_in++);
        tmp = delta(InReg, initOffset);
        initOffset = InReg;
        InReg = tmp;

        InReg = bitwise_xor(shift_r_arith(InReg, 31), shift_l(InReg, 1));
        //print_ints(InReg);
        accumulator = bitwise_or(accumulator, InReg);
    })

    return maxbitas32int(accumulator);
}

// Find maximum number of bits required to represent input in FOR encoding.
// minvalue and bits are return values, corresponding to the number of bits
// required when the frame of reference is set to minvalue
void simdmaxbitsFORwithmin(const uint32_t *in, uint32_t & bits, uint32_t & minvalue) {
    const vec *_in = (vec *)in;
    vec InReg;
    vec mins = splat(INT32_MAX);
    vec maxs = splat(0);
    int i;

    BP128_UNROLL_LOOP32(i, {
        InReg = load(_in++);
        mins = min(mins, InReg);
        maxs = max(maxs, InReg);
    })
    maxdiffbitsas32int(mins, maxs, bits, minvalue);
}

// Find maximum number of bits required to represent input in FOR encoding,
// given a known minvalue.
uint32_t simdmaxbitsFOR(const uint32_t minvalue, const uint32_t *in) {
    const vec *_in = (vec *)in;
    vec InReg;
    vec min = splat(minvalue);
    vec accumulator = splat(0);
    int i;

    BP128_UNROLL_LOOP32(i, {
        InReg = load(_in++);
        accumulator = bitwise_or(accumulator, sub(InReg, min));
    })
    return maxbitas32int(accumulator);
}



#define BP128_SWITCH_CASE(i, function, ...) \
    case i:                                 \
        function<i>(__VA_ARGS__);           \
        break;
#define BP128_SWITCH_CASE4(start, function, ...)          \
    BP128_SWITCH_CASE((start + 0), function, __VA_ARGS__) \
    BP128_SWITCH_CASE((start + 1), function, __VA_ARGS__) \
    BP128_SWITCH_CASE((start + 2), function, __VA_ARGS__) \
    BP128_SWITCH_CASE((start + 3), function, __VA_ARGS__)

#define BP128_SWITCH_CASE32(function, ...)        \
    BP128_SWITCH_CASE4(1, function, __VA_ARGS__)  \
    BP128_SWITCH_CASE4(5, function, __VA_ARGS__)  \
    BP128_SWITCH_CASE4(9, function, __VA_ARGS__)  \
    BP128_SWITCH_CASE4(13, function, __VA_ARGS__) \
    BP128_SWITCH_CASE4(17, function, __VA_ARGS__) \
    BP128_SWITCH_CASE4(21, function, __VA_ARGS__) \
    BP128_SWITCH_CASE4(25, function, __VA_ARGS__) \
    BP128_SWITCH_CASE4(29, function, __VA_ARGS__)


/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out".
 * The input values are masked to be less than 1<<bit. */
void simdpack(const uint32_t *in, uint32_t *out, const uint32_t bit) {
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
        BP128_SWITCH_CASE32(pack_mask, _in, _out)
    }
}

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out".
 * The input values are assumed to be less than 1<<bit. */
// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmask(const uint32_t *in, uint32_t *out, const uint32_t bit) {
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(pack_nomask, _in, _out)
//     }
// }

/* reads  "bit" 128-bit vectors from "in", writes  128 values to "out" */
void simdunpack(const uint32_t *in, uint32_t *out, const uint32_t bit) {
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
    case 0:
        unpack<0>(_in, _out);
        break;
        BP128_SWITCH_CASE32(unpack, _in, _out)
    }
}

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out"
   integer values should be in sorted order (for best results).
   The differences are masked so that only the least significant "bit" bits are
   used. */
void simdpackd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                   const uint32_t bit) {
    vec init = splat(initvalue);
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
        BP128_SWITCH_CASE32(packd1_mask, init, _in, _out)
    }
}

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out"
   integer values should be in sorted order (for best results).
   The difference values are assumed to be less than 1<<bit. */
// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmaskd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
//                            const uint32_t bit) {
//     vec init = splat(initvalue);
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(packd1_nomask, init, _in, _out)
//     }
// }

/* reads "bit" 128-bit vectors from "in", writes  128 values to "out" */
void simdunpackd1(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                     const uint32_t bit) {
    vec init = splat(initvalue);
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
    case 0:
        unpackd1<0, false>(init, _in, _out);
        break;
        BP128_SWITCH_CASE32(unpackd1_nozigzag, init, _in, _out)
    }
}

/* reads 128 values from "in", writes  "bit" 128-bit vectors to "out"
   integer values should be in nearly sorted order (for best results).
   The values are zigzag encoded, then masked so that only the least significant
   "bit" bits are used. 
   ZigZag encoding references: https://developers.google.com/protocol-buffers/docs/encoding?csw=1#signed-ints
   https://gist.github.com/lemire/b6437fbd193395d8e4ccac1a5b2e50cc*/
void simdpackd1z(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                   const uint32_t bit) {
    vec init = splat(initvalue);
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
        BP128_SWITCH_CASE32(packd1z_mask, init, _in, _out)
    }
}

/* reads "bit" 128-bit vectors from "in", writes  128 values to "out" */
void simdunpackd1z(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                     const uint32_t bit) {
    vec init = splat(initvalue);
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
    case 0:
        unpackd1<0, true>(init, _in, _out);
        break;
        BP128_SWITCH_CASE32(unpackd1_zigzag, init, _in, _out)
    }
}

void simdpackFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                    const uint32_t bit) {
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
        BP128_SWITCH_CASE32(packFOR_mask, initvalue, _in, _out)
    }
}

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// void simdpackwithoutmaskFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
//                     const uint32_t bit) {
//     const vec *_in = (vec *)in;
//     vec *_out = (vec *)out;
//     switch (bit) {
//         BP128_SWITCH_CASE32(packFOR_nomask, initvalue, _in, _out)
//     }
// }

void simdunpackFOR(uint32_t initvalue, const uint32_t *in, uint32_t *out,
                    const uint32_t bit) {
    vec init = splat(initvalue);
    const vec *_in = (vec *)in;
    vec *_out = (vec *)out;
    switch (bit) {
    case 0:
        unpackFOR<0>(init, _in, _out);
        break;
        BP128_SWITCH_CASE32(unpackFOR, init, _in, _out)
    }
}

} // end namespace BPCells