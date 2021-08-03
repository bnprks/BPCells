#include "bp128.h"


namespace BPCells {

// template <unsigned B>
// void print_vec(vec i) {
//     uint32_t array[4];
//     store((vec *) &array, i);
//     std::cout << std::bitset<B>(array[0]) << '\n';
//     std::cout << std::bitset<B>(array[1]) << '\n';
//     std::cout << std::bitset<B>(array[2]) << '\n';
//     std::cout << std::bitset<B>(array[3]) << '\n';
// }

// void print_ints(vec i) {
//     uint32_t array[4];
//     store((vec *) &array, i);
//     std::cout << array[0] << '\t';
//     std::cout << array[1] << '\t';
//     std::cout << array[2] << '\t';
//     std::cout << array[3] << '\n';
// }


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

template <unsigned B>
vec unpackd1(vec initOffset, const vec *in, vec *out) {
    unsigned int shift;
    vec InReg, OutReg;
    vec mask = (B == 32) ? splat(0xffffffff) : splat((1U << B) - 1);
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
        OutReg = prefixSum(OutReg, initOffset);
        initOffset = OutReg;
        store(out++, OutReg);
    })
    return initOffset;
}

template<>
vec unpackd1<0>(vec initOffset, const vec *in, vec *out) {
    int i;
    BP128_UNROLL_LOOP32(i, {
        store(out++, initOffset);
    })
    return initOffset;
}

template<>
vec unpackd1<32>(vec initOffset, const vec *in, vec *out) {
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

template <unsigned B, bool MASK>
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
void packd1<32, true>(vec initOffset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    unpack<32> (in, out);
}

template <>
void packd1<32, false>(vec initOffset, const vec *in, vec *out) {
    // To save a bit of computation time, we just do straight memory copy
    // on 32-bit packing
    unpack<32> (in, out);
}

// Temporarily removed since there's no need for unsafe operations when
// packing is already plenty fast enough. Unpacking speed is the main bottleneck.
// template <unsigned B>
// void packd1_nomask(vec initOffset, const vec *in, vec *out) {
//     packd1<B, false>(initOffset, in, out);
// }

template <unsigned B>
void packd1_mask(vec initOffset, const vec *in, vec *out) {
    packd1<B, true>(initOffset, in, out);
}

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


// Find maximum number of bits required to represent input in FOR encoding.
// minvalue and bits are return values, corresponding to the number of bits
// required when the frame of reference is set to minvalue
void simdmaxbitsFOR(const uint32_t *in, uint32_t & bits, uint32_t & minvalue) {
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
        unpackd1<0>(init, _in, _out);
        break;
        BP128_SWITCH_CASE32(unpackd1, init, _in, _out)
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



/************* TESTING FUNCTIONS FOR BITPACKING *************/

uint32_t random_uint32_t() {
    return (
        (uint32_t) (rand() & 255) +
        (uint32_t) ((rand() & 255) << 8) +
        (uint32_t) ((rand() & 255) << 16) +
        (uint32_t) ((rand() & 255) << 24)
    );
}
void fill_buf(uint32_t *buf, int bits) {
    uint32_t mask = (uint32_t) ((UINT64_C(1) << bits) - 1);
    
    for (int i = 0; i < 128; i++) {
        buf[i] = random_uint32_t() & mask;
    }
}
bool equal_buf(uint32_t *buf1, uint32_t *buf2) {
    for (int i = 0; i < 128; i++) {
        if (buf1[i] != buf2[i]) return false;
    }
    return true;
}

bool test_bitpacking() {
    
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t output_buf[128];

    // Test vanilla bitpacking
    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        if (equal_buf(input_buf, output_buf)) {
            printf("Input and output aren't different before test\n");
            return false;
        }
        int simd_bits = simdmaxbits(input_buf);
        if (simd_bits != bits) {
            printf("Wrong number of simdmaxbits\n");
            return false;
        }
        simdpack(input_buf, packed_buf, simd_bits);
        simdunpack(packed_buf, output_buf, simd_bits);
        if (!equal_buf(input_buf, output_buf)) {
            printf("Input buffer doesn't match output buffer\n");
            return false;
        }
    }
    
    // Test d1 bitpacking
    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        uint32_t start_val = rand() & ((2<<12) - 1);
        input_buf[0] += start_val;
        for (int i = 1; i < 128; i++)
            input_buf[i] += input_buf[i-1];

        if (equal_buf(input_buf, output_buf)) {
            printf("Input and output aren't different before test\n");
            return false;
        }
        int simd_bits = simdmaxbitsd1(start_val, input_buf);
        if (simd_bits != bits) {
            printf("Wrong number of simdmaxbits\n");
            return false;
        }
        simdpackd1(start_val, input_buf, packed_buf, simd_bits);
        simdunpackd1(start_val, packed_buf, output_buf, simd_bits);
        if (!equal_buf(input_buf, output_buf)) {
            printf("Input buffer doesn't match output buffer\n");
            return false;
        }
    }

    // Test FOR bitpacking
    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        uint32_t start_val = rand() & ((2<<12) - 1);
        
        for (int i = 0; i < 128; i++)
            input_buf[i] += start_val;

        if (equal_buf(input_buf, output_buf)) {
            printf("Input and output aren't different before test\n");
            return false;
        }
        uint32_t simd_bits, min_value;
        simdmaxbitsFOR(input_buf, simd_bits, min_value);
        uint32_t real_min = INT32_MAX;
        for (int i = 0; i < 128; i++) {
            real_min = std::min(real_min, input_buf[i]);
        }
        if (real_min != min_value) {
            printf("Wrong minimum value\n");
            return false;
        }
        if (simd_bits != bits) {
            printf("Wrong number of simdmaxbits\n");
            return false;
        }
        simdpackFOR(start_val, input_buf, packed_buf, simd_bits);
        simdunpackFOR(start_val, packed_buf, output_buf, simd_bits);
        if (!equal_buf(input_buf, output_buf)) {
            printf("Input buffer doesn't match output buffer\n");
            return false;
        }
    }
    return true;
}

} // end namespace BPCells