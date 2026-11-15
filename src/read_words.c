////////////////////////////////////////////////////////////////////////////
//                           **** WAVPACK ****                            //
//                  Hybrid Lossless Wavefile Compressor                   //
//              Copyright (c) 1998 - 2013 Conifer Software.               //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// read_words.c

// This module provides entropy word decoding functions using
// a variation on the Rice method.  This was introduced in version 3.93
// because it allows splitting the data into a "lossy" stream and a
// "correction" stream in a very efficient manner and is therefore ideal
// for the "hybrid" mode.  For 4.0, the efficiency of this method was
// significantly improved by moving away from the normal Rice restriction of
// using powers of two for the modulus divisions and now the method can be
// used for both hybrid and pure lossless encoding.

// Samples are divided by median probabilities at 5/7 (71.43%), 10/49 (20.41%),
// and 20/343 (5.83%). Each zone has 3.5 times fewer samples than the
// previous. Using standard Rice coding on this data would result in 1.4
// bits per sample average (not counting sign bit). However, there is a
// very simple encoding that is over 99% efficient with this data and
// results in about 1.22 bits per sample.

#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#include "wavpack_local.h"

#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#define USE_FAST_WORDS
#define get_words_lossless get_words_lossless_reference
#endif

#if defined (HAVE___BUILTIN_CTZ) || defined (_WIN64)
#define USE_CTZ_OPTIMIZATION    // use ctz intrinsic (or Windows equivalent) to count trailing ones
#elif defined(__WATCOMC__) && defined(__386__)
#define USE_CTZ_OPTIMIZATION
extern __inline uint32_t _bsf_watcom (uint32_t);
#pragma aux _bsf_watcom = \
  "bsf eax, eax" \
  parm [eax] nomemory \
  value [eax] \
  modify exact [eax] nomemory;
#else
#define USE_NEXT8_OPTIMIZATION  // optimization using a table to count trailing ones
#endif

#define USE_BITMASK_TABLES      // use tables instead of shifting for certain masking operations

///////////////////////////// local table storage ////////////////////////////

#ifdef USE_NEXT8_OPTIMIZATION
static const char ones_count_table [] = {
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,8
};
#endif

///////////////////////////// executable code ////////////////////////////////

static uint32_t __inline read_code (Bitstream *bs, uint32_t maxcode);

// Read the next word from the bitstream "wvbits" and return the value. This
// function can be used for hybrid or lossless streams, but since an
// optimized version is available for lossless this function would normally
// be used for hybrid only. If a hybrid lossless stream is being read then
// the "correction" offset is written at the specified pointer. A return value
// of WORD_EOF indicates that the end of the bitstream was reached (all 1s) or
// some other error occurred.

int32_t FASTCALL get_word (WavpackStream *wps, int chan, int32_t *correction)
{
    struct entropy_data *c = wps->w.c + chan;
    uint32_t ones_count, low, mid, high;
    int32_t value;
    int sign;

    if (!wps->wvbits.ptr)
        return WORD_EOF;

    if (correction)
        *correction = 0;

    if (!(wps->w.c [0].median [0] & ~1) && !wps->w.holding_zero && !wps->w.holding_one && !(wps->w.c [1].median [0] & ~1)) {
        uint32_t mask;
        int cbits;

        if (wps->w.zeros_acc) {
            if (--wps->w.zeros_acc) {
                c->slow_level -= (c->slow_level + SLO) >> SLS;
                return 0;
            }
        }
        else {
            for (cbits = 0; cbits < 33 && getbit (&wps->wvbits); ++cbits);

            if (cbits == 33)
                return WORD_EOF;

            if (cbits < 2)
                wps->w.zeros_acc = cbits;
            else {
                for (mask = 1, wps->w.zeros_acc = 0; --cbits; mask <<= 1)
                    if (getbit (&wps->wvbits))
                        wps->w.zeros_acc |= mask;

                wps->w.zeros_acc |= mask;
            }

            if (wps->w.zeros_acc) {
                c->slow_level -= (c->slow_level + SLO) >> SLS;
                CLEARA (wps->w.c [0].median);
                CLEARA (wps->w.c [1].median);
                return 0;
            }
        }
    }

    if (wps->w.holding_zero)
        ones_count = wps->w.holding_zero = 0;
    else {
#ifdef USE_CTZ_OPTIMIZATION
        while (wps->wvbits.bc < LIMIT_ONES) {
            if (++(wps->wvbits.ptr) == wps->wvbits.end)
                wps->wvbits.wrap (&wps->wvbits);

            wps->wvbits.sr |= *(wps->wvbits.ptr) << wps->wvbits.bc;
            wps->wvbits.bc += sizeof (*(wps->wvbits.ptr)) * 8;
        }

#ifdef _MSC_VER
        { unsigned long res; _BitScanForward (&res, (unsigned long)~wps->wvbits.sr); ones_count = (uint32_t) res; }
#elif defined(__WATCOMC__) && defined(__386__)
        ones_count = _bsf_watcom (~wps->wvbits.sr);
#else
        ones_count = __builtin_ctz (~wps->wvbits.sr);
#endif

        if (ones_count >= LIMIT_ONES) {
            wps->wvbits.bc -= ones_count;
            wps->wvbits.sr >>= ones_count;

            for (; ones_count < (LIMIT_ONES + 1) && getbit (&wps->wvbits); ++ones_count);

            if (ones_count == (LIMIT_ONES + 1))
                return WORD_EOF;

            if (ones_count == LIMIT_ONES) {
                uint32_t mask;
                int cbits;

                for (cbits = 0; cbits < 33 && getbit (&wps->wvbits); ++cbits);

                if (cbits == 33)
                    return WORD_EOF;

                if (cbits < 2)
                    ones_count = cbits;
                else {
                    for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                        if (getbit (&wps->wvbits))
                            ones_count |= mask;

                    ones_count |= mask;
                }

                ones_count += LIMIT_ONES;
            }
        }
        else {
            wps->wvbits.bc -= ones_count + 1;
            wps->wvbits.sr >>= ones_count + 1;
        }
#elif defined (USE_NEXT8_OPTIMIZATION)
        int next8;

        if (wps->wvbits.bc < 8) {
            if (++(wps->wvbits.ptr) == wps->wvbits.end)
                wps->wvbits.wrap (&wps->wvbits);

            next8 = (wps->wvbits.sr |= *(wps->wvbits.ptr) << wps->wvbits.bc) & 0xff;
            wps->wvbits.bc += sizeof (*(wps->wvbits.ptr)) * 8;
        }
        else
            next8 = wps->wvbits.sr & 0xff;

        if (next8 == 0xff) {
            wps->wvbits.bc -= 8;
            wps->wvbits.sr >>= 8;

            for (ones_count = 8; ones_count < (LIMIT_ONES + 1) && getbit (&wps->wvbits); ++ones_count);

            if (ones_count == (LIMIT_ONES + 1))
                return WORD_EOF;

            if (ones_count == LIMIT_ONES) {
                uint32_t mask;
                int cbits;

                for (cbits = 0; cbits < 33 && getbit (&wps->wvbits); ++cbits);

                if (cbits == 33)
                    return WORD_EOF;

                if (cbits < 2)
                    ones_count = cbits;
                else {
                    for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                        if (getbit (&wps->wvbits))
                            ones_count |= mask;

                    ones_count |= mask;
                }

                ones_count += LIMIT_ONES;
            }
        }
        else {
            wps->wvbits.bc -= (ones_count = ones_count_table [next8]) + 1;
            wps->wvbits.sr >>= ones_count + 1;
        }
#else
        for (ones_count = 0; ones_count < (LIMIT_ONES + 1) && getbit (&wps->wvbits); ++ones_count);

        if (ones_count >= LIMIT_ONES) {
            uint32_t mask;
            int cbits;

            if (ones_count == (LIMIT_ONES + 1))
                return WORD_EOF;

            for (cbits = 0; cbits < 33 && getbit (&wps->wvbits); ++cbits);

            if (cbits == 33)
                return WORD_EOF;

            if (cbits < 2)
                ones_count = cbits;
            else {
                for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                    if (getbit (&wps->wvbits))
                        ones_count |= mask;

                ones_count |= mask;
            }

            ones_count += LIMIT_ONES;
        }
#endif

        if (wps->w.holding_one) {
            wps->w.holding_one = ones_count & 1;
            ones_count = (ones_count >> 1) + 1;
        }
        else {
            wps->w.holding_one = ones_count & 1;
            ones_count >>= 1;
        }

        wps->w.holding_zero = ~wps->w.holding_one & 1;
    }

    if ((wps->wphdr.flags & HYBRID_FLAG) && !chan)
        update_error_limit (wps);

    if (ones_count == 0) {
        low = 0;
        high = GET_MED (0) - 1;
        DEC_MED0 ();
    }
    else {
        low = GET_MED (0);
        INC_MED0 ();

        if (ones_count == 1) {
            high = low + GET_MED (1) - 1;
            DEC_MED1 ();
        }
        else {
            low += GET_MED (1);
            INC_MED1 ();

            if (ones_count == 2) {
                high = low + GET_MED (2) - 1;
                DEC_MED2 ();
            }
            else {
                low += (ones_count - 2) * GET_MED (2);
                high = low + GET_MED (2) - 1;
                INC_MED2 ();
            }
        }
    }

    low &= 0x7fffffff;
    high &= 0x7fffffff;

    if (low > high)         // make sure high and low make sense
        high = low;

    mid = (high + low + 1) >> 1;

    if (!c->error_limit)
        mid = read_code (&wps->wvbits, high - low) + low;
    else while (high - low > c->error_limit) {
        if (getbit (&wps->wvbits))
            mid = (high + (low = mid) + 1) >> 1;
        else
            mid = ((high = mid - 1) + low + 1) >> 1;
    }

    sign = getbit (&wps->wvbits);

    if (bs_is_open (&wps->wvcbits) && c->error_limit) {
        value = read_code (&wps->wvcbits, high - low) + low;

        if (correction)
            *correction = sign ? (mid - value) : (value - mid);
    }

    if (wps->wphdr.flags & HYBRID_BITRATE) {
        c->slow_level -= (c->slow_level + SLO) >> SLS;
        c->slow_level += wp_log2 (mid);
    }

    return sign ? ~mid : mid;
}

// This is an optimized version of get_word() that is used for lossless only
// (error_limit == 0). Also, rather than obtaining a single sample, it can be
// used to obtain an entire buffer of either mono or stereo samples.

int32_t get_words_lossless (WavpackStream *wps, int32_t *buffer, int32_t nsamples)
{
    struct entropy_data *c = wps->w.c;
    uint32_t ones_count, low, high;
    Bitstream *bs = &wps->wvbits;
    int32_t csamples;
#ifdef USE_NEXT8_OPTIMIZATION
    int32_t next8;
#endif

    if (nsamples && !bs->ptr) {
        memset (buffer, 0, (wps->wphdr.flags & MONO_DATA) ? nsamples * 4 : nsamples * 8);
        return nsamples;
    }

    if (!(wps->wphdr.flags & MONO_DATA))
        nsamples *= 2;

    for (csamples = 0; csamples < nsamples; ++csamples) {
        if (!(wps->wphdr.flags & MONO_DATA))
            c = wps->w.c + (csamples & 1);

        if (wps->w.holding_zero) {
            wps->w.holding_zero = 0;
            low = read_code (bs, GET_MED (0) - 1);
            DEC_MED0 ();
            buffer [csamples] = (getbit (bs)) ? ~low : low;

            if (++csamples == nsamples)
                break;

            if (!(wps->wphdr.flags & MONO_DATA))
                c = wps->w.c + (csamples & 1);
        }

        if (wps->w.c [0].median [0] < 2 && !wps->w.holding_one && wps->w.c [1].median [0] < 2) {
            uint32_t mask;
            int cbits;

            if (wps->w.zeros_acc) {
                if (--wps->w.zeros_acc) {
                    buffer [csamples] = 0;
                    continue;
                }
            }
            else {
                for (cbits = 0; cbits < 33 && getbit (bs); ++cbits);

                if (cbits == 33)
                    break;

                if (cbits < 2)
                    wps->w.zeros_acc = cbits;
                else {
                    for (mask = 1, wps->w.zeros_acc = 0; --cbits; mask <<= 1)
                        if (getbit (bs))
                            wps->w.zeros_acc |= mask;

                    wps->w.zeros_acc |= mask;
                }

                if (wps->w.zeros_acc) {
                    CLEARA (wps->w.c [0].median);
                    CLEARA (wps->w.c [1].median);
                    buffer [csamples] = 0;
                    continue;
                }
            }
        }

#ifdef USE_CTZ_OPTIMIZATION
        while (bs->bc < LIMIT_ONES) {
            if (++(bs->ptr) == bs->end)
                bs->wrap (bs);

            bs->sr |= *(bs->ptr) << bs->bc;
            bs->bc += sizeof (*(bs->ptr)) * 8;
        }

#ifdef _MSC_VER
        { unsigned long res; _BitScanForward (&res, (unsigned long)~wps->wvbits.sr); ones_count = (uint32_t) res; }
#elif defined(__WATCOMC__) && defined(__386__)
        ones_count = _bsf_watcom (~wps->wvbits.sr);
#else
        ones_count = __builtin_ctz (~wps->wvbits.sr);
#endif

        if (ones_count >= LIMIT_ONES) {
            bs->bc -= ones_count;
            bs->sr >>= ones_count;

            for (; ones_count < (LIMIT_ONES + 1) && getbit (bs); ++ones_count);

            if (ones_count == (LIMIT_ONES + 1))
                break;

            if (ones_count == LIMIT_ONES) {
                uint32_t mask;
                int cbits;

                for (cbits = 0; cbits < 33 && getbit (bs); ++cbits);

                if (cbits == 33)
                    break;

                if (cbits < 2)
                    ones_count = cbits;
                else {
                    for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                        if (getbit (bs))
                            ones_count |= mask;

                    ones_count |= mask;
                }

                ones_count += LIMIT_ONES;
            }
        }
        else {
            bs->bc -= ones_count + 1;
            bs->sr >>= ones_count + 1;
        }
#elif defined (USE_NEXT8_OPTIMIZATION)
        if (bs->bc < 8) {
            if (++(bs->ptr) == bs->end)
                bs->wrap (bs);

            next8 = (bs->sr |= *(bs->ptr) << bs->bc) & 0xff;
            bs->bc += sizeof (*(bs->ptr)) * 8;
        }
        else
            next8 = bs->sr & 0xff;

        if (next8 == 0xff) {
            bs->bc -= 8;
            bs->sr >>= 8;

            for (ones_count = 8; ones_count < (LIMIT_ONES + 1) && getbit (bs); ++ones_count);

            if (ones_count == (LIMIT_ONES + 1))
                break;

            if (ones_count == LIMIT_ONES) {
                uint32_t mask;
                int cbits;

                for (cbits = 0; cbits < 33 && getbit (bs); ++cbits);

                if (cbits == 33)
                    break;

                if (cbits < 2)
                    ones_count = cbits;
                else {
                    for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                        if (getbit (bs))
                            ones_count |= mask;

                    ones_count |= mask;
                }

                ones_count += LIMIT_ONES;
            }
        }
        else {
            bs->bc -= (ones_count = ones_count_table [next8]) + 1;
            bs->sr >>= ones_count + 1;
        }
#else
        for (ones_count = 0; ones_count < (LIMIT_ONES + 1) && getbit (bs); ++ones_count);

        if (ones_count >= LIMIT_ONES) {
            uint32_t mask;
            int cbits;

            if (ones_count == (LIMIT_ONES + 1))
                break;

            for (cbits = 0; cbits < 33 && getbit (bs); ++cbits);

            if (cbits == 33)
                break;

            if (cbits < 2)
                ones_count = cbits;
            else {
                for (mask = 1, ones_count = 0; --cbits; mask <<= 1)
                    if (getbit (bs))
                        ones_count |= mask;

                ones_count |= mask;
            }

            ones_count += LIMIT_ONES;
        }
#endif

        low = wps->w.holding_one;
        wps->w.holding_one = ones_count & 1;
        wps->w.holding_zero = ~ones_count & 1;
        ones_count = (ones_count >> 1) + low;

        if (ones_count == 0) {
            low = 0;
            high = GET_MED (0) - 1;
            DEC_MED0 ();
        }
        else {
            low = GET_MED (0);
            INC_MED0 ();

            if (ones_count == 1) {
                high = low + GET_MED (1) - 1;
                DEC_MED1 ();
            }
            else {
                low += GET_MED (1);
                INC_MED1 ();

                if (ones_count == 2) {
                    high = low + GET_MED (2) - 1;
                    DEC_MED2 ();
                }
                else {
                    low += (ones_count - 2) * GET_MED (2);
                    high = low + GET_MED (2) - 1;
                    INC_MED2 ();
                }
            }
        }

        low += read_code (bs, high - low);
        buffer [csamples] = (getbit (bs)) ? ~low : low;
    }

    return (wps->wphdr.flags & MONO_DATA) ? csamples : (csamples / 2);
}

// Read a single unsigned value from the specified bitstream with a value
// from 0 to maxcode. If there are exactly a power of two number of possible
// codes then this will read a fixed number of bits; otherwise it reads the
// minimum number of bits and then determines whether another bit is needed
// to define the code.

static uint32_t __inline read_code (Bitstream *bs, uint32_t maxcode)
{
    unsigned long local_sr;
    uint32_t extras, code;
    int bitcount;

    if (maxcode < 2)
        return maxcode ? getbit (bs) : 0;

    bitcount = count_bits (maxcode);
#ifdef USE_BITMASK_TABLES
    extras = bitset [bitcount] - maxcode - 1;
#else
    extras = (1 << bitcount) - maxcode - 1;
#endif

    local_sr = bs->sr;

    while (bs->bc < bitcount) {
        if (++(bs->ptr) == bs->end)
            bs->wrap (bs);

        local_sr |= (long)*(bs->ptr) << bs->bc;
        bs->bc += sizeof (*(bs->ptr)) * 8;
    }

#ifdef USE_BITMASK_TABLES
    if ((code = local_sr & bitmask [bitcount - 1]) >= extras)
#else
    if ((code = local_sr & ((1 << (bitcount - 1)) - 1)) >= extras)
#endif
        code = (code << 1) - extras + ((local_sr >> (bitcount - 1)) & 1);
    else
        bitcount--;

    if (sizeof (local_sr) < 8 && bs->bc > sizeof (local_sr) * 8) {
        bs->bc -= bitcount;
        bs->sr = *(bs->ptr) >> (sizeof (*(bs->ptr)) * 8 - bs->bc);
    }
    else {
        bs->bc -= bitcount;
        bs->sr = local_sr >> bitcount;
    }

    return code;
}

#ifdef USE_FAST_WORDS

#undef get_words_lossless

#ifdef BITSTREAM_SHORTS
typedef uint16_t bitword;
#else
typedef unsigned char bitword;
#endif

// This version keeps the bitstream and entropy state in locals and decodes
// stereo pairs together. The bit-level sequence and all state transitions are
// identical to get_words_lossless_reference().

int32_t get_words_lossless (WavpackStream *wps, int32_t *buffer, int32_t nsamples)
{
    uint32_t ones_count, low, high;
    int32_t csamples = 0;
    int mono = (wps->wphdr.flags & MONO_DATA) != 0;
    int wrap_error = 0;

    if (nsamples && !wps->wvbits.ptr) {
        memset (buffer, 0, mono ? nsamples * 4 : nsamples * 8);
        return nsamples;
    }

    if (!mono)
        nsamples *= 2;

    bitword *ptr = wps->wvbits.ptr;
    bitword *end = wps->wvbits.end;
    uint64_t sr = wps->wvbits.sr;
    int bc = wps->wvbits.bc;

    uint32_t m0 = wps->w.c [0].median [0], m1 = wps->w.c [0].median [1], m2 = wps->w.c [0].median [2];
    uint32_t s0 = wps->w.c [1].median [0], s1 = wps->w.c [1].median [1], s2 = wps->w.c [1].median [2];
    uint32_t holding_one = wps->w.holding_one, zeros_acc = wps->w.zeros_acc;
    int holding_zero = wps->w.holding_zero;

#define REFILL_TO(limit) do { \
        while (bc < (limit)) { \
            if (++ptr == end) { \
                wps->wvbits.error = 1; \
                ptr = wps->wvbits.buf; \
                wrap_error = 1; \
            } \
            sr |= (uint64_t) *ptr << bc; \
            bc += sizeof (*ptr) * 8; \
        } \
    } while (0)

#define GETBIT(result) do { \
        if (bc) { \
            result = (uint32_t) sr & 1; \
            sr >>= 1; \
            --bc; \
        } \
        else { \
            if (++ptr == end) { \
                wps->wvbits.error = 1; \
                ptr = wps->wvbits.buf; \
                wrap_error = 1; \
            } \
            bc = sizeof (*ptr) * 8 - 1; \
            result = *ptr & 1; \
            sr = (uint64_t) (*ptr >> 1); \
        } \
    } while (0)

#define READ_SIGNED(maxcode, base, result) do { \
        uint32_t rs_max = (maxcode); \
        uint32_t rs_sign; \
        if (rs_max < 2) { \
            if (rs_max) { \
                uint32_t rs_value; \
                REFILL_TO (2); \
                rs_value = (uint32_t) sr & 1; \
                rs_sign = (uint32_t) (sr >> 1) & 1; \
                sr >>= 2; \
                bc -= 2; \
                result = rs_sign ? ~(rs_value + (base)) : rs_value + (base); \
            } \
            else { \
                REFILL_TO (1); \
                rs_sign = (uint32_t) sr & 1; \
                sr >>= 1; \
                --bc; \
                result = rs_sign ? ~(base) : (base); \
            } \
        } \
        else { \
            int rs_bitcount = count_bits (rs_max); \
            const uint32_t rs_extras = bitset [rs_bitcount] - rs_max - 1; \
            uint32_t rs_code, rs_value; \
            REFILL_TO (rs_bitcount + 1); \
            rs_code = (uint32_t) sr & bitmask [rs_bitcount - 1]; \
            { \
                uint32_t rs_take = rs_code >= rs_extras; \
                rs_value = rs_take ? ((rs_code << 1) - rs_extras + \
                    (uint32_t) ((sr >> (rs_bitcount - 1)) & 1)) : rs_code; \
                rs_bitcount -= !rs_take; \
            } \
            rs_sign = (uint32_t) (sr >> rs_bitcount) & 1; \
            sr >>= rs_bitcount + 1; \
            bc -= rs_bitcount + 1; \
            rs_value += (base); \
            result = rs_sign ? ~rs_value : rs_value; \
        } \
    } while (0)

#define FAST_GET_MED(med) (((med) >> 4) + 1)
#define INC_MED(med, div) ((med) += (((med) + (div)) / (div)) * 5)
#define DEC_MED(med, div) ((med) -= (((med) + ((div) - 2)) / (div)) * 2)

#define DECODE_SAMPLE(M0, M1, M2, S0, OUT) do { \
        int zero_sample = 0; \
        if (holding_zero) { \
            holding_zero = 0; \
            READ_SIGNED (FAST_GET_MED (M0) - 1, 0, OUT); \
            DEC_MED (M0, DIV0); \
        } \
        else { \
            if (M0 < 2 && !holding_one && S0 < 2) { \
                uint32_t mask; \
                int cbits; \
                if (zeros_acc) { \
                    if (--zeros_acc) { \
                        OUT = 0; \
                        zero_sample = 1; \
                    } \
                } \
                else { \
                    uint32_t bit; \
                    for (cbits = 0; cbits < 33; ++cbits) { \
                        GETBIT (bit); \
                        if (!bit) break; \
                    } \
                    if (cbits == 33) goto done; \
                    if (cbits < 2) \
                        zeros_acc = cbits; \
                    else { \
                        for (mask = 1, zeros_acc = 0; --cbits; mask <<= 1) { \
                            GETBIT (bit); \
                            if (bit) zeros_acc |= mask; \
                        } \
                        zeros_acc |= mask; \
                    } \
                    if (zeros_acc) { \
                        m0 = m1 = m2 = 0; \
                        s0 = s1 = s2 = 0; \
                        OUT = 0; \
                        zero_sample = 1; \
                    } \
                } \
            } \
            if (!zero_sample) { \
                uint32_t bit; \
                REFILL_TO (LIMIT_ONES); \
                ones_count = __builtin_ctz ((uint32_t) ~sr); \
                if (ones_count >= LIMIT_ONES) { \
                    bc -= LIMIT_ONES; \
                    sr >>= LIMIT_ONES; \
                    ones_count = LIMIT_ONES; \
                    while (ones_count < (LIMIT_ONES + 1)) { \
                        GETBIT (bit); \
                        if (!bit) break; \
                        ++ones_count; \
                    } \
                    if (ones_count == (LIMIT_ONES + 1)) goto done; \
                    if (ones_count == LIMIT_ONES) { \
                        uint32_t mask; \
                        int cbits; \
                        for (cbits = 0; cbits < 33; ++cbits) { \
                            GETBIT (bit); \
                            if (!bit) break; \
                        } \
                        if (cbits == 33) goto done; \
                        if (cbits < 2) \
                            ones_count = cbits; \
                        else { \
                            for (mask = 1, ones_count = 0; --cbits; mask <<= 1) { \
                                GETBIT (bit); \
                                if (bit) ones_count |= mask; \
                            } \
                            ones_count |= mask; \
                        } \
                        ones_count += LIMIT_ONES; \
                    } \
                } \
                else { \
                    bc -= ones_count + 1; \
                    sr >>= ones_count + 1; \
                } \
                low = holding_one; \
                holding_one = ones_count & 1; \
                holding_zero = ~ones_count & 1; \
                ones_count = (ones_count >> 1) + low; \
                if (ones_count == 0) { \
                    low = 0; \
                    high = FAST_GET_MED (M0) - 1; \
                    DEC_MED (M0, DIV0); \
                } \
                else { \
                    low = FAST_GET_MED (M0); \
                    INC_MED (M0, DIV0); \
                    if (ones_count == 1) { \
                        high = low + FAST_GET_MED (M1) - 1; \
                        DEC_MED (M1, DIV1); \
                    } \
                    else { \
                        low += FAST_GET_MED (M1); \
                        INC_MED (M1, DIV1); \
                        if (ones_count == 2) { \
                            high = low + FAST_GET_MED (M2) - 1; \
                            DEC_MED (M2, DIV2); \
                        } \
                        else { \
                            low += (ones_count - 2) * FAST_GET_MED (M2); \
                            high = low + FAST_GET_MED (M2) - 1; \
                            INC_MED (M2, DIV2); \
                        } \
                    } \
                } \
                READ_SIGNED (high - low, low, OUT); \
            } \
        } \
        ++csamples; \
    } while (0)

    if (mono) {
        while (csamples < nsamples)
            DECODE_SAMPLE (m0, m1, m2, s0, buffer [csamples]);
    }
    else {
        while (csamples + 1 < nsamples) {
            DECODE_SAMPLE (m0, m1, m2, s0, buffer [csamples]);
            DECODE_SAMPLE (s0, s1, s2, m0, buffer [csamples]);
        }

        if (csamples < nsamples)
            DECODE_SAMPLE (m0, m1, m2, s0, buffer [csamples]);
    }

done:
    if (wrap_error) {
        ptr = wps->wvbits.buf;
        bc = 0;
        sr = 0;
    }
    else {
        while (bc > 32) {
            --ptr;
            bc -= sizeof (*ptr) * 8;
        }

        sr &= ((uint64_t) 1 << bc) - 1;
    }

    wps->wvbits.ptr = ptr;
    wps->wvbits.bc = bc;
    wps->wvbits.sr = (uint32_t) sr;

    wps->w.c [0].median [0] = m0;
    wps->w.c [0].median [1] = m1;
    wps->w.c [0].median [2] = m2;
    wps->w.c [1].median [0] = s0;
    wps->w.c [1].median [1] = s1;
    wps->w.c [1].median [2] = s2;

    wps->w.holding_one = holding_one;
    wps->w.holding_zero = holding_zero;
    wps->w.zeros_acc = zeros_acc;

#undef REFILL_TO
#undef GETBIT
#undef READ_SIGNED
#undef FAST_GET_MED
#undef INC_MED
#undef DEC_MED
#undef DECODE_SAMPLE

    return mono ? csamples : (csamples / 2);
}

#endif
