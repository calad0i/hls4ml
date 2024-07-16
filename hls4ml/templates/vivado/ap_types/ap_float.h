#ifndef __AP_FLOAT_H__
#define __AP_FLOAT_H__

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

constexpr int ceillog2(int x) { return (x <= 2) ? 1 : 1 + ceillog2((x + 1) / 2); }

#include "ap_fixed.h"
#include <cassert>
#include <cmath>

template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
INLINE ap_fixed_base<_AP_W, _AP_I, false, _AP_Q, _AP_O, _AP_N>
abs(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> value) {
    if (value > 0) {
        return value;
    } else {
        return -value;
    }
}

template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
INLINE ap_uint<ceillog2(_AP_W)> msb_loc(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> value) {
    auto buf = abs(value);
    ap_uint<ceillog2(_AP_W)> msb = 0;
    for (int i = 0; i < _AP_W; i++) {
#pragma HLS UNROLL
        if (buf[i] && msb < i) {
            msb = i;
        }
    }
    return msb;
}

template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
auto rt_floorlog2(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> value)
    -> ap_int_base<ceillog2(MAX(_AP_I - _AP_S, _AP_W - _AP_I + 1)) + (_AP_W - _AP_I > 0), (_AP_W - _AP_I > 0)> {
    // Runtime floorlog2 for fixed point numbers
    auto msb = msb_loc(value);
    ap_int_base<6, 1> r = msb;
    return r - (_AP_W - _AP_I);
}

template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
auto rt_ceillog2(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> value) -> decltype(rt_floorlog2(value)) {
    // Runtime ceillog2 for fixed point numbers
    auto msb = msb_loc(value);
    ap_ufixed<_AP_W, _AP_I> buf;
    if (value > 0) {
        buf = value;
    } else {
        buf = -value;
    }
    bool is_pow2 = buf << (_AP_W - msb) == 0;
    return msb + !is_pow2 - (_AP_W - _AP_I);
}

template <int M, int E, int E0 = 0> struct ap_float {

    bool is_negative;
    typedef ap_ufixed<M, 0, AP_RND, AP_SAT> mantissa_t;
    typedef ap_fixed<E, E, AP_RND, AP_SAT> exponent_t;
    typedef ap_ufixed<M + 1, 1> opr_mantissa_t;
    mantissa_t mantissa;
    exponent_t exponent;
    INLINE ap_float() {}
    INLINE ap_float(mantissa_t mantissa, exponent_t exponent) : mantissa(mantissa), exponent(exponent) {}

    template <typename T> INLINE ap_float(T value) { *this = ap_fixed<32, 16>(value); }

    template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
    INLINE int operator=(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> value) {
        is_negative = value < 0;
        auto abs_value = abs(value);
        exponent = rt_floorlog2(abs_value) - E0; // E0 is the offset
        mantissa = (abs_value >> (ap_int<E>(exponent) + E0)) - 1;
        return 0;
    }

    INLINE ap_fixed<M + 2 + (1 << E) - 1, 1 + (1 << (E - 1)) - 1 + E0> to_ap_fixed() const {
        if (exponent == -(1 << (E - 1))) {
            return 0;
        }
        ap_ufixed<M + (1 << E), (1 << (E - 1))> _result = mantissa + 1;
        ap_ufixed<M + (1 << E), (1 << (E - 1)) + E0> result = 0;
        _result = _result << exponent;
        result.range() = _result.range();
        if (is_negative) {
            result = -result;
        }
        return result;
    }

    INLINE double to_float() const { return to_ap_fixed().to_float(); }
    INLINE operator float() const { return to_float(); }
    INLINE operator double() const { return to_float(); }

    template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
    auto operator*(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> other)
        -> ap_ufixed<M + (1 << E) - 1 + _AP_W, 1 + (1 << (E - 1)) - 1 + E0 + _AP_I> {
        if (exponent == -(1 << (E - 1))) {
            return 0;
        }
        opr_mantissa_t mantissa = this->mantissa + 1;
        auto shift = exponent + E0;
        constexpr int E_max = (1 << (E - 1)) - 1 + E0;
        constexpr int E_min = -(1 << (E - 1)) + E0;
        constexpr int I = 1 + E_max + _AP_I;
        constexpr int W = M + E_max - E_min + _AP_W;
        ap_fixed<W, I> result_fixed = mantissa * other;
        if (is_negative) {
            result_fixed = -result_fixed;
        }
        result_fixed <<= shift;
        return result_fixed;
    }

    template <int _M, int _E, int _E0> ap_float<_M, _E, _E0> operator*(ap_float<_M, _E, _E0> other) {
        assert(0) && "ap_float can only be multiplied by ap_fixed";
    }
};

template <int M, int E, int E0, int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
auto operator*(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> a,
               ap_float<M, E, E0> b) -> ap_fixed<M + (1 << E) - 1 + _AP_W, 1 + (1 << (E - 1)) - 1 + E0 + _AP_I> {
    return b * a;
}

#ifndef __SYNTHESIS__
template <int M, int E, int E0> int operator>>(std::istringstream s, ap_float<M, E, E0> &b) {
    std::string str;
    s >> str;
    b = std::stof(str);
    return 0;
}
#endif

#endif