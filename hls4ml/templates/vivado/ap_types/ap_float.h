#include "ap_fixed.h"
#include <cassert>
#include <cmath>

template <int M, int E, int E0 = 0> struct ap_float {
    typedef ap_fixed<M, 1> mantissa_t;
    typedef ap_int<E> exponent_t;
    mantissa_t mantissa;
    exponent_t exponent;

    INLINE ap_float() {}
    INLINE ap_float(mantissa_t mantissa, exponent_t exponent) : mantissa(mantissa), exponent(exponent) {}
    INLINE ap_float(float value) { *this = value; }
    INLINE int operator=(float value) {
        int exponent = 0;
        float mantissa = frexpf(value, &exponent);
        if (mantissa == -0.5) {
            mantissa *= 2;
            exponent -= 1;
        }
        this->mantissa = ap_fixed<M, 1, AP_RND, AP_SAT>(mantissa);
        this->exponent = exponent - E0;
        return 0;
    }

    INLINE ap_fixed<M + (1 << E) - 1, 1 + (1 << (E - 1)) - 1 + E0> to_ap_fixed() const {
        ap_fixed<M + (1 << E) - 1, 1 + (1 << (E - 1)) - 1> _result = mantissa;
        ap_fixed<M + (1 << E) - 1, 1 + (1 << (E - 1)) - 1 + E0> result = 0;
        _result = _result << exponent;
        result.range() = _result.range();
        return result;
    }

    INLINE double to_float() const { return to_ap_fixed().to_float(); }
    INLINE operator float() const { return to_float(); }
    INLINE operator double() const { return to_float(); }

    template <int _AP_W, int _AP_I, bool _AP_S, ap_q_mode _AP_Q, ap_o_mode _AP_O, int _AP_N>
    auto operator*(ap_fixed_base<_AP_W, _AP_I, _AP_S, _AP_Q, _AP_O, _AP_N> other)
        -> ap_fixed<M + (1 << E) - 1 + _AP_W, 1 + (1 << (E - 1)) - 1 + E0 + _AP_I> {
        auto shift = exponent + E0;
        constexpr int E_max = (1 << (E - 1)) - 1 + E0;
        constexpr int E_min = -(1 << (E - 1)) + E0;
        constexpr int I = 1 + E_max + _AP_I;
        constexpr int W = M + E_max - E_min + _AP_W;
        ap_fixed<W, I> result_fixed = mantissa * other;
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