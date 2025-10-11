#pragma once
#include "instance.h"

namespace ratio_detail 
{
	constexpr double EPS = 1e-9;   // decision epsilon
	constexpr double TAU = 0.01;   // ~1% baseline
	constexpr double K_NEG = 0.5;
	constexpr double CAP = 5.0;

	inline double posneg(double d_norm, double ds) {
		if (d_norm > EPS) {
			// guard denominator so /fp:fast can’t make it tiny
			double denom = d_norm + TAU;
			if (denom < TAU + EPS) denom = TAU + EPS;
			return ds / denom;
		}
		else if (d_norm < -EPS) {
			double benefit = std::min(CAP, (-d_norm) / (TAU + EPS));
			return ds * (1.0 + K_NEG * benefit);
		}
		else {
			// treat “almost zero” as neutral
			return ds;
		}
}

	inline double norm(double delta, double bound) {
		// prevent divide by ~0 if a bound is bogus
		return delta / std::max(bound, 1e-9);
	}
}

#if defined(__GNUC__) || defined(__clang__)
#define FORCE_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

FORCE_INLINE double ratio_scorediff(double /*dt*/, double ds, double /*dw*/, double /*dv*/,
    double /*t_max*/, double /*w_max*/, double /*v_max*/) noexcept {
    return ds;
}

FORCE_INLINE double ratio_scorediff_time(double dt, double ds, double /*dw*/, double /*dv*/,
    double t_max, double /*w_max*/, double /*v_max*/) noexcept {
    using namespace ratio_detail;
    return posneg(norm(dt, t_max), ds);
}

FORCE_INLINE double ratio_scorediff_weight(double /*dt*/, double ds, double dw, double /*dv*/,
    double /*t_max*/, double w_max, double /*v_max*/) noexcept {
    using namespace ratio_detail;
    return posneg(norm(dw, w_max), ds);
}

FORCE_INLINE double ratio_scorediff_volume(double /*dt*/, double ds, double /*dw*/, double dv,
    double /*t_max*/, double /*w_max*/, double v_max) noexcept {
    using namespace ratio_detail;
    return posneg(norm(dv, v_max), ds);
}

using RatioFn = double(*)(double, double, double, double, double, double, double);
