/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "bqvec/VectorOpsComplex.h"
#include "bqvec/VectorOps.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>
#include <iostream>

using namespace breakfastquay;

using namespace std;

BOOST_AUTO_TEST_SUITE(TestVectorOpsComplex)

#ifdef USE_APPROXIMATE_ATAN2
static const double eps = 5.0e-3;
#else
#ifdef USE_SINGLE_PRECISION_COMPLEX
static const double eps = 1.0e-7;
#else
static const double eps = 1.0e-14;
#endif
#endif
    
#define COMPARE_N(a, b, n) \
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i] - b[cmp_i], eps); \
    }

#define COMPARE_NC(a, b, n) \
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i] - b[cmp_i], (bq_complex_element_t) eps); \
    }

#define COMPARE_CPLX_N(a, b, n)						\
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i].re - b[cmp_i].re, (bq_complex_element_t) eps); \
        BOOST_CHECK_SMALL(a[cmp_i].im - b[cmp_i].im, (bq_complex_element_t) eps); \
    }

BOOST_AUTO_TEST_CASE(add)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_t b[] = { { -1.0, 3.0 }, { -4.5, 0.0 } };
    bq_complex_t expected[] = { { 0.0, 5.0 }, { -1.5, -4.0 } };
    v_add(a, b, 2);
    COMPARE_CPLX_N(a, expected, 2);
}

BOOST_AUTO_TEST_CASE(add_with_gain)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_t b[] = { { -1.0, 3.0 }, { -4.5, 0.0 } };
    bq_complex_t expected[] = { { -0.5, 6.5 }, { -3.75, -4.0 } };
    v_add_with_gain(a, b, (bq_complex_element_t) 1.5, 2);
    COMPARE_CPLX_N(a, expected, 2);
}

BOOST_AUTO_TEST_CASE(multiply)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_t b[] = { { -1.0, 3.0 }, { -4.5, 0.0 } };
    bq_complex_t expected[] = { { -7.0, 1.0 }, { -13.5, 18.0 } };
    v_multiply(a, b, 2);
    COMPARE_CPLX_N(a, expected, 2);
}

BOOST_AUTO_TEST_CASE(multiply_to)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_t b[] = { { -1.0, 3.0 }, { -4.5, 0.0 } };
    bq_complex_t o[2];
    bq_complex_t expected[] = { { -7.0, 1.0 }, { -13.5, 18.0 } };
    v_multiply_to(o, a, b, 2);
    COMPARE_CPLX_N(o, expected, 2);
}

BOOST_AUTO_TEST_CASE(cartesian_to_magnitudes_bq)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_element_t o[2];
    bq_complex_element_t expected[] = { sqrt(5.0), 5.0 };
    v_cartesian_to_magnitudes(o, a, 2);
    COMPARE_NC(o, expected, 2);
}

BOOST_AUTO_TEST_CASE(cartesian_to_magnitudes)
{
    double re[] = { 1.0, 3.0 };
    double im[] = { 2.0, -4.0 };
    double o[2];
    double expected[] = { sqrt(5.0), 5.0 };
    v_cartesian_to_magnitudes(o, re, im, 2);
    COMPARE_N(o, expected, 2);
}

BOOST_AUTO_TEST_CASE(cartesian_interleaved_to_magnitudes)
{
    double a[] = { 1.0, 2.0, 3.0, -4.0 };
    double o[2];
    double expected[] = { sqrt(5.0), 5.0 };
    v_cartesian_interleaved_to_magnitudes(o, a, 2);
    COMPARE_N(o, expected, 2);
}

BOOST_AUTO_TEST_CASE(cartesian_to_polar_bq)
{
    bq_complex_t a[] = { { 0.0, 0.0 }, { 1.0, 1.0 }, { 0.0, -1.0 } };
    bq_complex_element_t mo[3], po[3];
    bq_complex_element_t me[] = { 0.0, sqrt(2.0), 1.0 };
    bq_complex_element_t pe[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    v_cartesian_to_polar(mo, po, a, 3);
    COMPARE_NC(mo, me, 3);
    COMPARE_NC(po, pe, 3);
}

BOOST_AUTO_TEST_CASE(cartesian_to_polar_interleaved_bq)
{
    bq_complex_t a[] = { { 0.0, 0.0 }, { 1.0, 1.0 }, { 0.0, -1.0 } };
    bq_complex_element_t o[6];
    bq_complex_element_t e[] = { 0.0, 0.0, sqrt(2.0), M_PI / 4.0, 1.0, -M_PI * 0.5 };
    v_cartesian_to_polar_interleaved(o, a, 3);
    COMPARE_NC(o, e, 6);
}

BOOST_AUTO_TEST_CASE(cartesian_to_polar)
{
    double re[] = { 0.0, 1.0, 0.0 };
    double im[] = { 0.0, 1.0, -1.0 };
    double mo[3], po[3];
    double me[] = { 0.0, sqrt(2.0), 1.0 };
    double pe[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    v_cartesian_to_polar(mo, po, re, im, 3);
    COMPARE_N(mo, me, 3);
    COMPARE_N(po, pe, 3);
}

BOOST_AUTO_TEST_CASE(cartesian_to_polar_interleaved_inplace)
{
    double a[] = { 0.0, 0.0, 1.0, 1.0, 0.0, -1.0 };
    double e[] = { 0.0, 0.0, sqrt(2.0), M_PI / 4.0, 1.0, -M_PI * 0.5 };
    v_cartesian_to_polar_interleaved_inplace(a, 3);
    COMPARE_N(a, e, 6);
}

BOOST_AUTO_TEST_CASE(cartesian_interleaved_to_polar)
{
    double a[] = { 0.0, 0.0, 1.0, 1.0, 0.0, -1.0 };
    double mo[3], po[3];
    double me[] = { 0.0, sqrt(2.0), 1.0 };
    double pe[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    v_cartesian_interleaved_to_polar(mo, po, a, 3);
    COMPARE_N(mo, me, 3);
    COMPARE_N(po, pe, 3);
}

BOOST_AUTO_TEST_CASE(polar_to_cartesian_bq)
{
    bq_complex_element_t m[] = { 0.0, sqrt(2.0), 1.0 };
    bq_complex_element_t p[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    bq_complex_t o[3];
    bq_complex_t e[] = { { 0.0, 0.0 }, { 1.0, 1.0 }, { 0.0, -1.0 } };
    v_polar_to_cartesian(o, m, p, 3);
    COMPARE_CPLX_N(o, e, 3);
}

BOOST_AUTO_TEST_CASE(polar_to_cartesian_interleaved_bq)
{
    bq_complex_t a[] = { { 0.0, 0.0 }, { 1.0, 1.0 }, { 0.0, -1.0 } };
    bq_complex_element_t o[6];
    bq_complex_element_t e[] = { 0.0, 0.0, sqrt(2.0), M_PI / 4.0, 1.0, -M_PI * 0.5 };
    v_cartesian_to_polar_interleaved(o, a, 3);
    COMPARE_NC(o, e, 6);
}

BOOST_AUTO_TEST_CASE(polar_to_cartesian)
{
    double m[] = { 0.0, sqrt(2.0), 1.0 };
    double p[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    double ro[3], io[3];
    double re[] = { 0.0, 1.0, 0.0 };
    double ie[] = { 0.0, 1.0, -1.0 };
    v_polar_to_cartesian(ro, io, m, p, 3);
    COMPARE_N(ro, re, 3);
    COMPARE_N(io, ie, 3);
}

BOOST_AUTO_TEST_CASE(polar_to_cartesian_interleaved_inplace)
{
    double a[] = { 0.0, 0.0, sqrt(2.0), M_PI / 4.0, 1.0, -M_PI * 0.5 };
    double e[] = { 0.0, 0.0, 1.0, 1.0, 0.0, -1.0 };
    v_polar_interleaved_to_cartesian_inplace(a, 3);
    COMPARE_N(a, e, 6);
}

BOOST_AUTO_TEST_CASE(polar_to_cartesian_interleaved)
{
    double m[] = { 0.0, sqrt(2.0), 1.0 };
    double p[] = { 0.0, M_PI / 4.0, -M_PI * 0.5 };
    double o[6];
    double e[] = { 0.0, 0.0, 1.0, 1.0, 0.0, -1.0 };
    v_polar_to_cartesian_interleaved(o, m, p, 3);
    COMPARE_N(o, e, 6);
}

BOOST_AUTO_TEST_SUITE_END()

