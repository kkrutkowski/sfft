#ifndef MINTRIG_HPP
#define MINTRIG_HPP

#include <immintrin.h>

//Most ninja tricks used here:
//http://fastcpp.blogspot.fr/2011/03/changing-sign-of-float-values-using-sse.html
//http://www.songho.ca/misc/sse/sse.html
//http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/
//http://www.masmforum.com/board/index.php?PHPSESSID=786dd40408172108b65a5a36b09c88c0&topic=9515.0
//http://cbloomrants.blogspot.fr/2010/11/11-20-10-function-approximation-by_20.html
//http://assemblyrequired.crashworks.org/2009/10/16/timing-square-root/
//http://nghiaho.com/?p=997
//http://www.researchgate.net/publication/3321724_Efficient_approximations_for_the_arctangent_function
//http://www.ganssle.com/approx/approx.pdf
//http://forum.allaboutcircuits.com/newsgroups/viewtopic.php?t=68185

const float invtwopi=0.1591549f;
const float twopi=6.283185f;
const float threehalfpi=4.7123889f;
const float pi=3.141593f;
const float halfpi=1.570796f;
const float quarterpi=0.7853982f;

static const __m256 AVX_SIGNMASK =  _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

//SCALAR
namespace FTA{
    __m256 sqrt_ps(__m256 squared);
    __m256 cos_52s_ps(__m256 x);
    __m256 cos_ps(__m256 angle);
    __m256 sin_ps(__m256 angle);
    void sincos_ps(__m256, __m256*, __m256*);
}



__m256 FTA::sqrt_ps(__m256 squared){
    return _mm256_sqrt_ps(squared);
}

// FMA implementation
__m256 FTA::cos_52s_ps(__m256 x) {
    const __m256 c1 = _mm256_set1_ps(0.9999932946f);
    const __m256 c2 = _mm256_set1_ps(-0.4999124376f);
    const __m256 c3 = _mm256_set1_ps(0.0414877472f);
    const __m256 c4 = _mm256_set1_ps(-0.0012712095f);
    __m256 x2 = _mm256_mul_ps(x, x);

    // Using FMA instructions for more efficient computation
    __m256 result = _mm256_fmadd_ps(c4, x2, c3); // c3 + c4 * x2
    result = _mm256_fmadd_ps(result, x2, c2);    // c2 + (c3 + c4 * x2) * x2
    result = _mm256_fmadd_ps(result, x2, c1);    // c1 + (c2 + (c3 + c4 * x2) * x2) * x2

    return result;
}

__m256 FTA::cos_ps(__m256 angle){
    //clamp to the range 0..2pi

    //take absolute value
    angle=_mm256_andnot_ps(AVX_SIGNMASK,angle);
    //fmod(angle,twopi)
    angle=_mm256_sub_ps(angle,_mm256_mul_ps(_mm256_floor_ps(_mm256_mul_ps(angle,_mm256_set1_ps(invtwopi))),_mm256_set1_ps(twopi)));
    //angle = _mm256_fnmadd_ps(_mm256_floor_ps(_mm256_mul_ps(angle, _mm256_set1_ps(invtwopi))), _mm256_set1_ps(twopi), angle); //seemingly slower than AVX

    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(pi),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ), AVX_SIGNMASK));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(twopi),angle))));

    __m256 result=FTA::cos_52s_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_LT_OQ)), AVX_SIGNMASK));
    return result;
}

__m256 FTA::sin_ps(__m256 angle){
    return FTA::cos_ps(_mm256_sub_ps(_mm256_set1_ps(halfpi),angle));
}


void FTA::sincos_ps(__m256 angle, __m256 *sin, __m256 *cos){
    __m256 anglesign=_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(AVX_SIGNMASK,angle));
    //clamp to the range 0..2pi

    //take absolute value
    angle=_mm256_andnot_ps(AVX_SIGNMASK,angle);
    //fmod(angle,twopi)
    angle=_mm256_sub_ps(angle,_mm256_mul_ps(_mm256_floor_ps(_mm256_mul_ps(angle,_mm256_set1_ps(invtwopi))),_mm256_set1_ps(twopi)));
    // angle = _mm256_fnmadd_ps(_mm256_floor_ps(_mm256_mul_ps(angle, _mm256_set1_ps(invtwopi))), _mm256_set1_ps(twopi), angle); //seemingly slower than AVX


    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(pi),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ),AVX_SIGNMASK));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(twopi),angle))));

    __m256 result=FTA::cos_52s_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_LT_OQ)), AVX_SIGNMASK));
    *cos=result;

    __m256 sinmultiplier=_mm256_mul_ps(anglesign,_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ),AVX_SIGNMASK)));
    *sin=_mm256_mul_ps(sinmultiplier,FTA::sqrt_ps(_mm256_fnmadd_ps(result, result, _mm256_set1_ps(1.f))));

    return;
}


#endif // MINTRIG_HPP
