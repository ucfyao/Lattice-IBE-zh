/***
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-09 12:39:55
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-09 12:40:09
 * @FilePath     : /IBBE/params.h
 * @Description  :
 * @Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */
#ifndef IBBE_PARAMS_H
#define IBBE_PARAMS_H

#include <math.h>
#include <complex.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>

using namespace NTL;
using std::complex;

//=====================================================================================
// These are the parameters you need to change
// N0 is the degree of the polynomial ring used. N0 must be a power of 2!
// q0 is the modulus w.r.t. whom the integers are reduced. We suggest to take q0 prime
// polynomial ring 多项式环，为抽象代数招知识
//=====================================================================================
#define N0 512
// 1<<30 表示1左移30位,每左移一位乘以2,所以就是1*2^30=1073741824。表示将一个运算对象的各二进制位全部左移若干位(左边的二进制位丢弃,右边补0)。
#define q0 (1<<30)

// 定义常量q1，转换类型。conv<T>,  T is the return type。
const ZZ q1 = conv<ZZ>(q0);

// 在不同的平台进行兼容
//#ifdef  USE_FLOAT128
//    typedef __float128 RR_t;
//#else
    typedef long double RR_t;
    // 复数类。
    typedef complex<RR_t> CC_t;
//#endif

// Data structure of the Master Secret Key
typedef struct
{
    ZZX PrK[4];
    CC_t PrK_fft[4][N0];
    RR_t GS_Norms[2*N0];
    RR_t sigma;
    RR_t B[2*N0][2*N0];
    RR_t Bstar[2*N0][2*N0];
} MSK_Data;

// Data Structure of Master Public Key
typedef struct
{
    ZZ_pX h;
    CC_t h_FFT[N0];
} MPK_Data;


//==============================================================================
//          Seed for the RNG
// 不同架构下的获取时钟指令有一些差异。像 i386 以及 arm 平台，所以，我们需要变更 rdtsc_time。
//==============================================================================

// intel 80386 架构
#ifdef __i386
extern __inline__ uint64_t rdtsc(void) {
  uint64_t x;
  __asm__ volatile ("rdtsc" : "=A" (x));
  return x;
}
//  amd x86_64架构
#elif defined __amd64
extern __inline__ uint64_t rdtsc(void) {
// 在AMD64 架构下的 rdtsc 指令的时间统计。
// 从 CPU 寄存器中直接获取 walllock 的汇编指令 RDTSC
// 拿到的 walllock 是 自1970.1.1 到现在 的CPU周期数，对于 GHZ 的CPU 来说单位就是ns。
  uint64_t a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d<<32) | a;
}
#elif defined(__aarch64__) // arm 64为架构
    {
        uint64_t t;

        __asm__ volatile("mrs %0,  cntvct_el0" : "=r"(t));
        return (t);
    }
#else
    return (0);
#endif

// Useful constants up to ~500 bits
const long double sigma_1= 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023L;//sqrt(1/(2*log(2)))
const long double log_2 = 0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516L;
const RR_t Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081L;
const RR_t PiPrime = 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263518268551099195455583724299621273062L; //1/sqrt(2*Pi)
const RR_t LDRMX = ((RR_t)RAND_MAX );
const CC_t ii(0, 1);
const CC_t omega   = exp( ii*(Pi/N0));
const CC_t omega_1 = exp(-ii*(Pi/N0));

#endif
