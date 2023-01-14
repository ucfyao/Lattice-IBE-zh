#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "FFT.h"
#include "params.h"

using namespace std;
using namespace NTL;

// 门限密码学中的离散傅里叶变换
//  f_fft 为 complex<RR_t> CC_t 复数类型。 f 为 long double 类型。N 为无符号长整数类型。w0 为 复数类型。
void FFTStep(CC_t * const f_fft, RR_t const * const f, const unsigned long N, const CC_t w0)
{
    // 根据N是否为1进行处理
    if(N==1)
    {
        // 如果为1，进行赋值
        f_fft[0] = f[0];
    }
    else
    {
        if(N==2)
        {
            // 如果 N 为2，数组进行前2位赋值
            f_fft[0] = f[0] + ii*f[1];
            f_fft[1] = f[0] - ii*f[1];
        }
        else
        {
            // 如果 n%2的结果为0 进行错误提示。
            assert(N%2==0);
            // 初始化2个 RR_t类型数组。此处为2分算法。
            RR_t f0[N0/2], f1[N0/2];
            // 初始化
            CC_t f0_fft[N0/2], f1_fft[N0/2], wk, w02;
            // 初始化无符号整数 k。
            unsigned int k;
            // 进行数组拆分
            for(k=0; k<(N/2); k++)
            {
                f0[k] = f[2*k];
                f1[k] = f[2*k+1];
            }

            w02 = w0*w0;
            wk = w0;
            // 递归调用，进行迭代
            FFTStep(f0_fft, f0, N/2, w02);
            FFTStep(f1_fft, f1, N/2, w02);

            // 将计算结果合并到数组 f_fft
            for(k=0; k<N; k++)
            {
                f_fft[k] = f0_fft[k%(N/2)] + wk*f1_fft[k%(N/2)];
                wk *= w02;
            }
        }
    }
}

// 逆FFT步骤
void ReverseFFTStep(CC_t * const f, CC_t const * const f_fft, const unsigned long N, const CC_t w0)
{
    if(N!=2)
    {
        // 断言，如果 N%2==0 进行报错提示。
        assert(N%2==0);
        // 初始化 CC_t类型数组 f0[N0/2], f1[N0/2];
        CC_t f0[N0/2], f1[N0/2];
        // 初始化CC_t类型数组 f0_fft[N0/2], f1_fft[N0/2],变量  w02, wk
        CC_t f0_fft[N0/2], f1_fft[N0/2], w02, wk;
        // 初始化无符号整型 k
        unsigned int k;

        w02 = w0*w0;
        wk = w0;

        // 进行翻转计算
        for(k=0; k<N/2; k++)
        {
            f0_fft[k] = (f_fft[k] + f_fft[k+(N/2)])*0.5l;
            f1_fft[k] = wk*(f_fft[k] - f_fft[k+(N/2)])*0.5l;
            wk *= w02;
        }
        // 2分计算
        ReverseFFTStep(f0, f0_fft, (N/2), w02);
        ReverseFFTStep(f1, f1_fft, (N/2), w02);

        // 合并数组
        for(k=0; k<N/2; k++)
        {
            f[2*k] = f0[k];
            f[2*k+1] = f1[k];
        }
    }
    else
    {
        // 赋值
        f[0] = (f_fft[0] + f_fft[1])*0.5l;
        f[1] = (f_fft[0] - f_fft[1])*(-0.5l*ii);
    }
}

// 复数FFT的逆变换，接受 double 类型数据
void MyRealReverseFFT(double * const f, CC_t const * const f_fft)
{
    // 初始化N0长度的复数类数组 fprime
    CC_t fprime[N0];
    // 初始化无符号整型 i
    unsigned int i;
    // 进行逆计算
    ReverseFFTStep(fprime, f_fft, N0, omega_1);

    // 循环赋值
    for(i=0; i<N0; i++)
    {
        // real返回实数部分
        f[i] = (fprime[i]).real();
    }
}

// 复数FFT的逆变换，接受 long int 类型数据
void MyIntReverseFFT(long int * const f, CC_t const * const f_fft)
{
    CC_t fprime[N0];
    unsigned int i;

    ReverseFFTStep(fprime, f_fft, N0, omega_1);
    for(i=0; i<N0; i++)
    {
        f[i] = ((long int) round( fprime[i].real() ) );
    }
}

// 将int 类型的进行类型转换后再进行FFT
void MyIntFFT(CC_t * f_FFT, const long int * const f)
{
    // 初始化RR_t类型变量
    RR_t f_double[N0];
    // 初始化int变量
    unsigned int i;

    // 进行类型转换
    for(i=0; i<N0; i++)
    {
        f_double[i] = ( RR_t ( f[i] ) );
    }

    // 进行FFT
    FFTStep(f_FFT, f_double, N0, omega);
}

// 将ZZX转为FFT。
// FFT 快速傅里叶变换(Fast Fourier Transform, FFT) 在多项式相乘的过程中可以减少运算，提高程序的效率。
// FFT就是能够在O(nlogb)时间内将系数表示法转化为点值表示法，相乘后再将点值表示法转为系数表示法，然后实现卷积。
void ZZXToFFT(CC_t * f_FFT, const ZZX f)
{
    // 初始化  long double 类型的 N0 长度数组 f_double
    RR_t f_double[N0];
    // 初始化无符号整型 i
    unsigned int i;

    // 断言。如果 多项式f的 degree 为 N0-1，进行错误提示
    assert(deg(f)==N0-1);
    // 断言，如果f系数的最大NumBits小于900，进行错误提示
    assert(MaxBits(f)<900);
    // 循环进行数据类型转换为RR_t类型。
    for(i=0; i<N0; i++)
    {
        f_double[i] = ( RR_t ( conv<double>(f[i]) ) );
    }

    FFTStep(f_FFT, f_double, N0, omega);
}


void FFTToZZX(ZZX& f, CC_t const * const f_FFT)
{
    double f_double[N0];
    unsigned int i;
    MyRealReverseFFT(f_double, f_FFT);

    f.SetLength(N0);
    for(i=0; i<N0; i++)
    {
        f[i] = conv<ZZ>(round(f_double[i]));
    }

}
