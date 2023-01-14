/*
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-11 10:26:14
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-12 13:03:54
 * @FilePath     : /IBBE/Random.cc
 * @Description  :
 * Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "Random.h"
#include "params.h"

using namespace std;
using namespace NTL;

// 生成随机向量
vec_ZZ RandomVector()
{
    // 初始化向量
    vec_ZZ w;
    // 初始化整型变量
    unsigned int i;
    // 设置长度
    w.SetLength(N0);

    // 随机赋值
    for(i=0; i<N0; i++)
    {
        w[i] = conv<ZZ>(rand()) % q1;
    }
    return w;
}

//==============================================================================
//Generates a random polynomial of fixed degree
//生成固定次数的随机多项式
//==============================================================================
ZZX RandomPoly(const unsigned int degree)
{
    unsigned int i;
    ZZX f;
    f.SetLength(degree+1);
    for(i=0; i<=degree; i++)
    {
        // 随机数
        f[i] = rand();
    }
    return f;
}

//==============================================================================
//Generates a random polynomial of fixed degree and "approximately" fixed squared norm
//==============================================================================
ZZX RandomPolyFixedSqNorm(const ZZ& SqNorm, const unsigned int degree)
{
    // 无符号整数
    unsigned int i;
    // 初始化大整数
    ZZ SqNorm0, Ratio;
    // 初始化多项式
    ZZX f;
    // 设置多项式的长度
    // all vector elements are initialized to zero
    f.SetLength(degree+1);

    // sqrt 求（SqNorm/degree+1）根
    RR_t sigma = sqrt( ( (double) conv<double>(SqNorm)/(degree+1) ) );

    // 获取随机数，赋值到多项式中。
    for(i=0; i<=degree; i++)
    {
        f[i] = conv<ZZ>(Sample3(sigma));
    }
    // 按位或并赋值。
    f[degree] |= 1;
    return f;
}
