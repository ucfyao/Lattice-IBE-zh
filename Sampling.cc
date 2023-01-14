/*
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-11 10:32:28
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-11 14:41:03
 * @FilePath     : /IBBE/Sampling.cc
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

#include "Sampling.h"
#include "params.h"

using namespace std;
using namespace NTL;


//==============================================================================
// Takes in input a random value and samples from distribution D_{\sigma_2}^+,
// Samples an element in Z+ with probability proportionnal to 2^{-x^2}
// 对Z+中的元素采样，其概率比例为2^{-x^2}
//==============================================================================
unsigned int Sample0(unsigned long alea)
{
    // alea 按位与1UL；
    // 1UL---无符号长整型1。 如果没有UL后缀,则系统默认为 int类型,即,有符号整形。UL：unsigned long
    // 0&0=0;   0&1=0;    1&0=0;     1&1=1;
    // 即：两位同时为“1”，结果才为“1”，否则为0
    // 例如：3&5  即 0000 0011 & 0000 0101 = 0000 0001 因此，3&5的值得1。
    if((alea & 1UL) == 0UL)
    {
        return 0;
    }
    // 初始化无符号整型i
    unsigned int i;
    // 初始化无符号整型k
    unsigned int k = 1;
    // 初始化无符号长整型 mask
    unsigned long mask = 0;
    // 初始化无符号长整型 aux
    unsigned long aux;

    // 最多递归循环1000次
    for( i=1; i<1000;)
    {
        // 按二进制位进行 “与” 运算。两位同时为“1”，结果才为“1”，否则为0。
        // 特殊用途，进行清零。如果想将一个单元清零，即使其全部二进制位为0，只要与一个各位都为零的数值相与，结果为零。
        aux = (alea & mask);
        // 二进制位全部右移k位
        alea = (alea>>k);

        // 如果alea和mask相等，进行重新计算
        if(aux)
        {
            return Sample0(alea);
        }
        else
        {
            // 使用与运算。两位同时为“1”，结果才为“1”，否则为0。
            // 没有任何一个二进制位为 1，也就是 == 0 的话，返回 i
            if((alea & 1UL) == 0UL)
            {
                return i;
            }
        }
        // i进行 +1
        i++;
        // k 进行 +2
        k += 2;
        // | 位或运算。只要有一个为1，其值为1。
        // 0|0=0；   0|1=1；   1|0=1；    1|1=1；
        // 例如:3|5　即 0000 0011 | 0000 0101 = 0000 0111   因此，3|5的值得7。　
        mask = (mask << 2) | 6UL;
    }
    cout << "ERROR" << endl;
    return 999999;
}

//==============================================================================
// Samples from distribution D_{k\sigma_2}^+, ie
// Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2}
// 对Z+中的元素采样，其概率比例为2^{-（x/k）^2}
//==============================================================================
unsigned int Sample1(const unsigned int k)
{
    // 初始化无符号整型变量 x, y, z
    unsigned int x, y, z;
    // 初始化无符号长整型随机数 alea
    unsigned long alea = rand();

    // 使用 Sample0 计算随机采样
    x = Sample0(alea);
    // 产生 0~k 这k个整数中的一个随机整数。
    y = rand() % k;
    // 计算z
    z = k*x + y;
    // 计算w
    RR_t w = y*( (z<<1) - y );
    // 初始化变量 borne
    // exp(n)表示e的n次方
    RR_t borne =  LDRMX / exp( w*log_2/(k*k) );

    // 初始化无符号长整型随机数
    alea = rand();

    // 如果 alea > borne 重新进行计算，不然返回z。
    if(alea>borne)
    {
        return Sample1(k);
    }
    else
    {
        return z;
    }
    cout << "ERROR" << endl;
    return 999999;
}

//==============================================================================
// Samples from distribution D_{k\sigma_2}, ie
// Samples an element in Z with probability proportionnal to 2^{-(x/k)^2}
//==============================================================================
signed int Sample2(const unsigned int k)
{
    // 初始化有符号整型 signe
    signed int signe;
    // 初始化有符号整型 x
    signed int x;
    // 初始化无符号long类型 alea
    // rand() 会返回一随机数值,范围在 0 至 RAND_MAX 间
    unsigned long alea = rand();

    // 死循环
    while(1)
    {
        // 调用 Sample1 函数，生成随机数x
        x = Sample1(k);
        // 如果 x不为0 或者 alea 位与1后为1 进行计算。
        if( (x!=0) || ((alea & 1)==1) )
        {
            alea >>= 1;
            signe = 1 - 2*(alea & 1);
            x *= signe;
            return x;
        }
        // 将 alea  >> 右移1 位，假设k=2，换算成二进制是0010，右移一位变成了0001。
        // 低位舍去，高位补零。
        // 等同于 x = x/2。
        alea >>= 1;
    }
}

//==============================================================================
// Samples from distribution D_{sigma}, ie
// Samples an element in Z with probability proportionnal to e^{-x^2/2*(sigma^2)}
//==============================================================================
signed int Sample3(const RR_t sigma128)
{
    // 初始化有符号整型 x
    signed int x;
    // 初始化 double类型变量  alea, borne
    double alea, borne;

    // 初始化 long double（RR_t） 类型常量 sigma
    const RR_t sigma = sigma128;

    // 计算无符号 long 类型的常量 k
    // ceil 向上取整
    const unsigned long k = ( (unsigned long) ceil( (RR_t) sigma/sigma_1 ) );

    // 死循环
    while(1)
    {
        // 通过 Sample2 函数计算 x
        x = Sample2(k);
        // 获取 0 - RAND_MAX之间的随机数，初始化为 long double RR_t 类型。
        alea = ((RR_t)rand()) / LDRMX;
        // 赋值double类型的变量 borne。计算以自然常数e为底的指数。即:exp(n)表示e的n次方。
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) ) );
        // 断言。如果 borne<=1，代码会终止运行，并且会把源文件、错误的代码、以及行号都输出来。
        assert(borne <= 1);
        // 如果 alea<borne 返回 x
        if(alea < borne)
        {
            return x;
        }
    }
}


//==============================================================================
// Samples from distribution D_{c,sigma}, ie
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}
//以与e^{-（c-x）^2/2*（sigma^2）}成比例的概率采样Z中的元素
//==============================================================================
signed int Sample4(RR_t c, RR_t sigma)
{
    // 初始化 RR_t 类型变量 alea, borne;
    RR_t alea, borne;
    // 初始化有符号整型 x
    signed int x;
    // 初始化无符号整型 coin
    unsigned int coin;
    // 初始化有符号整型常量 intc
    const signed int intc = ( (signed int) floor(c) );
    // 初始化RR_t常量 fracc
    const RR_t fracc = c - intc;
    // 获取一个随机数
    coin = rand();

    // 初始化RR_t类型的常量 denom
    const RR_t denom = 1/(2*sigma*sigma);

    while(1)
    {
        // 调用 Sample3 函数获取 x
        x = Sample3(sigma);
        // x 增加 coin & 1。
        //  coin & 1。两位同时为“1”，结果才为“1”，否则为0
        x += (coin&1);
        // abs 求整型数据的绝对值。
        // 如果 |x| > 8, 输出 x
        if(abs(x)>8){cout << x << endl;}
        // coin 右移一位
        coin >>= 1;
        // 计算 borne
        borne = exp(-(x-fracc)*(x-fracc)*denom)/ ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );
        // 断言。如果 borne<1，代码会终止运行，并且会把源文件、错误的代码、以及行号都输出来。
        assert(borne<1);
        // 计算随机数 alea
        alea = ( (RR_t)rand() ) / LDRMX;
        // 如果 alea < borne，返回 x + intc
        if(alea < borne)
        {
            return (x + intc);
        }
    }
}