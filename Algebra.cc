/*
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-09 12:39:48
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-14 20:47:37
 * @FilePath     : /IBBE/Algebra.cc
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

#include "Algebra.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"

using namespace std;
using namespace NTL;

// 二次封装求模函数
// void trunc(ZZX& x, const ZZX& a, long m);
//  x = a % X^m
ZZX FastMod(const ZZX& f)
{
    //! - (f >> N0) 啥作用，强制转换类型吗？
    return (trunc(f, N0) - (f >> N0));
}

ZZX Cyclo()
{
    // ZZX 初始化。
    // 结构实例：[2 10 14 6]
    ZZX phi0;
    // all vector elements are initialized to zero
    phi0.SetLength(N0+1);
    phi0[0] = 1;
    phi0[N0] = 1;
    return phi0;
}

// 初始化一个N0长度的常量多项式 phi
const ZZX phi = Cyclo();

//==============================================================================
//Computes the squared norm of a polynomial f
// 计算多项式f的平方范数
//==============================================================================
ZZ SquaredNorm(const ZZX& f, const unsigned int degree)
{
    // 初始化无符号整型 i
    unsigned int i;
    // 初始化大整数 somme
    ZZ somme;

    // 循环计算 开方总和
    for(i=0; i<=degree; i++)
    {
        somme += sqr(f[i]);
    }
    return somme;
}


//==============================================================================
//Computes f(1/x) mod (x^N + 1)
//If f = a0 + a1*x + ... + a_{N-1}*x^{N-1}, then
//Reverse(f) = a0 + a_{N-1}*x + ... + a1*x^{N-1}
// 反转计算
//==============================================================================
ZZX Reverse(const ZZX& f)
{
    // 断言，多项式 f 的次数 >= 0。次数一般为最大n次项
    assert(deg(f) >= 0);
    // 断言，多项式 f 的次数 < N0
    assert(deg(f) < N0);

    // 初始化一个多项式 fb
    ZZX fb;
    // 初始化无符号整型 i
    unsigned int i;
    // 初始化长度
    fb.SetLength(N0);
    // 赋值
    fb[0] = f[0];
    // ! 为啥2次赋值
    fb.SetLength(N0);

    // 从 N0-f 多项式的次数开始进行计算
    for(i = N0 - deg(f); i < N0; i++)
    {
        fb[i] = -f[N0-i];
    }
    fb[0] = f[0];
    return fb;
}

//==============================================================================
//Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
//计算多项式k，使得（F，G）<-（F，B）-k*（F，G）最小化（F，C）的大小
//==============================================================================
ZZX ReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, unsigned int & mb)
{
    // 初始化无符号整型 i
    unsigned int i;
    // 初始化大整数 a
    ZZ a;
    // 初始化多项式变量 fb, gb, num, den, iden, iphi, k
    ZZX fb, gb, num, den, iden, iphi, k;

    // 进行反转
    fb = Reverse(f);
    gb = Reverse(g);

    // 计算 trunc(fb*F + gb*G, N0)
    num = FastMod(fb*F + gb*G);
    // 计算 trunc(f*fb + g*gb, N0)
    den = FastMod(f*fb + g*gb);
    //获取最大的数
    mb = MaxBits(num);

    XGCD(a, iden, iphi, den, phi);

    // 计算 trunc(num*iden, N0)
    k = FastMod(num*iden);
    // 设置长度
    k.SetLength(N0);
    // 循环进行除运算
    for(i=0; i<N0; i++)
    {
        k[i] /= a;
    }

    return k;
}


//==============================================================================
//Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
//==============================================================================
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G)
{
    unsigned int i;
    ZZX k;
    CC_t f_FFT[N0], g_FFT[N0], F_FFT[N0], G_FFT[N0], num_FFT[N0], den_FFT[N0], k_FFT[N0];

    assert(MaxBits(f)<900);
    ZZXToFFT(f_FFT, f);

    assert(MaxBits(g)<900);
    ZZXToFFT(g_FFT, g);

    assert(MaxBits(F)<900);
    ZZXToFFT(F_FFT, F);

    assert(MaxBits(G)<900);
    ZZXToFFT(G_FFT, G);

    for(i=0; i<N0; i++)
    {
        num_FFT[i] = f_FFT[N0-1-i]*F_FFT[i] + g_FFT[N0-1-i]*G_FFT[i];
        den_FFT[i] = f_FFT[N0-1-i]*f_FFT[i] + g_FFT[N0-1-i]*g_FFT[i];
        k_FFT[i] = num_FFT[i]/den_FFT[i];
    }

    FFTToZZX(k, k_FFT);
    return k;
}

//==============================================================================
//Returns the anticircular matrix associated to integral polynomial f and integer N
//返回与整数多项式f和整数N相关的反圆形矩阵
//==============================================================================
mat_ZZ AnticircularMatrix(const ZZX& f)
{
    // 初始化无符号整型
    unsigned int i,j;
    // 初始化整型 df
    int df;
    // 初始化大整数矩阵
    mat_ZZ M;
    // 设置大小
    M.SetDims(N0, N0);
    // 获取多项式次数
    df = deg(f);
    // 多项式为空，直接返回空矩阵
    if(df==-1)
    {
        return M;
    }
    // 初始化无符号整型dfu
    unsigned dfu;
    // 赋值
    dfu = ((unsigned) df);
    // 如果 dfu >=N0 输出
    if(dfu>=N0)
    {
        cout << "df = " << dfu << endl;
        cout << "f = " << f << endl;
    }

    // 断言，dfu<N0进行报错提示
    assert(dfu<N0);

    //循环赋值
    for(i=0; i<N0; i++)
    {
        for(j=i; ((j<=dfu+i) && (j<N0)); j++)
        {
            M[i][j] = f[j-i];
        }
        for(j=0; (j+N0)<=(dfu+i); j++)
        {
            M[i][j] = -f[j-i+N0];
        }
    }
    return M;
}


//==============================================================================
//Generates a basis from the double pair (f,g), (F,G) and N
// 生成基
//This basis has the form :
//    |f g|
//M = |F G|
//==============================================================================
mat_ZZ BasisFromPolynomials(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G)
{
    // 初始化无符号整型  i,j;
    unsigned int i,j;
    // 初始化大整数矩阵 A, M
    mat_ZZ A,M;
    // 设置矩阵大小
    M.SetDims(2*N0, 2*N0);

    // 反圆形矩阵，生成f区域
    A = AnticircularMatrix(f);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i][j] = A[i][j];
    }}

    // 反圆形矩阵，生成g区域
    A = AnticircularMatrix(g);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i][j+N0] = A[i][j];
    }}

    // 反圆形矩阵，生成F区域
    A = AnticircularMatrix(F);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i+N0][j] = A[i][j];
    }}

    // 反圆形矩阵，生成G区域
    A = AnticircularMatrix(G);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i+N0][j+N0] = A[i][j];
    }}
    // 返回矩阵
    return M;
}

//==============================================================================
//Computes the Inverse of f (mod phi) (mod q)
// 计算倒数
//==============================================================================
ZZ_pX Inverse(const ZZX& f)
{
    // 初始化 integers mod p
    ZZ_p::init(q1);
    // 初始化多项式 rho_f, iphi
    ZZX rho_f, iphi;
    // 初始化大整数 Res_f
    ZZ Res_f;
    // 初始化 mod p的变量 Res_f_1
    ZZ_p Res_f_1;

    XGCD(Res_f, rho_f, iphi, f, phi, 0);
    // Res_f_1 = 1/Res_f
    inv(Res_f_1, conv<ZZ_p>(Res_f));
    // 断言，如果不 == 1报错提示
    assert(Res_f_1*conv<ZZ_p>(Res_f) == 1);
    // 返回
    return ( Res_f_1 * conv<ZZ_pX>(rho_f) );
}

//==============================================================================
//Computes h = g/f (mod phi) (mod q)
// 计算hash函数
//==============================================================================
ZZ_pX Quotient(const ZZX& f, const ZZX& g)
{
    // 初始化参数  f_1, g0, h0, phi0
    ZZ_pX f_1, g0, h0, phi0;
    // 求倒
    f_1 = Inverse(f);
    // g 转换为 ZZ_pX类型
    g0 = conv<ZZ_pX>(g);
    // phi 转换为 ZZ_pX类型
    phi0 = conv<ZZ_pX>(phi);
    // 计算h0
    h0 = (f_1*g0)%phi0;
    return h0;
}

//==============================================================================
//Verifies that for a parameter N, polynomials f, g are a valid semi-basis for building a NTRU lattice.
//验证参数N，多项式f，g是构建NTRU晶格的有效半基。
//If PGCD!=1, then (f,g) isn't a valid pair
//==============================================================================
void ValidPair(ZZ& PGCD, ZZ& Alpha, ZZ& Beta, ZZX& rho_f, ZZX& rho_g, const ZZX& f, const ZZX& g)
{
    // 初始化多项式 Res_fx, Res_gx, iphi
    ZZX Res_fx, Res_gx, iphi;
    // 初始化大整数  Res_f, Res_g
    ZZ Res_f, Res_g;

    // !看不懂，不知道干啥的
    // Res_f = resultant of f and phi;
    // if Res_f != 0, then computes rho_f and iphi such that: f*rho_f + phi*iphi = Res_f; otherwise rho_f and iphi not affected.
    // if !deterministic, then resultant computation may use a randomized
    // strategy that errs with probability no more than 2^{-80}.
    XGCD(Res_f, rho_f, iphi, f, phi, 0);

    // 计算两个整数最大公约数GCD(greatest common divisor)
    if(GCD(Res_f, q1)!=1)
    {
        PGCD = 0;
    }
    else
    {    XGCD(Res_g, rho_g, iphi, g, phi, 0);
         XGCD(PGCD, Alpha, Beta, Res_f, Res_g);
    }
}


//==============================================================================
//Computes the Gram-Schmidt norm of the basis B generated from f,g
// Gram-Schmidt 正交化
//==============================================================================
void GS_Norm(const ZZX fx, const ZZX gx, int& flag)
{
    // 初始化无符号整型 i
    unsigned int i;
    // 初始化doubel 类型 acc, acc3, Fred[N0], Gred[N0];
    double acc, acc3, Fred[N0], Gred[N0];
    // 初始化 CC_t复数类类型 f[N0], g[N0], F[N0], G[N0];
    CC_t f[N0], g[N0], F[N0], G[N0];

    // acc 赋值 0
    acc = 0;

    // 进行计算
    for(i=0; i < N0; i++)
    {
        acc += conv<double>(fx[i]*fx[i] + gx[i]*gx[i]);
    }

    // sqrt()函数返回数字 acc 的平方根
    acc = sqrt(acc);

    // 进行傅里叶变换。
    ZZXToFFT(f, fx);
    ZZXToFFT(g, gx);

    // 进行计算
    for(i=0; i<N0; i++)
    {
        F[i] = f[i]/(f[i]*f[N0-1-i]+g[i]*g[N0-1-i]);
        G[i] = g[i]/(f[i]*f[N0-1-i]+g[i]*g[N0-1-i]);
    }
    // 获取复数的实数部分
    MyRealReverseFFT(Fred, F);
    MyRealReverseFFT(Gred, G);

    // 赋值为0
    acc3 = 0;
    // 进行循环计算
    for(i=0; i<N0; i++)
    {
        acc3 += Fred[i]*Fred[i] + Gred[i]*Gred[i];
    }

    // acc3 = 1<<30 * acc3 的平方根。
    acc3 = q0 * sqrt(acc3);

    // 当 acc3<acc 的时候，将 flag 设置为1
    if(acc3<acc)
    {
        flag = 1;
    }
}


//==============================================================================
// Generates a secret basis (f,g),(F,G) from the parameters N,q,Norme
// This bases generates a NTRU lattice
// 生成一个 NTRU 格
//==============================================================================
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme)
{
    // 初始化整型变量 i
    int i;
    // 初始化多项式
    ZZX rho_f, rho_g, k, aux, fb, gb, num;
    // 初始 "big integers": signed, arbitrary length integers.
    ZZ PGCD, Alpha, Beta;

    // 初始化一个标签为0
    int flag = 0;

    while( (PGCD!=1) || (flag==0) )
    {
        flag = 1;
        // 生成固定次数和“近似”固定平方范数的随机多项式 f
        f = RandomPolyFixedSqNorm(Norme, N0-1);
        // 生成固定次数和“近似”固定平方范数的随机多项式 g
        g = RandomPolyFixedSqNorm(Norme, N0-1);
        // 计算由f，g生成的基B的Gram-Schmidt范数
        GS_Norm(f, g, flag);

        // !验证基
        ValidPair(PGCD, Alpha, Beta, rho_f, rho_g, f, g);
    }

    // 计算F
    F = -q1*Beta*rho_g;
    // 计算G
    G = q1*Alpha*rho_f;

    // 设置多项式f 长度
    f.SetLength(N0);
    // 设置多项式g 长度
    g.SetLength(N0);

    // 初始化无符号整型 mb
    unsigned int mb;

    // 计算 Reduce F and G: F ← F − k ∗ f, G ← G − k ∗ g
    k = ReductionCoefficient(f, g, F, G, mb);
    // 当 多项式的次数 >=0 的时候进行计算
    while(deg(k)>=0)
    {
        i++;
        // 计算 trunc(F-k*f, N0)
        F = FastMod(F - k*f);
        // 计算 trunc(G - k*g, N0)
        G = FastMod(G - k*g);

        // f 进行反转
        fb = Reverse(f);

        // f 进行反转
        gb = Reverse(g);

        // 计算 trunc(fb*F + gb*G, N0)
        num = FastMod(fb*F + gb*G);

        // 返回num系数的最大NumBits
        mb = MaxBits(num);

        // 计算 Reduce F and G: F ← F − k ∗ f, G ← G − k ∗ g
        k = ReductionCoefficient(f, g, F, G, mb);
        // 去除尾部所有的0,
        k.normalize();
    }
    // 计算 trunc((f*G - g*F)
    aux = FastMod(f*G - g*F);
    // 断言，如果 aux[0]==q1 或者  aux 次数为0 进行错误提示
    assert(aux[0]==q1);
    assert(deg(aux)==0);

    // cout << aux << endl;

    // 初始化长度
    aux.SetLength(N0);
    // cout << aux << endl;

}

// 向量点积
RR_t DotProduct(const RR_t * x1, const RR_t * x2)
{
    // 初始化无符号整型 i
    unsigned int i;
    // 初始化 RR_t类型 rep
    RR_t rep = 0;

    // 进行连加计算
    for(i=0; i<2*N0; i++)
    {
        rep += x1[i]*x2[i];
    }
    return rep;
}


void Rotate(RR_t * const dest, RR_t const * const src)
{
    unsigned int i;
    for(i=0; i<N0-1; i++)
    {
        dest[i+1] = src[i];
        dest[N0+i+1] = src[N0+i];
    }
    dest[0] = -src[N0-1];
    dest[N0] = -src[2*N0-1];
}



void ClassicMGS(RR_t Bstar[2*N0][2*N0], const RR_t B[2*N0][2*N0])
{
    RR_t SquareNorm[2*N0], aux[2*N0];
    unsigned int i,j,k;

    SquareNorm[0] = DotProduct(B[0], B[0]);
    for(j=0; j<2*N0; j++)
    {

        Bstar[0][j] = B[0][j];
    }

    for(i=1; i<2*N0; i++)
    {
        for(k=0; k<2*N0; k++)
        {
            Bstar[i][k] = B[i][k];
        }
        for(j=0; j<i; j++)
        {
            aux[j]= DotProduct(Bstar[i], Bstar[j]) / SquareNorm[j];
        }
        for(k=0; k<2*N0; k++)
        {
            for(j=0; j<i; j++)
            {
                Bstar[i][k] -= aux[j]*Bstar[j][k];
            }
        }
        SquareNorm[i] = DotProduct(Bstar[i], Bstar[i]);
    }
}

void FastMGS(RR_t Bst[2*N0][2*N0], const RR_t B[2*N0][2*N0])
{
    // 初始化 long double类型的变量 v[2*N0], v1[2*N0], C_k, D_k, C_ko, D_ko, aux
    RR_t v[2*N0], v1[2*N0], C_k, D_k, C_ko, D_ko, aux;
    //RR_t C[2*N0], D[2*N0];
    // 初始化无符号整型 j, k;
    unsigned int j, k;

    // 输出空行
    cout << endl;

    //Reducing first vector (obvious)
    for(j=0; j<2*N0; j++)
    {
        Bst[0][j] = B[0][j];
    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {
        v[j] = Bst[0][j+1];
        v[j+N0] = Bst[0][j+1+N0];
    }

    v[N0-1] = -Bst[0][0];
    v[2*N0-1] = -Bst[0][N0];

    for(j=0; j<2*N0; j++)
    {
        v1[j] = v[j];
    }

    //Initialising recurring variables
    C_k = DotProduct(Bst[0], v);
    //初始化重复变量
    D_k = DotProduct(v, v);

    //C[0] = C_k;
    //D[0] = D_k;
    //CD[0] = C[0]/D[0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=1; k<N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }

        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
        //printf ("C[%d]= %Lf		", k, C_k);
        //printf ("D[%d]= %Lf\n", k, D_k);
    }

    //Reducing second half!
    //cout << "aux = " << (1<<10)/D[N0-1] << endl;
    for(j=0; j<N0; j++)
    {    Bst[N0][N0+j] = Bst[N0-1][N0-1-j]*q0/D_k;
         Bst[N0][j] = -Bst[N0-1][2*N0-1-j]*q0/D_k;    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {    v[j] = Bst[N0][j+1];
         v[j+N0] = Bst[N0][j+1+N0];    }
    v[N0-1] = -Bst[N0][0];
    v[2*N0-1] = -Bst[N0][N0];

    for(j=0; j<2*N0; j++)
    {    v1[j] = v[j];    }


    //Initialising recursive variables
    C_k = DotProduct(Bst[N0], v1);
    D_k = DotProduct(Bst[N0], Bst[N0]);

    //C[N0] = C_k;
    //D[N0] = D_k;
    //CD[N0] = C[N0]/D[N0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=N0+1; k<2*N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }
        //SquareNorm[k] = SquareNorm[k-1] - aux*aux*sqnorm_v;


        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
    }
}
