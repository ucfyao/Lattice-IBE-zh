/*
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-11 10:18:13
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-14 21:01:33
 * @FilePath     : /IBBE/Scheme.cc
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
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    // 初始化大整数 ZZ 类型的变量 SqNorm
    ZZ SqNorm;
    // 初始化多项式 ZZX 类型变量 f,g,F,G
    ZZX f,g,F,G;
    // 计算 SqNorm = 1.36*q0/2
    SqNorm = conv<ZZ>(1.36*q0/2);
    // 生成格
    GenerateBasis(f, g, F, G, SqNorm);

    // 赋值到 PrivateKey
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
        // all vector elements are initialized to zero
        PrivateKey[i].SetLength(N0);
    }

    // 生成 PublicKey
    PublicKey = Quotient(f, g);
}


//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    // 初始化变量
    ZZX f,g,F,G;

    // 进行赋值
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;
    // 生成基
    B = BasisFromPolynomials(g, f, G, F);
}

// GPV签名算法。
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD)
{

    // 初始化 i
    int i;
    // 初始化 i
    unsigned j;
    // 初始化 long double 类型 ci[2*N0], zi, cip, sip, aux;
    RR_t ci[2*N0], zi, cip, sip, aux;

    // 进行赋值
    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }
    // !为什么空循环
    for(j=0; j<2*N0; j++)
    {

    }

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        // 向量点积 H(id)
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        // 随机取样
        zi = Sample4(cip, sip*PiPrime);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(MSKD->B)[i][j];
        }
    }

    // SKid ← s2
    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}

//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================

// 填充 MSK 数据
void CompleteMSK(MSK_Data * MSKD, ZZX * MSK)
{
    // 初始化无符号整数 i, j
    unsigned int i, j;
    // 初始化矩阵 B0
    mat_ZZ B0;

    // cout << MSKD << endl;

    // 赋值 MSK 数据
    for(i=0; i<4; i++)
    {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }
    cout << MSKD << endl;

    // 生成基B0
    CompletePrivateKey(B0, MSK);

    // 进行类型转换赋值
    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            MSKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    // for 1次。。？
    for(i=0; i<1; i++)
    {
        // 生成 Bstar
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    // GS_Norms 赋值
    for(i=0; i<2*N0; i++)
    {
        MSKD->GS_Norms[i] = sqrt( DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]) );
    }

    // sigma 赋值
    MSKD->sigma = 2*MSKD->GS_Norms[0];

}

// 填充 MPK 数据
void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK)
{
    MPKD->h = MPK;
    ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK));
}

// 提取 User secret key
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD)
{
    // 初始化循环变量 i
    unsigned int i;
    // 初始化 long double 类型
    RR_t c[2*N0], sk[2*N0], sigma;
    // 初始化多项式
    ZZX f,g,aux;

    // 赋值f
    f = MSKD -> PrK[0];
    // 赋值g
    g = MSKD -> PrK[1];
    // 赋值sigma
    sigma = MSKD->sigma;

    // 设置多项式长度
    SK_id[0].SetLength(N0);
    // 设置多项式长度
    SK_id[1].SetLength(N0);

    // cout << id << endl;
    // 将id 的类型转换为 RR_t
    for(i=0; i<N0 ;i++)
    {
        c[i] = ((RR_t) conv<double>(id[i])) ;
        c[i+N0] = 0;
    }
    // cout << c << endl;
    // 使用gpv算法生成sk
    GPV(sk, c, sigma, MSKD);

    // 将skid进行存储
    for(i=0; i<N0; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N0] = - sk[i+N0];
    }

    for(i=0; i<N0; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N0];
    }

}

// 对生成的 skid进行校验
unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD)
{
    // 初始化无符号整型
    unsigned int i;
    // 初始化多项式
    ZZX f,g,t,aux;

    // 进行赋值
    f = MSKD -> PrK[0];
    // 进行赋值
    g = MSKD -> PrK[1];

    // 转换id类型为多项式
    t = conv<ZZX>(id);
    // 生成辅助变量
    aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;

    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    // 进行判断
    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2) doesn't verify the required equality [ (s1 - t)*f + g*s2 = 0 ] !\nActually, (s1 - t)*f + g*s2 = " << aux << endl << endl;
    }

    // 返回
    return IsZero(aux);
}

// 消息加密
void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0], const MPK_Data * const MPKD)
{

    // 初始化无符号long类型变量
    unsigned long i;
    // 初始化long变量
    long r[N0], e1[N0], e2[N0];
    // 初始化复数类变量
    CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0];

    // 循环赋值随机数
    for(i=0; i<N0; i++)
    {
        e1[i] = (rand()%3) - 1;
        e2[i] = (rand()%3) - 1;
        r[i] = (rand()%3) - 1;
    }

    // fft 转换
    MyIntFFT(r_FFT, r);
    MyIntFFT(t_FFT, id0);

    for(i=0; i<N0; i++)
    {
        aux1_FFT[i] = r_FFT[i] * ((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i] * t_FFT[i];
    }

    // 进行逆变换
    MyIntReverseFFT(C[0], aux1_FFT);
    MyIntReverseFFT(C[1], aux2_FFT);

    // 密文计算
    for(i=0; i<N0; i++)
    {
        C[0][i] = (C[0][i] + e1[i]               + q0/2) % q0 - (q0/2);
        C[1][i] = (C[1][i] + e2[i] + (q0/2)*m[i] + q0/2) % q0 - (q0/2);
        // cout << "IBE_Encrypt C_0_i:"<< i << ":" << C[0][i] << endl;
        // cout << "IBE_Encrypt C_1_i:"<< i  << ":" << C[1][i] << endl;
    }
}

// 消息解密
void IBE_Decrypt(long message[N0], const long C[2][N0], const CC_t * const SKid_FFT)
{
    // 初始化无符号整数
    unsigned int i;
    // 初始化复数类
    CC_t c0_FFT[N0], aux_FFT[N0];

    // 进行FFT
    MyIntFFT(c0_FFT, C[0]);

    // 循环计算
    for(i=0; i<N0; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKid_FFT[i];
    }

    // 逆变换
    MyIntReverseFFT(message, aux_FFT);

    // 解密信息
    // 异或运算的意思是求两个运算分量相应位值是否相异,相异的为1,相同的为0
    for(i=0; i<N0; i++)
    {
        message[i] = C[1][i] - message[i];
        message[i] = ((unsigned long)(message[i] ))%q0;
        message[i] = (message[i] + (q0>>2) )/(q0>>1);
        message[i] %= 2;
    }

}



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD)
{
    // 初始化时间变量  t1, t2
    clock_t t1, t2;
    // 初始化 float 变量
    float diff;

    // 初始化循环变量 i
    unsigned int i;
    // 初始化向量id
    vec_ZZ id;
    // 初始化多项式
    ZZX SK_id[2];

    // 开始计时
    t1 = clock();

    // 开始输出
    cout << "0%" << flush;
    for(i=0; i < nb_extr; i++)
    {
        // 生成随机id
        id = RandomVector();
        // 生成 User secret key
        IBE_Extract(SK_id, id, MSKD);

        // 计算进度
        if((i+1) % (nb_extr/10) == 0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }
    // cout << "\nSKid" << endl;
    // cout << SK_id[0] << endl;
    // cout << SK_id[1] << endl;
    // 结束
    t2 = clock();
    // 计算时差
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
}


// 加密
void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    // 初始化时间变量
    clock_t t1, t2;
    // 初始化浮点类型变量 diff
    double diff;
    // 初始化无符号整型
    unsigned int i,j;
    // 初始化向量；
    vec_ZZ id;
    // 初始化多项式 SK_id，w
    ZZX SK_id[2], w;
    // 初始化复数类
    CC_t SKid_FFT[N0];
    // 初始化消息，解密
    long int message[N0], decrypted[N0];
    // 初始化身份，密文
    long int identity[N0], Ciphertext[2][N0];

    // 随机生成id
    id = RandomVector();

    // 提取skid
    IBE_Extract(SK_id, id, MSKD);

    // 验证key的有效性
    IBE_Verify_Key(SK_id, id, MSKD);

    // 傅里叶变换
    ZZXToFFT(SKid_FFT, SK_id[1]);
    // 将类型转换为long int
    for(i=0; i<N0; i++)
    {
        identity[i] = conv<long int>(id[i]);
    }

    // 开始时间
    t1 = clock();

    cout << "0%" << flush;
    // 进行多次加密解密
    for(i=0; i<nb_cryp; i++)
    {
        // 产生随机消息
        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }
        // 生成的随机信息
        // cout<< "random message:" << *message << endl;

        // 消息加密
        IBE_Encrypt(Ciphertext, message, identity, MPKD);

        // cout<< "random message Encrypt:" << *message << endl;
        // cout<< "random message Encrypt Ciphertext:" << *Ciphertext << endl;

        // 消息解密
        IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

        // cout << "random message Decrypt:" << *decrypted << endl;


        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    // 记录结束时间
    t2 = clock();
    // 计算时间差
    diff = ((double)t2 - (double)t1)/1000000.0l;
    // 格式化输出
    cout << "\n\nIt took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
    cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
    cout << "That's " << (diff/nb_cryp)*1000*1024/N0 << " milliseconds per encryption+decryption per Kilobit." << endl << endl;
}

// 测试
void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD)
{
    // 初始化无符号整型变量
    unsigned int i, rep;
    // 初始化向量
    vec_ZZ id;
    // 初始化多项式变量
    ZZX SK_id[2];

    // 进行赋值
    rep = 0;

    // 开始输出
    cout << "0%" << flush;

    for(i=0; i<nb_extr; i++)
    {
        // 生成随机id
        id = RandomVector();
        // 生成 User secret key
        IBE_Extract(SK_id, id, MSKD);
        // 进行验证
        rep += IBE_Verify_Key(SK_id, id, MSKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_extr << " extractions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
}

// 测试密钥对
void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    // 初始化无符号整型 i, j, rep;
    unsigned int i, j, rep;
    // 初始化向量id；
    vec_ZZ id;
    // 初始化多项式SK_id[2], m, w
    ZZX SK_id[2], m, w;
    // 初始化复数类
    CC_t SKid_FFT[N0];
    // 初始化id，密码本
    long int id0[N0], Ciphertext[2][N0];
    // 初始化消息，解密
    long int message[N0], decrypted[N0];

    // 随机生成id
    id = RandomVector();
    // 提取skid
    IBE_Extract(SK_id, id, MSKD);

    // 验证key的有效性
    IBE_Verify_Key(SK_id, id, MSKD);
    ZZXToFFT(SKid_FFT, SK_id[1]);

    rep = 0;

    for(i=0; i<N0; i++)
    {
        id0[i] = conv<long int>(id[i]);
    }

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {
        // 产生随机消息
        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }
        // 消息加密
        IBE_Encrypt(Ciphertext, message, id0, MPKD);
        // 消息解密
        IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

        // 进行信息对比
        for(j=0; j<N0; j++)
        {
            if(message[j] != decrypted[j])
            {
                cout << "ERROR : Dec(Enc(m)) != m " << endl;
                rep++;
                break;
            }
        }

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    cout << endl;
    // 输出结果
    if(rep == 0)
    {    cout << endl << nb_cryp << " encryptions+decryptions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_cryp << " encryptions+decryptions failed miserabily!" << endl << endl;    }
}
