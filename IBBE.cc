/*
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-09 11:10:14
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-14 20:37:26
 * @FilePath     : /IBBE/IBBE.cc
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
// #include <NTL/ZZXFactoring.h>

#include "params.h"
#include "Algebra.h"
#include "Scheme.h"
// #include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"

using namespace std;
using namespace NTL;

int main(){
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept for our lattice-based IBBE scheme\n";
    cout << "The private key generation depends on the lattice basis delegation without increasing the dimension.";
    cout << "\n=======================================================================\n\n";

    // univariate polynomials over ZZ，初始化为具有4个元素的数组
    ZZX MSK[4];
    // univariate polynomials over ZZ_p
    ZZ_pX phiq, MPK;
    // 无符号整型
    unsigned int i;
    // float 4字节	数字介于 ±3.4E-38 和 ±3.4E38 之间，
    float diff;

    // new 申请一个新的对象，返回值为申请空间的对应数据类型的地址。
    MSK_Data * MSKD = new MSK_Data;
    MPK_Data * MPKD = new MPK_Data;

    // 计时函数、unsigned long 类型
    clock_t t1, t2;

    // 初始化一个ZZX的类，所使用的多项式环的次数为常量N0。
    const ZZX phi = Cyclo();

    // 可以看结构。
    // 输入 [2 10 14 6]
    // 因子：2
    // 多项式：[[[1 3] 1] [[1 1] 2]]
    // 公式：2 + 10*X + 14*x^2 +6*X^3 = 2 * (1 + 3*X) * (1 + X)^2

    // Vec< Pair< ZZX, long > > factors;
    // ZZ c;
    // factor(c, factors, phi);

    // cout << c << "\n";
    // cout << factors << "\n";

    // initialisation of rand
    srand(rdtsc());

    cout << "N = " << N0 << endl;
    cout << "q = " << q0 << endl;

    // The class ZZ_p is used to represent integers mod p.
    // The modulus p may be any positive integer, not necessarily prime.
    // Objects of the class ZZ_p are represented as a ZZ in the range 0..p-1.
    ZZ_p::init(q1);

    // The class zz_p is used to represent integers mod p, where 1 <= p <  NTL_SP_BOUND.
    // Note that NTL_SP_BOUND is usually 2^30 on 32-bit machines and 2^50 on 64-bit machines.
    // The modulus p may be any positive integer, not necessarily prime.
    // Objects of the class zz_p are represented as a long in the range 0..p-1.
    zz_p::init(q0);

    // 转换为 ZZ_pX 类型
    phiq = conv<ZZ_pX>(phi);
    // modular 运算
    ZZ_pXModulus PHI(phiq);

    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";

    // 获取开始时间 t1
    t1 = clock();

    // 生成 MPK, MSK
    for(i=0; i<1; i++)
    {
        Keygen(MPK, MSK);
    }
    // 填充msk的数据
    CompleteMSK(MSKD, MSK);

    // 填充mpk的数据
    CompleteMPK(MPKD, MPK);

    // 获取结束时间 t2
    t2 = clock();
    // 进行时间计算
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;


    //==============================================================================
    //Key extraction bench and encryption/decryption bench
    // 测试使用的时间
    //==============================================================================
    // const unsigned int nb_extrb = 100;
    // const unsigned int nb_crypb = 1000;

    // cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
    // cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
    // Extract_Bench(nb_extrb, MSKD);

    // // 加密
    // cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
    // cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
    // Encrypt_Bench(nb_crypb, MPKD, MSKD);



    ///==============================================================================
    // Key extraction test and encryption/decryption test
    // 加上了验证是否正确。
    //==============================================================================
    const unsigned int nb_extrt = 100;
    const unsigned int nb_crypt = 100;

    cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
    cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
    Extract_Test(nb_extrt, MSKD);

    cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
    cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
    Encrypt_Test(nb_crypt, MPKD, MSKD);

    free(MSKD);
    free(MPKD);
    return 0;

}