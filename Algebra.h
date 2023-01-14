/***
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-09 12:46:31
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-09 12:46:42
 * @FilePath     : /IBBE/Algebra.h
 * @Description  :
 * @Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */

// ifndef的用法在于避免重复包含和编译.
// 第一行ifndef回首先判断是否已经定义，如果已经定义将直接跳到#endif，否则将执行#define命令。
#ifndef IBBE_ALGEBRA_H
// define是C/C++中的宏定义，常用#define来定义常量。
#define IBBE_ALGEBRA_H

#include "params.h"

ZZX Cyclo();
ZZX FastMod(const ZZX& f);
ZZ SquaredNorm(const ZZX& f, const unsigned int degree);
void ValidPair(ZZ& PGCD, ZZ& Alpha, ZZ& Beta, ZZX& rho_f, ZZX& rho_g, const ZZX& f, const ZZX& g);
ZZX Reverse(const ZZX& f);
ZZX ReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, unsigned int & mb);
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G);
mat_ZZ AnticircularMatrix(const ZZX& f);
mat_ZZ BasisFromPolynomials(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G);
ZZ_pX Inverse(const ZZX& f);
ZZ_pX Quotient(const ZZX& f, const ZZX& g);
void GS_Norm(const ZZX fx, const ZZX gx, int& flag);
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme);
RR_t DotProduct(const RR_t * x1, const RR_t * x2);
void Rotate(RR_t * const dest, RR_t const * const src);
void ClassicMGS(RR_t Bstar[2*N0][2*N0], const RR_t B[2*N0][2*N0]);
void FastMGS(RR_t Bst[2*N0][2*N0], const RR_t B[2*N0][2*N0]);

#endif
