/***
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-11 10:28:59
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-11 10:29:07
 * @FilePath     : /IBBE/Random.h
 * @Description  :
 * @Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */

#ifndef IBBE_RANDOM_H
#define IBBE_RANDOM_H

#include "params.h"
#include "Sampling.h"

vec_ZZ RandomVector();
ZZX RandomPoly(const unsigned int degree);
ZZX RandomPolyFixedSqNorm(const ZZ& SqNorm, const unsigned int degree);

#endif
