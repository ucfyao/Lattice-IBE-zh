/***
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-14 20:57:34
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-14 20:57:43
 * @FilePath     : /IBBE/io.h
 * @Description  :
 * @Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */

#ifndef IBBE_IO_H
#define IBBE_IO_H

#include "params.h"


void PrintDoubleTab(const unsigned int N, RR_t const * const Tab);
void PrintDoubleMatrix(const unsigned int N, const RR_t Tab[2*N0][2*N0]);
void PrintIntTab(const unsigned int N, long int const * const Tab);
void PrintFFT(CC_t const * const f_fft);


#endif
