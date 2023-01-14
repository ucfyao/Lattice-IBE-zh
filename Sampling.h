/***
 * @Author       : yaozihao yaozihao@yaozihao.com
 * @Date         : 2023-01-11 10:32:51
 * @LastEditors  : yaozihao yaozihao@yaozihao.com
 * @LastEditTime : 2023-01-11 10:32:58
 * @FilePath     : /IBBE/Sampling.h
 * @Description  :
 * @Copyright (c) 2023 by yaozihao email: yaozihao@yaozihao.com, All Rights Reserved.
 */

#ifndef IBBE_SAMPLING_H
#define IBBE_SAMPLING_H

#include "params.h"


unsigned int Sample0(unsigned long alea);
unsigned int Sample1(const unsigned int k);
signed int Sample2(const unsigned int k);
signed int Sample3(const RR_t sigma128);
signed int Sample4(RR_t c, RR_t sigma);


#endif
