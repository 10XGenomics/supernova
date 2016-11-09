// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_MISASSEMBLY_H
#define TENX_MISASSEMBLY_H

#include "CoreTools.h"
#include "math/HoInterval.h"

void Misassembly( const SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view,
     int64_t& total_err_dis_num, int64_t& total_err_dis_den,
     int64_t& total_err_ori_num, int64_t& total_err_ori_den,
     int64_t& total_err_ord_num, int64_t& total_err_ord_den );

#endif
