// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#include "10X/Martian.h"
// order matters below -- gAlertTypes must be before the instance
const vec<std::string> Martian::gAlertTypes{"exit","alarm","log_info","log_warn"};
// static singleton instance
Martian Martian::gInst;
