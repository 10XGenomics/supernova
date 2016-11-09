///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Thread.h
 *
 *  Created on: Jul 31, 2014
 *      Author: tsharpe
 */
#ifndef SYSTEM_THREAD_H_
#define SYSTEM_THREAD_H_

int getTid();
bool isMainThread();
bool inParallelSection();

#endif /* SYSTEM_THREAD_H_ */
