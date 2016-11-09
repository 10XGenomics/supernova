///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * LargeKDispatcher.h
 *
 *  Created on: May 14, 2013
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_LARGEKDISPATCHER_H_
#define PATHS_LONG_LARGEKDISPATCHER_H_

#include "system/System.h"

class BigK
{
    static int constexpr gK[] =
    {20, 24, 28, 32, 40, 48, 60, 72, 80, 84, 88, 100, 128, 144, 172, 200, 320, 368, 400, 
     500, 544, 640, 720, 1000, 1200, 1600, 2000, 4000, 10000};
    // Case 1:  DANGER, WILL ROBINSON.  If you add a value to the array, you
    //          must also add another case statement at the end of the switch.
    //          If you forget, you won't get a dispatch on the last K value,
    //          you'll get a fatal err at runtime.  This would be bad.
    //          Keep the array size and the switch size in sync!!!
    // Case 2:  If you remove a value from the array, you'll have to remove
    //          the final case statement.  But if you forget, the compiler will
    //          warn you because you'll be referring to a non-existent array
    //          element.  So you'll fix it.
    // Case 3:  If you modify a value in the array, everything is fine.
    //          Nothing else changes.

public:

    // iterator pair over allowable values
    static int const* begin() { return gK; }
    static int const* end() { return gK+sizeof(gK)/sizeof(gK[0]); }

    // requirements:  a functor templated on K.  it may have any arg list,
    // and any return type, including void.
    // I.e.,
    // template <int K> struct FunctorT
    // { return_t operator()( arg1_t, arg2_t, etc ); };
    //
    // what it does:  instantiate and call the functor for the value of K (which
    // must be one of the allowable values) that you have passed as an arg.
    //
    // how you invoke it:
    // return_t ret = BigK::dispatch<FunctorT>(k, arg1, arg2, etc);
    //
    //
    //
    // NOTE(tim): Intel compiler fails to infer variadic template arguments
    // However, the variable return functionality was never actually used
    // anywhere, much easier to just make this function return void
    template <template <int K> class C, typename... Args>
    static void dispatch( int k, Args&&... args )
    {
        // NOTE(tim): return removed from each line of this switch, e.g.:
        // case gK[0]: return C<gK[0]>()(std::forward<Args>(args)...);
        switch (k)
        {
        case gK[0]: C<gK[0]>()(std::forward<Args>(args)...); break;
        case gK[1]: C<gK[1]>()(std::forward<Args>(args)...); break;
        case gK[2]: C<gK[2]>()(std::forward<Args>(args)...); break;
        case gK[3]: C<gK[3]>()(std::forward<Args>(args)...); break;
        case gK[4]: C<gK[4]>()(std::forward<Args>(args)...); break;
        case gK[5]: C<gK[5]>()(std::forward<Args>(args)...); break;
        case gK[6]: C<gK[6]>()(std::forward<Args>(args)...); break;
        case gK[7]: C<gK[7]>()(std::forward<Args>(args)...); break;
        case gK[8]: C<gK[8]>()(std::forward<Args>(args)...); break;
        case gK[9]: C<gK[9]>()(std::forward<Args>(args)...); break;
        case gK[10]: C<gK[10]>()(std::forward<Args>(args)...); break;
        case gK[11]: C<gK[11]>()(std::forward<Args>(args)...); break;
        case gK[12]: C<gK[12]>()(std::forward<Args>(args)...); break;
        case gK[13]: C<gK[13]>()(std::forward<Args>(args)...); break;
        case gK[14]: C<gK[14]>()(std::forward<Args>(args)...); break;
        case gK[15]: C<gK[15]>()(std::forward<Args>(args)...); break;
        case gK[16]: C<gK[16]>()(std::forward<Args>(args)...); break;
        case gK[17]: C<gK[17]>()(std::forward<Args>(args)...); break;
        case gK[18]: C<gK[18]>()(std::forward<Args>(args)...); break;
        case gK[19]: C<gK[19]>()(std::forward<Args>(args)...); break;
        case gK[20]: C<gK[20]>()(std::forward<Args>(args)...); break;
        case gK[21]: C<gK[21]>()(std::forward<Args>(args)...); break;
        case gK[22]: C<gK[22]>()(std::forward<Args>(args)...); break;
        case gK[23]: C<gK[23]>()(std::forward<Args>(args)...); break;
        case gK[24]: C<gK[24]>()(std::forward<Args>(args)...); break;
        case gK[25]: C<gK[25]>()(std::forward<Args>(args)...); break;
        case gK[26]: C<gK[26]>()(std::forward<Args>(args)...); break;
        case gK[27]: C<gK[27]>()(std::forward<Args>(args)...); break;
        case gK[28]: C<gK[28]>()(std::forward<Args>(args)...); break;
        default: FatalErr("Illegal value " << k << " for K in BigK dispatcher.");
        }
    }
};

#endif /* PATHS_LONG_LARGEKDISPATCHER_H_ */
