///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Some read extension/scoring code pulled out of the anonymous namespace
// in BuildReadQGraph.

// Currently just static functions, but packaged together here for
// organizational purposes, and because I could see this expanding.
// At the very least, maybe some of the heuristics could be tunable.


#include "Basevector.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"


#ifndef EXTENDREADPATH_H_
#define EXTENDREADPATH_H_


class ExtendReadPath
{
public:
    // Paramters:
    //
    // mPenaltyDecay -  The scoring code keeps a running penalty that
    // builds up with errors in close proximity, so errors close
    // together are worse than errors spread throughout.  It's not
    // clear that that's ideal, but that's the current scoring model
    // and it was designed to deal with particular situations.
    //
    // The extra penalty from a "recent" error decreases by the fraction mPenaltyDecay
    // with each matching base seen.  So a mPenaltyDecay of 1.0 would cause a single
    // correct base to erase any previous penalty.  Note that this is not quite
    // the same as removing the penalty altogether.
    //
    // mMapQ2    - Q2 bases are assigned this value
    // mLeftOverPenalty - an additional score penalty for each base "left over"
    //           in the read -- i.e. the read is longer than the edge.
    //
    ExtendReadPath( HyperBasevector const& hbv, vec<int> const* to_left = nullptr,
            vec<int> const* to_right = nullptr, bool debug = false ) :
                    mHBV(hbv),
                    mpToLeft(to_left),
                    mpToRight(to_right),
                    mDebug(debug),
                    mManagedLeft(false),
                    mManagedRight(false),
                    mPenaltyDecay(0.2),
                    mMapQ2(20),
                    mLeftOverPenalty(10) {}

    ~ExtendReadPath() {
        if ( mManagedLeft  ) delete mpToLeft;
        if ( mManagedRight ) delete mpToRight;
    }

    void setPenaltyDecay(double p) {mPenaltyDecay=p;}
    double getPenaltyDecay() const { return mPenaltyDecay; }

    void setMapQ2(qual_t q) { mMapQ2 = q; }
    qual_t getMapQ2() const { return mMapQ2; }

    void setLeftOverPenalty(qual_t q) { mLeftOverPenalty = q; }
    qual_t getLeftOverPenalty() const { return mLeftOverPenalty; }

    void attemptLeftRightExtension( ReadPath& path,
                                    basevector const& bases, qualvector const& quals );

    bool attemptLeftwardExtension( ReadPath& path,
                                    basevector const& bases, qualvector const& quals );

    bool attemptRightwardExtension( ReadPath& path,
                                    basevector const& bases, qualvector const& quals );


    //////////////////////////////////////////////////////////////////////////
    // COMPATIBILITY INTERFACES
    // old compatibility interfaces use an anonymous object
    //////////////////////////////////////////////////////////////////////////
    static void attemptLeftRightExtension(  ReadPath& path, basevector const& bases,
                                            qualvector const& quals, HyperBasevector const& hbv,
                                            vec<int> const& to_left, vec<int> const& to_right,
                                            const bool debug = false ) {
        ExtendReadPath( hbv, &to_left, &to_right, debug).attemptLeftRightExtension(path,bases,quals);
    }

    static bool attemptLeftwardExtension(  ReadPath& path, basevector const& bases,
                                           qualvector const& quals, HyperBasevector const& hbv,
                                           vec<int> const& to_left, const bool debug = false ) {
        return ExtendReadPath( hbv, &to_left, nullptr, debug).attemptLeftwardExtension(path,bases,quals);
    }
  
    static bool attemptRightwardExtension( ReadPath& path, basevector const& bases,
                                           qualvector const& quals, HyperBasevector const& hbv,
                                           vec<int> const& to_right, const bool debug = false ) {
        return ExtendReadPath( hbv, nullptr, &to_right, debug).attemptRightwardExtension(path,bases,quals);
    }

private:
    // Manage mpToLeft and mpToRight so that they get created on-demand, if not
    // passed-in.  This is because the object may be instantiated with only
    // leftward or rightward (but not both) used.  This is to support the legacy
    // interfaces.
    vec<int> const& getToLeft() {
        if ( !mpToLeft ) {
            vec<int>* to_left = new vec<int>();
            mHBV.ToLeft(*to_left);
            mpToLeft = to_left;
            mManagedLeft = true;
        }
        return *mpToLeft;
    }
    vec<int> const& getToRight() {
        if ( !mpToRight ) {
            vec<int>* to_right = new vec<int>();
            mHBV.ToRight(*to_right);
            mpToRight = to_right;
            mManagedRight = true;
        }
        return *mpToRight;
    }

    HyperBasevector const& mHBV;
    vec<int> const* mpToLeft;
    vec<int> const* mpToRight;
    bool mDebug;

    bool mManagedLeft;
    bool mManagedRight;

    double mPenaltyDecay;
    qual_t mMapQ2;
    qual_t mLeftOverPenalty;
};


#endif /* EXTENDREADPATH_H_ */
