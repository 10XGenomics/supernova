///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/FillPairs.h"
#include "MapReduceEngine.h"
#include "ParallelVecUtilities.h"
#include "feudal/HashSet.h"
#include "kmers/KMer.h"
#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "kmers/naif_kmer/KmerMap.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/DiscovarTools.h"
#include "system/WorklistN.h"

namespace
{

unsigned const K = 60;
typedef KMer<K> Kmer;
typedef HashSet<Kmer,Kmer::Hasher> Dict;

class KHash
{
public:
  KHash( Dict* pDict, unsigned minFreq ) : mpDict(pDict), mMinFreq(minFreq) {}

  template <class OItr>
  void map( vecbvec::const_iterator itr, OItr oItr )
  { Kmer::kmerize(itr->begin(),itr->end(),oItr); }

  void reduce( Kmer const* beg, Kmer const* end )
  { if ( end-beg >= mMinFreq ) mpDict->add(*beg); }

  Kmer* overflow( Kmer* beg, Kmer* end )
  { return beg+std::min(unsigned(end-beg),mMinFreq); }

private:
  Dict* mpDict;
  unsigned mMinFreq;
};
typedef MapReduceEngine<KHash,Kmer,Kmer::Hasher> MRE;

void TrimReads( vecbvec const& reads, const int minFreq, vecbvec& trimmed )
{
    unsigned const COVERAGE = 50;
    size_t nKmers = reads.getKmerCount(K);

    // build a dictionary of the good kmers
    Dict dict(nKmers/COVERAGE);
    KHash impl(&dict,minFreq);
    MRE mre(impl);
    if ( !mre.run(nKmers,reads.begin(),reads.end()) )
        FatalErr("Failed to trim reads.  Ran out of key space in MRE.");

    // trim each read at end of 1st bad kmer
    trimmed = reads;
    parallelForBatch(trimmed.begin(),trimmed.end(),10000,
        [&dict]( vecbvec::iterator readsItr )
        {
            bvec& bv = *readsItr;
            if ( bv.size() < K )
                return;
            bvec::const_iterator beg(bv.cbegin());
            Kmer kkk(beg);
            Kmer krc(kkk); krc.rc();
            Kmer const* pEntry;
            if ( K&1 ) pEntry = dict.lookup(kkk.isRev() ? krc : kkk);
            else pEntry = dict.lookup(krc < kkk ? krc : kkk);
            if ( !pEntry )
            {
                bv.clear(); // kill the entire read if 1st kmer is bad
                return;
            }
            bvec::const_iterator itr(beg+K);
            bvec::const_iterator end(bv.cend());
            while ( itr != end )
            {
                kkk.toSuccessor(*itr);
                krc.toPredecessor(*itr ^ 3);
                if ( K&1 ) pEntry = dict.lookup(kkk.isRev() ? krc : kkk);
                else pEntry = dict.lookup(krc < kkk ? krc : kkk);
                if ( !pEntry )
                {
                    bv.resize(itr-beg); // itr points to end of bad kmer
                    break;
                }
                ++itr;
            }
        });
}

template<int K> void TrimReadsOld( const vecbasevector& bases,
                                    const int min_freq,
                                    const int NUM_THREADS,
                                    vecbasevector& basesx )
{
     ForceAssertLe(K, 60);      // only implmemented for K<=60 right now; need a wider Kmer below for higher K.

     typedef Kmer60 Kmer_t;
     typedef KmerKmerFreq<Kmer_t> KmerRec_t;            // Kmer frequency db
     vec<KmerRec_t> kmer_vec;

     // check that there are reads of longer than K bases
     size_t longs = 0, shorts = 0;
     for ( auto iter = bases.begin(); iter != bases.end(); ++iter ) {
         if ( iter->size() >= K ) longs++;
         else shorts++;
     }
     if ( longs == 0 ) {
         ostringstream err;
         err <<  "we found that all " << shorts << " reads were shorter than K=" << K;
         err << " and we therefore cannot continue.";
         DiscovarTools::ExitShortReads(err.str());
     }

     // calculate the kmer frequency db, thresholding at min_freq
     Validator valid( min_freq, 0 );
     KernelKmerStorer<KmerRec_t> storer( bases, K, &kmer_vec, &valid );
     naif_kmerize( &storer, NUM_THREADS, false );

     // populate the output before trimming
     basesx.assign(bases.begin(), bases.end());

     // calculate and apply trim if possible
     if ( kmer_vec.size() > 0 ) {       //edge case -- we may come back with no suprathreshold kmers

         // kmer map hash into the kmer db
         KmerMap<KmerRec_t> kmer_map( kmer_vec );

         vec<int> trim_to( bases.size( ), 0 );

         for ( size_t j = 0; j < bases.size(); j++ ) {      // for each read
             const basevector& b = bases[j];

             SubKmers<basevector, Kmer_t> subs( K, b );    // kmerize read b
             trim_to[j] = b.size();                        // default trim is full-sized
             bool first = true;
             for ( ; subs.not_done(); subs.next() ) {                   // for each kmer
                 if ( !kmer_map( subs.canonical() ).is_valid_kmer() ) {
                     // first kmer is special, because the error could be
                     // anywhere within it.  So we truncate the read to
                     // zero if the first kmer was low-freq.  After that,
                     // the newly introduced base is the culprit.
                     trim_to[j] = first ? 0 : subs.index_stop() - 1;    // found missing kmer, therefore low-frequency
                     break;
                 }
                 first = false;
             }
         }

         // trim reads
         for ( int64_t i = 0; i < (int64_t) bases.size( ); i++ )
              basesx[i].resize( trim_to[i] );
    }

}

} // close anonymous namespace

void FillPairs( const vecbasevector& bases, const PairsManager& pairs,
                    const int MIN_FREQ, vecbasevector& filled, bool newMethod,
                    bool useOldLRPMethod )
{
     // Build truncated reads.
     vecbasevector basesx;
//std::cout << Date() << ": FillPairs -- Trim reads." << std::endl; //TODO: remove
     if ( newMethod )
         TrimReads( bases, MIN_FREQ, basesx );
     else
         TrimReadsOld<60>( bases, MIN_FREQ, getConfiguredNumThreads(), basesx );

     // Path the trimmed reads.
//std::cout << Date() << ": FillPairs -- Path reads." << std::endl; //TODO: remove
     unsigned const COVERAGE = 50u;
     int verbosity = 0;
     HyperKmerPath h;
     vecKmerPath paths, paths_rc;
     HyperBasevector hb;
     LongReadsToPaths( basesx, K, COVERAGE, verbosity, useOldLRPMethod,
                         &hb, &h, &paths, &paths_rc );
     vecKmerPath hpaths;
     vec<tagged_rpint> hpathsdb;
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
          hpaths.push_back_reserve( h.EdgeObject(e) );
     CreateDatabase( hpaths, hpathsdb );

     // Close pairs.  Code copied from LoadAndCorrect.cc, which was copied from 
     // LongHyper.cc.

//std::cout << Date() << ": FillPairs -- Close pairs." << std::endl; //TODO: remove
     filled.resize(0);
     filled.resize( bases.size( ) );
     #pragma omp parallel for
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
     {    const int id2 = pairs.getPartnerID(id1);
          if ( id2 < id1 ) continue;
          vec< vec<int> > u(2);
          vec<int> left(2);

          for ( int pass = 0; pass < 2; pass++ )
          {    const KmerPath& p = ( pass == 0 ? paths[id1] : paths_rc[id2] );
               vec< triple<ho_interval,int,ho_interval> > M, M2;
               int rpos = 0;
               for ( int j = 0; j < p.NSegments( ); j++ )
               {    const KmerPathInterval& I = p.Segment(j);
                    vec<longlong> locs;
                    Contains( hpathsdb, I, locs );
                    for ( int l = 0; l < locs.isize( ); l++ )
                    {    const tagged_rpint& t = hpathsdb[ locs[l] ];
                         int hid = t.PathId( );
                         if ( hid < 0 ) continue;
                         longlong hpos = I.Start( ) - t.Start( );
                         longlong start = Max( I.Start( ), t.Start( ) );
                         longlong stop = Min( I.Stop( ), t.Stop( ) );
                         longlong hstart = start - t.Start( );
                         for ( int r = 0; r < t.PathPos( ); r++ )
                              hstart += hpaths[hid].Segment(r).Length( );
                         longlong hstop = hstart + stop - start;
                         longlong rstart = rpos + start - I.Start( );
                         longlong rstop = rstart + stop - start;
                         M.push( ho_interval( rstart, rstop ), hid,
                              ho_interval( hstart, hstop ) );    }
                    rpos += I.Length( );    }    
               Bool bad = False;
               for ( int i = 0; i < M.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < M.isize( ); j++ )
                    {    if ( M[j].first.Start( ) != M[j-1].first.Stop( ) + 1 ) 
                              break;
                         if ( M[j].second != M[j-1].second ) break;
                         if ( M[j].third.Start( ) != M[j-1].third.Stop( ) + 1 ) 
                              break;    }
                    u[pass].push_back( M[i].second );
                    Bool incomplete = False;
                    if ( i > 0 && M[i].third.Start( ) > 0 ) incomplete = True;
                    if ( j < M.isize( ) && M[j-1].third.Stop( )
                         != hpaths[ M[i].second ].KmerCount( ) - 1 )
                    {    incomplete = True;
                         bad = True;    }
                    if ( i == 0 && j == M.isize( ) && !incomplete )
                    {    i = j - 1;
                         continue;    }
                    int last = ( i == 0 ? -1 : M2.back( ).first.Stop( ) );
                    if ( M[i].first.Start( ) > last + 1 ) bad = True;
                    M2.push( ho_interval( M[i].first.Start( ),
                         M[j-1].first.Stop( ) ), M[i].second,
                         ho_interval( M[i].third.Start( ), M[j-1].third.Stop( ) ) );
                    if ( j == M.isize( ) && M[j-1].first.Stop( ) 
                         < p.KmerCount( ) - 1 )
                    {    bad = True;    }
                    i = j - 1;    }
               if (bad) u[pass].clear( );
               if ( u[pass].nonempty( ) )
                    left[pass] = M.front( ).third.Start( );    }
          if ( u[0].solo( ) && u[1].solo( ) && u[0][0] == u[1][0] )
          {    basevector b1 = basesx[id1], b2 = basesx[id2];
               b2.ReverseComplement( );
               int offset = left[1] - left[0];
               if ( offset >= 0 )
               {    basevector b( hb.EdgeObject( u[0][0] ),
                         left[0], left[1] - left[0] + b2.size( ) );
                    filled[id1] = b;
                    b.ReverseComplement( );
                    filled[id2] = b;    }    }    }    }
