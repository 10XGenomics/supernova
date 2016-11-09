///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PreCorrectAlt1.cc
 * \author tsharpe
 * \date Nov 29, 2012
 *
 * \brief
 */

#include "paths/long/PreCorrectAlt1.h"
#include "kmers/ReadPatherDefs.h"
#include "system/SysConf.h"
#include "system/WorklistN.h"

namespace
{
    unsigned const K = 31;

    class ReadFixer
    {
    public:
        ReadFixer( vecbvec* pReads, KmerDict<K> const& dict, int verbosity )
        : mpReads(pReads), mDict(dict),
          mDotter((pReads->size()+BATCH_SIZE-1)/BATCH_SIZE),
          mVerbosity(verbosity)
        { mCertified.resize(mpReads->size(),true);
          Spectrum<K> spectrum(dict);
          if ( mVerbosity >= 3 )
              std::cout << spectrum << std::endl;
          mCertThreshold = spectrum.getFirstTroughIndex();
          mCertThreshold = mCertThreshold + mCertThreshold/2; }

        // compiler-supplied destructor is OK

        size_t getNBatches() const { return mDotter.getNBatches(); }
        size_t getNReads() const { return mpReads->size(); }
        bool isCertified( size_t idx ) const { return mCertified[idx]; }

        class FixProc
        {
        public:
            FixProc( ReadFixer* pFixer ) : mpFixer(pFixer) {}

            // compiler-supplied copying and destructor are OK

            void operator()( size_t batchId )
            { using std::min;
              size_t nReads = mpFixer->getNReads();
              size_t start = min(batchId*BATCH_SIZE,nReads);
              size_t end = min(start+BATCH_SIZE,nReads);
              for ( ; start != end; ++start )
                mpFixer->fixRead(mpFixer->mpReads[0][start],start);
              mpFixer->mDotter.batchDone(); }

        private:
            ReadFixer* mpFixer;
        };
        friend struct FixProc;

        void fixAll( unsigned nThreads )
        { std::cout << Date() << ": Correcting " << getNReads() << " reads in "
                    << getNBatches() << " batches of " << BATCH_SIZE << '.'
                    << std::endl;
          parallelFor(0ul,mDotter.getNBatches(),FixProc(this),nThreads); }

    private:
        ReadFixer( ReadFixer const& ); // unimplemented -- no copying
        ReadFixer& operator=( ReadFixer const& ); // unimplemented -- no copying

        void fixRead( bvec& read, size_t readId )
        { if ( read.size() < K ) return;
          if ( mVerbosity >= 2 ) dumpMultiplicities(read,readId);
          if ( fixRead(read.begin(),read.end(),readId) )
              fixRead(read.rcbegin(),read.rcend(),readId); }

        void dumpMultiplicities( bvec const& read, size_t readId )
        {
            std::cout << "Read " << readId << std::endl;
            typedef bvec::const_iterator Itr;
            Itr beg(read.begin());
            Itr end(read.end());
            KMer<K> kkk(beg);
            KDef const* pDef = mDict.lookup(kkk);
            size_t count = pDef ? pDef->getCount() : 0;
            std::cout << count;
            size_t col = 1;
            for ( Itr itr(beg+K); itr != end; ++itr )
            {
                kkk.toSuccessor(*itr);
                pDef = mDict.lookup(kkk);
                count = pDef ? pDef->getCount() : 0;
                std::cout << ' ' << count;
                if ( !(++col % 40) )
                    std::cout << '\n';
            }
            std::cout << std::endl;
        }

        template <class Itr>
        bool fixRead( Itr beg, Itr const& end, size_t readId );

        static size_t const BATCH_SIZE = 10000;

        vecbvec* mpReads;
        vec<bool> mCertified;
        KmerDict<K> const& mDict;
        Dotter mDotter;
        unsigned mCertThreshold;
        int mVerbosity;
    };

    template <class Itr>
    bool ReadFixer::fixRead( Itr beg, Itr const& end, size_t readId )
    {
        unsigned const BAD_FACTOR = 5;
        unsigned const GOOD_FACTOR = 2;

        bool result = false;
        KMer<K> kkk(beg);
        KDef const* pDef = mDict.lookup(kkk);
        size_t count = pDef ? pDef->getCount() : 0;
        if ( count < mCertThreshold )
            mCertified[readId] = false;
        for ( Itr itr(beg+K); itr != end; ++itr )
        {
            kkk.toSuccessor(*itr);
            pDef = mDict.lookup(kkk);
            size_t nextCount = pDef ? pDef->getCount() : 0;

            if ( nextCount > BAD_FACTOR &&
                    count < mCertThreshold && BAD_FACTOR*count < nextCount )
                result = true; // reverse cliff detected

            // look for a place where the k-mer multiplicity falls off a cliff
            // the final base in that kmer is probably an error
            if ( count > BAD_FACTOR &&
                    nextCount < mCertThreshold && BAD_FACTOR*nextCount < count )
            {
                size_t bestCount = 0;
                unsigned bestCode = 0;
                size_t nextBestCount = 0;
                for ( unsigned code = 0; code != 4; ++code )
                {
                    kkk.setBack(code);
                    pDef = mDict.lookup(kkk);
                    if ( pDef )
                    {
                        size_t fixCount = pDef->getCount();
                        if ( fixCount > bestCount )
                        {
                            nextBestCount = bestCount;
                            bestCount = fixCount;
                            bestCode = code;
                        }
                        else if ( fixCount > nextBestCount )
                            nextBestCount = fixCount;
                    }
                }

                // if there's only one conceivable choice for the final base
                // in the kmer -- i.e., all others look like errors --
                // then correct the read to that base.
                if ( GOOD_FACTOR*bestCount >= count &&
                        BAD_FACTOR*nextBestCount < count )
                {
                    if ( mVerbosity > 1 )
                        std::cout << "read[" << readId << "][" << itr.pos()
                                    << "] sub " << Base::val2Char(bestCode)
                                    << " for " << Base::val2Char(*itr)
                                    << std::endl;
                    itr.set(bestCode);
                    nextCount = bestCount;
                }
                else if ( mVerbosity )
                    std::cout << "read[" << readId << "][" << itr.pos()
                                << "] uncorrected:" << count << ' '
                                << nextCount << ' ' << bestCount << ' '
                                << nextBestCount << std::endl;

                kkk.setBack(*itr);
            }
            count = nextCount;
            if ( count < mCertThreshold )
                mCertified[readId] = false;
        }
        return result;
    }
}

void precorrectAlt1( vecbvec* pReads, unsigned COVERAGE,
                        int VERBOSITY, unsigned NUM_THREADS )
{
    NUM_THREADS = boundNumThreads(NUM_THREADS);

    size_t estDictSize = 10*pReads->SizeSum()/COVERAGE;
    KmerDict<K> dict(estDictSize);
    if ( VERBOSITY )
        std::cout << "Creating dictionary." << std::endl;
    dict.process(*pReads,VERBOSITY,false,NUM_THREADS,10000);
    if ( VERBOSITY )
        std::cout << "There are " << dict.size()
                    << " dictionary entries (expected "
                    << estDictSize << ")." << std::endl;

    ReadFixer fixer(pReads,dict,VERBOSITY);
    fixer.fixAll(VERBOSITY?1:NUM_THREADS);
}
