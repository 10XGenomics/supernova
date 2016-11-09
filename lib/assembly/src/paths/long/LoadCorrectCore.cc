///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "paths/FindErrorsCore.h"
#include "paths/HyperKmerPath.h"
#include "kmers/KmerSpectrumCore.h"
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "lookup/LibInfo.h"
#include "lookup/SAM2CRD.h"
#include "paths/MergeReadSetsCore.h"
#include "paths/long/Correct1.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/CorrectPairs1.h"
#include "paths/long/DataSpec.h"
#include "paths/long/FillPairs.h"
#include "paths/long/Heuristics.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/PreCorrectAlt1.h"
#include "paths/long/PreCorrectOldNew.h"
#include "paths/long/DiscovarTools.h"
#include "random/Shuffle.h"
#include <numeric>
#include <type_traits>

void SamIAm( const int i, const String& getsam, const String& TMP, bool keepLocs, 
     String const& dexterLibs, Bool PF_ONLY, const Bool KEEP_NAMES )
{
    std::ofstream ofs((TMP+"/SAMProcessing"+ToString(i)+".log").c_str());
    ofs << "Running: " << getsam << endl;
    Logger logger(ofs);
    procbuf buf(getsam.c_str(),std::ios::in);
    istream is(&buf);
    SAM::SAMFile samfile(is,logger);
    vecbasevector seqs;
    vecqualvector quals;
    vec<look_align_x> alns;
    vec<pairinfo> pairsVec;
    vecString namesv;
    vecString libNames;
    bool use_OQ = true;
    vec<Bool> first_in_pair;
    SAM2CRD(samfile, seqs, quals, alns, pairsVec, namesv, first_in_pair, libNames,
                false, true, use_OQ, PF_ONLY, false, KEEP_NAMES );
    if (KEEP_NAMES)
    {    for ( int64_t i = 0; i < (int64_t) namesv.size( ); i++ )
         {    if ( first_in_pair[i] ) namesv[i] += ".1";
              else namesv[i] += ".2";    }    }
    int samtoolsResult = buf.close();
    if ( samtoolsResult )
    {
        ofs << "Error: samtools returned with non-zero status="
                << samtoolsResult << '\n';
        ofs.close();
        DiscovarTools::ExitSamtoolsFailed();
    }

    String outHead = TMP+'/'+ToString(i);
    if ( keepLocs )
    {
        Ofstream( aout, outHead + ".qltout");
        Ofstream( mapout, outHead + ".mapq");
        for ( auto itr=alns.begin(),end=alns.end(); itr != end; ++itr )
        {
            itr->PrintParseable(aout);
            mapout << (int)itr->mapQ() << "\n";
        }
    }

    int const SEP = -15;
    int const DEV = 12;
    int const NOMINAL_READ_LEN = 251;
    PairsManager pairs(seqs.size());
    if ( dexterLibs.empty() )
        for ( auto itr(libNames.begin()),end(libNames.end()); itr!=end; ++itr )
            pairs.addLibrary(SEP, DEV, *itr);
    else
    {
        LibInfoDB libInfoDB(dexterLibs);
        for ( auto itr(libNames.begin()),end(libNames.end()); itr!=end; ++itr )
        {
            LibInfo const* pInfo = libInfoDB.getInfo(*itr);
            if ( pInfo )
            {
                int sep = pInfo->mMean-2*NOMINAL_READ_LEN;
                pairs.addLibrary(sep,pInfo->mStdDev,*itr);
            }
            else
            {
                ofs << "Warning: Cannot find entry for library: "
                        << (*itr) << " in library information file: "
                        << dexterLibs << '\n';
                ofs << "Using default SEP and DEV.\n";
                pairs.addLibrary(SEP, DEV, *itr);
            }
        }
    }

    for ( auto itr(pairsVec.begin()),end(pairsVec.end()); itr != end; ++itr )
        pairs.addPairToLib(itr->readID1,itr->readID2,itr->libraryID,false);

    pairs.Write( outHead + ".pairs" );
    seqs.WriteAll(outHead + ".fastb");
    quals.WriteAll(outHead + ".qualb");
    if (KEEP_NAMES) namesv.WriteAll( outHead + ".names" );
    ofs.close();
}

void MergeReadSets( const vec<String>& heads, const String& TMP,
     const long_logging& logc, const Bool KEEP_NAMES )
{    if ( heads.empty( ) )
     {    cout << "\nInternal error, heads is empty." << endl;
          cout << "Please send us a report.\n" << endl;
          _exit(1);    }
     else if ( heads.solo( ) )
     {    System( "cp " + TMP + "/" + heads[0] + ".fastb "
               + TMP + "/frag_reads_orig.fastb" );
          System( "cp " + TMP + "/" + heads[0] + ".qualb "
               + TMP + "/frag_reads_orig.qualb" );
          System( "cp " + TMP + "/" + heads[0] + ".pairs "
               + TMP + "/frag_reads_orig.pairs" );
          if (KEEP_NAMES)
          {    System( "cp " + TMP + "/" + heads[0] + ".names "
                    + TMP + "/frag_reads_orig.names" );    }
          const String& qltout = TMP + "/" + heads[0] + ".qltout";
          if ( logc.KEEP_LOCS ) {
              if ( IsRegularFile( qltout ) )
                  System( "cp " + qltout + " " +
                        TMP + "/frag_reads_orig.qltout" );
              else
                  cout << "warning: missing aligns file: " << qltout << endl;
          }
     }
     else
     {    String head_out = TMP + "/frag_reads_orig";
          vec<String> headsp( heads.size( ) );
          for ( int i = 0; i < heads.isize( ); i++ )
               headsp[i] = TMP + "/" + heads[i];
          if (CheckFileSetExists( headsp, ".fastb", true ) == false)
               FatalErr("Fastb file missing - see above for details");
          vec<size_t> sizes;
          for ( size_t i = 0; i < heads.size( ); i++ )
               sizes.push_back(MastervecFileObjectCount( headsp[i] + ".fastb") );
          size_t size_sum = BigSum(sizes);
          MergeFeudal( head_out, headsp, ".fastb", sizes );
          MergeFeudal( head_out, headsp, ".qualb", sizes );
          if (KEEP_NAMES) MergeFeudal( head_out, headsp, ".names", sizes );
          PairsManager pairs;
          for ( size_t i = 0; i < heads.size( ); ++i )
          {    PairsManager pairs_loc;
               pairs_loc.Read( headsp[i] + ".pairs" );
               ForceAssertEq( pairs_loc.nReads( ), sizes[i] );
               pairs.Append(pairs_loc);    }
          pairs.Write( head_out + ".pairs" );
          if ( logc.KEEP_LOCS && CheckFileSetExists( headsp, ".qltout", true ) )
               MergeQltout( head_out, headsp, sizes );
      }  }

void SelectRandom( const String& TMP, const double SELECT_FRAC,
     const long_logging& logc, const long_data_spec& spec )
{
     String READS_IN = TMP + "/frag_reads_orig";
     String READS_OUT = TMP + "/frag_reads_orig";
     double FRAC = SELECT_FRAC;
     unsigned int SEED = spec.SELECT_SEED;
     Bool USE_LOCS = logc.KEEP_LOCS;
     Bool SHUFFLE = spec.SHUFFLE_INPUT;
     String in_bases_file = READS_IN + ".fastb";
     String in_quals_file = READS_IN + ".qualb";
     String in_names_file = READS_IN + ".names";
     String in_pairs_file = READS_IN + ".pairs";
     String in_qlt_file   = READS_IN + ".qltout";
     String out_bases_file  = READS_OUT + ".fastb";
     String out_quals_file  = READS_OUT + ".qualb";
     String out_names_file  = READS_OUT + ".names";
     String out_select_file = READS_OUT + ".select";
     String out_pairs_file  = READS_OUT + ".pairs";
     String out_qlt_file    = READS_OUT + ".qltout";
     Ofstream( log, TMP + "/SelectRandomPairs.log" );
     ForceAssertGt( FRAC, 0 );
     size_t n_reads = MastervecFileObjectCount( in_bases_file );
     vec<size_t> select;
     vec<longlong> maps_to( n_reads, -1 );
     {    log << Date( ) << ": loading pairing info" << endl;
          PairsManager pairs( in_pairs_file );
          size_t n_pairs = pairs.nPairs();
          log << Date( ) << ": loading bases" << endl;
          vecbvec bases( in_bases_file );
          log << Date( ) << ": shuffling pairs ids" << endl;
          vec<uint64_t> shuffled;
          if (SHUFFLE)
          {    if ( n_pairs > 0 )
                    Shuffle64( (uint64_t)n_pairs, shuffled, (uint64_t)SEED );    }
          else shuffled = vec<uint64_t>( n_pairs, vec<uint64_t>::IDENTITY );

          size_t n_keepers = 0;
          if (FRAC == 1.0)
          {    n_keepers = (n_pairs);
               size_t tot_length = 0;
               for (size_t ii = 0; ii < n_keepers; ii++)
               {    tot_length += (bases[pairs.ID1(shuffled[ii])].size() +
                         bases[pairs.ID2(shuffled[ii])].size());    }
               log << endl
                    << "total length of selected reads: " << tot_length << endl
                    << "pairs in input:                 " << n_pairs << endl
                    << "pairs selected:                 " << n_keepers << endl
                    << endl;    }
          else
          {    size_t required_length = 0;
               if (FRAC > 0)
               {    required_length
                         = uint64_t(round(double(bases.sumSizes()) * FRAC));    }
               size_t tot_length = 0;
               for (size_t ii=0; ii<n_pairs; ii++)
               {    if ( tot_length >= required_length ) break;
                    tot_length += bases[ pairs.ID1( shuffled[ii] ) ].size( );
                    tot_length += bases[ pairs.ID2( shuffled[ii] ) ].size( );
                    n_keepers = ii+1;    }
               log << endl
                    << "total length of selected reads: " << tot_length << endl
                    << "required total length:          " << required_length << endl
                    << "pairs in input:                 " << n_pairs << endl
                    << "pairs selected:                 " << n_keepers << endl
                    << endl;    }

          log << Date( ) << ": saving pairs" << endl;
          PairsManager sel_pairs( n_keepers * 2 );
          for( size_t lib_idx=0;lib_idx<pairs.nLibraries();++lib_idx){
              sel_pairs.addLibrary( pairs.getLibrarySep(lib_idx)
                                  , pairs.getLibrarySD(lib_idx)
                                  , pairs.getLibraryName(lib_idx));
          }
          for (size_t ii=0; ii<n_keepers; ii++)
          {    size_t id1 = select.isize();
               select.push_back( pairs.ID1( shuffled[ii] ) );
               maps_to[ pairs.ID1( shuffled[ii] ) ] = id1;
               size_t id2 = select.isize( );
               select.push_back( pairs.ID2( shuffled[ii] ) );
               maps_to[ pairs.ID2( shuffled[ii] ) ] = id2;
               sel_pairs.addPair( id1, id2, pairs.sep( shuffled[ii] ),
                    pairs.sd( shuffled[ii] ),
                    pairs.libraryName( shuffled[ii] ) );    }
          sel_pairs.Write( out_pairs_file );
          log << Date( ) << ": saving correspondence map (select)" << endl;
          BinaryWriter::writeFile( out_select_file.c_str(), select );
          log << Date( ) << ": saving bases" << endl;
          {    IncrementalWriter<bvec> sel_bases(out_bases_file.c_str());
               for (size_t ii=0; ii < select.size(); ii++)
                     sel_bases.add( bases[select[ii]] );    }    }

     if (IsRegularFile(in_quals_file))
     {    log << Date( ) << ": loading quals" << endl;
          vecqvec quals(in_quals_file);
          log << Date( ) << ": saving quals" << endl;
          IncrementalWriter<qvec> sel_quals(out_quals_file.c_str());
          for (size_t ii=0; ii < select.size(); ii++)
            sel_quals.add( quals[select[ii]] );
          sel_quals.close();    }

     if (IsRegularFile(in_names_file))
     {    log << Date( ) << ": loading names" << endl;
          vecString names(in_names_file);
          log << Date( ) << ": saving names" << endl;
          IncrementalWriter<String> sel_names(out_names_file.c_str());
          for (size_t ii=0; ii < select.size(); ii++)
            sel_names.add( names[select[ii]] );
          sel_names.close();    }

     // Select aligns (this could be improved: no need to load all aligns).
     if ( USE_LOCS && IsRegularFile(in_qlt_file) )
     {    log << Date() << ": loading alignments" << endl;
          vec<look_align> aligns_in;
          LoadLookAligns( in_qlt_file, aligns_in );
          log << Date( ) << ": selecting and saving alignments" << endl;
          vec<look_align> aligns_out;
          aligns_out.reserve( aligns_in.size( ) );
          for( size_t ii=0; ii<aligns_in.size( ); ii++)
          {    if ( maps_to[ aligns_in[ii].query_id ] < 0 ) continue;
               look_align al = aligns_in[ii];
               al.query_id = maps_to[ aligns_in[ii].query_id ];
               aligns_out.push_back( al );    }
          sort( aligns_out.begin( ), aligns_out.end( ) );
          WriteLookAligns( out_qlt_file, aligns_out );     }
     log << Date( ) << ": done selecting random pairs" << endl;

     int64_t n = MastervecFileObjectCount( TMP + "/frag_reads_orig.fastb" );
     if ( n == 0 ) DiscovarTools::ExitNoReads( );    }

void PopulateSpecials( const vecbasevector& creads, const PairsManager& pairs,
     const vecbasevector& creads_done, const vec<Bool>& done,
     const VecEFasta& corrected, const int NUM_THREADS, vec<Bool>& special,
     const long_logging& logc )
{
     if (logc.STATUS_LOGGING) cout << Date( ) << ": computing kmers_plus" << endl;
     const int M = 40;
     const int min_strong = 5;

     // only implmemented for K<=60 right now; need a wider Kmer below for higher K.
     ForceAssertLe(M, 60);

     typedef Kmer60 Kmer_t;
     typedef KmerKmerFreq<Kmer_t> KmerRec_t;            // Kmer + frequency
     vec<KmerRec_t> kmer_vec;

     // calculate the kmer frequency db, thresholding at min_freq
     Validator valid( min_strong, 0 );
     KernelKmerStorer<KmerRec_t> storer( creads, M, &kmer_vec, &valid );
     naif_kmerize( &storer, NUM_THREADS, false );

     // this is a kludge, but for now we're just going to convert the
     // kmers to the style expected by the code below.
     vec< kmer<M> > kmers;
     for ( size_t i = 0; i < kmer_vec.size(); ++i ) {
         basevector bv(M);
         for ( int j = 0; j < M; ++j )
              bv.Set(j, kmer_vec[i][j] );
         kmer<M> x(bv);
         kmers.push_back(x);
         x.ReverseComplement();
         kmers.push_back(x);
     }
     vec<KmerRec_t>().swap( kmer_vec );                 // return memory for kmer_vec.
     std::sort( kmers.begin(), kmers.end() );           // sort for later lookup, serial to avoid Parallel sort memory penalty

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": computing right extensions" << endl;
     const int min_ext = 200;
     vec<Bool> right_ext( kmers.size( ), False );
     #pragma omp parallel for
     for ( size_t id = 0; id < corrected.size( ); id++ )           // calculating right_ext[i] for kmers[i]
     {    vec<basevector> v;
          corrected[id].ExpandTo(v);
          if ( done[id] ) v.push_back( creads_done[id] );
          for ( int j = 0; j < v.isize( ); j++ )
          {    kmer<M> x;
               for ( int s = 0; s <= v[j].isize( ) - M; s++ )
               {    int ext = v[j].isize( ) - s;
                    x.SetToSubOf( v[j], s );
                    int64_t p;
                    if ( ext >= min_ext )
                    {    p = BinPosition( kmers, x );
                         if ( p >= 0 && !right_ext[p] )
                              right_ext[p] = True;    }
                    ext = s + M;
                    if ( ext >= min_ext )
                    {    x.ReverseComplement( );
                         p = BinPosition( kmers, x );
                         if ( p >= 0 && !right_ext[p] )
                              right_ext[p] = True;    }    }    }    }

     if (logc.STATUS_LOGGING) cout << Date( ) << ": finding specials" << endl;
     special.resize( creads.size(), False );
     vec< kmer<M> > fails;
     for ( int64_t i = 0; i < kmers.jsize( ); i++ )
          if ( !right_ext[i] ) fails.push_back( kmers[i] );           // fails is list of kmers for which !right_ext[]
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
     {    const int64_t idp = pairs.getPartnerID(id);
          kmer<M> x;
          for ( int s = 0; s <= creads[id].isize( ) - M; s++ )
          {    x.SetToSubOf( creads[id], s );                         // kmerize read and populate "special"
               if ( BinMember( fails, x ) )
               {
                    special[id] = True;
                    special[idp] = True;    }
               if ( s + M >= min_ext )
               {    x.ReverseComplement( );
                    if ( BinMember( fails, x ) )
                    {
                         special[id] = True;
                         special[idp] = True;    }    }    }    }
}

// to switch between VirtualMasterVec<bvec> and vecbasevector. The former cannot be const&
template<typename T>
struct ZeroCorrectedQuals_impl{
    static void do_it(T oreads, vecbvec const& creads, vecqvec* pQuals){
        vecqvec& cquals = *pQuals;
        ForceAssertEq(oreads.size(),creads.size());
        ForceAssertEq(oreads.size(),cquals.size());
        auto iOBV = oreads.begin();
        auto iCBV = creads.begin();
        auto iEnd = cquals.end();
        for ( auto iQV = cquals.begin(); iQV != iEnd; ++iOBV,++iCBV,++iQV )
        {
            bvec const& bvOrig = *iOBV;
            bvec const& bvCorr = *iCBV;
            qvec& qv = *iQV;
            ForceAssertEq(bvOrig.size(), bvCorr.size());
            ForceAssertEq(bvOrig.size(), qv.size());
            auto iO = bvOrig.cbegin();
            auto iC = bvCorr.cbegin();
            for ( auto iQ = qv.begin(), iE = qv.end(); iQ != iE; ++iQ,++iO,++iC )
                if ( *iO != *iC )
                    *iQ = 0;
        }
    }
};
// zero all quality scores associated with corrections (for reads in a file)
void ZeroCorrectedQuals( String const& readsFile, vecbvec const& creads,
                            vecqvec* pQuals )
{
    VirtualMasterVec<bvec> oreads(readsFile);
    ZeroCorrectedQuals_impl<typename std::add_lvalue_reference<decltype(oreads)>::type>::do_it(oreads,creads,pQuals);
}
// zero all quality scores associated with corrections (for reads already in memory)
void ZeroCorrectedQuals( vecbasevector const& oreads, vecbvec const& creads,
                            vecqvec* pQuals )
{
    ZeroCorrectedQuals_impl<decltype(oreads)>::do_it(oreads,creads,pQuals);
}

void CapQualityScores( vecqualvector& cquals, const vec<Bool>& done )
{    const int cap_radius = 4;
     for ( int64_t id = 0; id < (int64_t) cquals.size( ); id++ )
     {    if ( done[id] ) continue;
          vec<int> q( cquals[id].size( ), 1000000 );
          for ( int j = 0; j < (int) cquals[id].size( ); j++ )
          {    int start = Max( 0, j - cap_radius );
               int stop = Min( (int) cquals[id].size( ) - 1, j + cap_radius );
               for ( int l = start; l <= stop; l++ )
                    q[j] = Min( q[j], (int) cquals[id][l] );    }
          for ( int j = 0; j < (int) cquals[id].size( ); j++ )
               cquals[id][j] = q[j];    }    }

void CorrectionSuite( const String& TMP, const long_heuristics& heur,
     const long_logging& logc, const long_logging_control& log_control,
     vecbasevector& creads, VecEFasta& corrected, vec<int>& cid,
     vec<pairing_info>& cpartner, const uint NUM_THREADS, const String& EXIT,
     const double clock, bool useOldLRPMethod )
{
    LongProtoTmpDirManager tmp_mgr(TMP);
    CorrectionSuite( tmp_mgr, heur, logc, log_control, creads, corrected, cid, cpartner, NUM_THREADS, EXIT, clock, useOldLRPMethod );
} ;

void CorrectionSuite( LongProtoTmpDirManager& tmp_mgr, const long_heuristics& heur,
     const long_logging& logc, const long_logging_control& log_control,
     vecbasevector& creads,
     VecEFasta& corrected, vec<int>& cid, vec<pairing_info>& cpartner,
     const uint NUM_THREADS, const String& EXIT, const double clock,
     bool useOldLRPMethod )
{
     // Run Correct1.

     vec<int> trace_ids, precorrect_seq;
     ParseIntSet( logc.TRACE_IDS, trace_ids );

     vec<int> trace_pids;
     ParseIntSet( logc.TRACE_PIDS, trace_pids );
     for ( int i = 0; i < trace_pids.isize( ); i++ )
     {    int64_t pid = trace_pids[i];
          trace_ids.push_back( 2*pid, 2*pid + 1 );    }

     ParseIntSet( "{" + heur.PRECORRECT_SEQ + "}", precorrect_seq, false );
     const int max_freq = heur.FF_MAX_FREQ;

     const String sFragReadsOrig="frag_reads_orig";
     const String sFragReadsMod0="frag_reads_mod0";
     const bool bOrgReadsInMem = tmp_mgr[sFragReadsOrig].inMem();

     {    double bclock = WallClockTime( );
          vecqualvector cquals;
          if( bOrgReadsInMem ){
              creads = tmp_mgr[sFragReadsOrig].reads();
              cquals = tmp_mgr[sFragReadsOrig].quals();
          }
          else{
              creads.ReadAll( tmp_mgr.dir() + "/frag_reads_orig.fastb" );
              cquals.ReadAll( tmp_mgr.dir() + "/frag_reads_orig.qualb" );
          }
          size_t nReads = creads.size();
          ForceAssertEq(nReads,cquals.size());
          size_t nBases = 0, qualSum = 0;
          for ( qvec const& qv : cquals )
          {   nBases += qv.size();
              qualSum = std::accumulate(qv.begin(),qv.end(),qualSum);   }
          ForceAssertEq(nBases,creads.SizeSum());

          if (logc.MIN_LOGGING)
          {    cout << Date( ) << ": there are " << ToStringAddCommas(nReads)
                    << " reads" << endl;
               cout << Date( ) << ": mean read length = " << setiosflags(ios::fixed)
                    << setprecision(1) << double(nBases)/nReads
                    << resetiosflags(ios::fixed) << endl;
               cout << Date( ) << ": mean base quality = " << setiosflags(ios::fixed)
                    << setprecision(1) << double(qualSum)/nBases
                    << resetiosflags(ios::fixed) << endl;    }

          if (logc.STATUS_LOGGING) ReportPeakMem( "start precorrection" );
          if ( heur.PRECORRECT_ALT1 )
              precorrectAlt1(&creads);
          else if (heur.PRECORRECT_OLD_NEW)
              PreCorrectOldNew( &creads, cquals, trace_ids );
          else
          {   PC_Params pcp;
              const int K_PC = 25;
              KmerSpectrum kspec(K_PC);
              pre_correct_parallel( pcp, K_PC, &creads, &cquals, &kspec,
                                        -1, NUM_THREADS );    }

          if( bOrgReadsInMem ){
              ZeroCorrectedQuals(tmp_mgr[sFragReadsOrig].reads(),creads,&cquals);
          }
          else{
              ZeroCorrectedQuals(tmp_mgr.dir()+"/frag_reads_orig.fastb",creads,&cquals);
          }


          if ( logc.DUMP_PRE )
          {    creads.WriteAll( tmp_mgr.dir() + "/frag_reads_pre.fastb" );
               cquals.WriteAll( tmp_mgr.dir() + "/frag_reads_pre.qualb" );    }
          if (logc.STATUS_LOGGING) ReportPeakMem("precorrection done");
          REPORT_TIME( bclock, "used in initial precorrection" );

          if ( (*log_control.G).size( ) > 0 )
          {    double true_size = 0;
               for ( int g = 0; g < (int) (*log_control.G).size( ); g++ )
               {    true_size += (*log_control.G)[g].size( )
                         * (*log_control.ploidy)[g];    }
               cout << Date( ) << ": nominal coverage = "
                    << setiosflags(ios::fixed) << setprecision(1)
                    << double(nBases)/true_size
                    << resetiosflags(ios::fixed) << "x" << endl;    }
          if ( EXIT == "NOMINAL_COV" ) Done(clock);

          PairsManager const& pairs = tmp_mgr[sFragReadsOrig].pairs();
//          pairs.Read( TMP + "/frag_reads_orig.pairs" );
          pairs.makeCache( );

          // Carry out initial pair filling.

          vecbasevector creads_done;
          vec<Bool> to_edit( nReads, True );
          vec<Bool> done( nReads, False );
          if (heur.CORRECT_PAIRS)
          {    double fclock = WallClockTime( );
               creads_done = creads;
               vecbasevector filled;
               if (logc.STATUS_LOGGING)
                    cout << Date( ) << ": start initial pair filling" << endl;
               const int MIN_FREQ = 5;
               FillPairs( creads, pairs, MIN_FREQ, filled, heur.FILL_PAIRS_ALT, useOldLRPMethod );
               REPORT_TIME( fclock, "used in FillPairs" );
               double f2clock = WallClockTime( );
               int64_t fill_count = 0;
               for ( int64_t id = 0; id < (int64_t) filled.size( ); id++ )
               {    if ( filled[id].size( ) == 0 ) continue;
                    fill_count++;
                    int n = creads[id].size( );
                    creads_done[id] = filled[id];
                    cquals[id].resize(0);
                    cquals[id].resize( filled[id].size( ), 40 );
                    creads[id] = creads_done[id];
                    if ( n < creads[id].isize( ) )
                    {    cquals[id].resize(n);
                         if ( pairs.getPartnerID(id) >= id ) creads[id].resize(n);
                         else
                         {    creads[id].SetToSubOf( creads[id],
                                   creads[id].isize( ) - n, n );    }    }

                    done[id] = True;
                    if ( pairs.getPartnerID(id) < id )
                    {    creads_done[id].resize(0);    }
                    to_edit[id] = False;    }
               if (logc.STATUS_LOGGING)
               {    cout << Date( ) << ": "
                         << PERCENT_RATIO( 3, fill_count, (int64_t) filled.size( ) )
                         << " of pairs filled" << endl;    }
               REPORT_TIME( f2clock, "used in filling tail" );    }

          // New precorrection.

          vec<int> trim_to;
          if (heur.CORRECT_PAIRS)
          {    double mclock = WallClockTime( );
               if (logc.STATUS_LOGGING)
                    cout << Date( ) << ": begin new precorrection" << endl;

               // Cap quality scores.

               CapQualityScores( cquals, done );
               REPORT_TIME( mclock, "used capping" );

               // Do precorrection.

               for ( int j = 0; j < precorrect_seq.isize( ); j++ )
               {    Correct1Pre( tmp_mgr.dir(), precorrect_seq[j], max_freq, creads, cquals,
                         pairs, to_edit, trim_to, trace_ids, logc, heur );    }
               /*
               Correct1( 40, max_freq, creads, cquals, pairs, to_edit, trim_to,
                    trace_ids, log_control, logc );
               */

               double nclock = WallClockTime( );
               if(bOrgReadsInMem){
                   tmp_mgr[sFragReadsMod0].reads(true)=creads;
                   tmp_mgr[sFragReadsMod0].quals(true)=cquals;
                   tmp_mgr[sFragReadsMod0].pairs(true);

               }
               else{
                   creads.WriteAll( tmp_mgr.dir() + "/frag_reads_mod0.fastb" );
                   cquals.WriteAll( tmp_mgr.dir() + "/frag_reads_mod0.qualb" );
               }

               // Path the precorrected reads.

               if (logc.STATUS_LOGGING)
                    cout << Date( ) << ": pathing precorrected reads" << endl;
               unsigned const COVERAGE = 50u;
               const int K2 = 80; // SHOULD NOT BE HARDCODED!
               vecbasevector correctedv(creads);
               for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
                    correctedv[id].resize( trim_to[id] );
               HyperBasevector hb;
               HyperKmerPath h;
               vecKmerPath paths, paths_rc;
               LongReadsToPaths( correctedv, K2, COVERAGE, logc.verb[ "LRP" ],
                                  useOldLRPMethod, &hb, &h, &paths, &paths_rc );
               vecKmerPath hpaths;
               vec<tagged_rpint> hpathsdb;
               for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
                    hpaths.push_back_reserve( h.EdgeObject(e) );
               CreateDatabase( hpaths, hpathsdb );
               if (logc.STATUS_LOGGING) cout << Date( ) << ": done" << endl;

               // Close pairs that we're done with.  Code copied with minor
               // changes from LongHyper.cc.  Should be completely rewritten.

               if (logc.STATUS_LOGGING)
                    cout << Date( ) << ": initially closing pairs" << endl;
               #pragma omp parallel for
               for ( int64_t id1 = 0; id1 < (int64_t) creads.size( ); id1++ )
               {    if ( done[id1] ) continue;
                    const int id2 = pairs.getPartnerID(id1);
                    if ( id2 < id1 ) continue;
                    vec< vec<int> > u(2);
                    vec<int> left(2);

                    for ( int pass = 0; pass < 2; pass++ )
                    {    const KmerPath& p
                              = ( pass == 0 ? paths[id1] : paths_rc[id2] );
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
                              {    if ( M[j].first.Start( )
                                        != M[j-1].first.Stop( ) + 1 )
                                   {    break;    }
                                   if ( M[j].second != M[j-1].second ) break;
                                        if ( M[j].third.Start( )
                                        != M[j-1].third.Stop( ) + 1 )
                                   {    break;    }    }
                              u[pass].push_back( M[i].second );
                              Bool incomplete = False;
                              if ( i > 0 && M[i].third.Start( ) > 0 )
                                   incomplete = True;
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
                                   ho_interval( M[i].third.Start( ),
                                   M[j-1].third.Stop( ) ) );
                              if ( j == M.isize( ) && M[j-1].first.Stop( )
                                   < p.KmerCount( ) - 1 )
                              {    bad = True;    }
                              i = j - 1;    }
                         if (bad) u[pass].clear( );
                         if ( u[pass].nonempty( ) )
                         {    left[pass] = M.front( ).third.Start( );    }    }
                    if ( u[0].solo( ) && u[1].solo( ) && u[0][0] == u[1][0] )
                    {    int b1siz = correctedv[id1].isize();
                         int b2siz = correctedv[id2].isize();
                         int offset = left[1] - left[0];
                         if ( b1siz == creads[id1].isize( )
                              && b2siz == creads[id2].isize( ) && offset >= 0 )
                         {    auto beg = hb.EdgeObject(u[0][0]).begin()+left[0];
                              auto end = beg+(left[1]-left[0]+b2siz);
                              creads_done[id1].assign(beg,end);
                              creads_done[id2] = creads_done[id1];
                              creads_done[id2].ReverseComplement();

                              creads[id1] = creads_done[id1];
                              creads[id1].resize( b1siz );
                              creads[id2] = creads_done[id2];

                              creads[id2].SetToSubOf( creads[id2],
                                   creads[id2].isize( ) - b2siz, b2siz );
                              cquals[id1].resize(0);
                              cquals[id1].resize( creads[id1].size( ), 40 );
                              cquals[id2].resize(0);
                              cquals[id2].resize( creads[id2].size( ), 40 );

                              done[id1] = done[id2] = True;
                              creads_done[id2].resize(0);
                              to_edit[id1] = False;
                              to_edit[id2] = False;    }    }    }
               if (logc.STATUS_LOGGING)
               {    cout << Date( ) << ": "
                         << PERCENT_RATIO( 3, Sum(done), done.isize( ) )
                         << " of pairs were preclosed" << endl;    }
               REPORT_TIME( nclock, "used in main precorrection tail" );    }

          if (heur.CORRECT_PAIRS)
          {    corrected.clear().resize( creads.size( ) );
               CorrectPairs1( tmp_mgr.dir(), 40, max_freq, creads, cquals, pairs, to_edit,
                    trace_ids, heur, log_control, logc, corrected );
               for ( size_t id = 0; id < corrected.size( ); id++ )
               {    if ( corrected[id].size( ) > 0 )
                    {    to_edit[id] = False;
                         const int64_t idp = pairs.getPartnerID(id);
                         to_edit[idp] = False;    }    }

               if (heur.CP2)
               {
               double cp2_clock = WallClockTime( );
               vec<Bool> special;
               PopulateSpecials( creads, pairs, creads_done, done, corrected,
                    NUM_THREADS, special, logc );

               for ( size_t id = 0; id < corrected.size( ); id++ )
                    if ( !special[id] ) to_edit[id] = False;
               if (logc.STATUS_LOGGING)
               {    cout << Date( ) << ": second pass of CorrectPairs1 to use "
                         << Sum(to_edit)/2 << " pairs" << endl;    }

               if ( logc.verb[ "CP2" ] >= 1 )
               {    cout << "\nCP2, using ids:\n";
                    for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
                         if ( to_edit[id] ) cout << id << endl;
                    cout << endl;    }

               long_heuristics heur2(heur);
               // heur2.CP_MIN_GLUE = 5;
               heur2.CP_MIN_GLUE = 15;
               heur2.CP_MINQ_FLOOR = 0;
               heur2.CP_RAISE_ZERO = True;
               heur2.CP_MAX_QDIFF = 25.0;

               if (logc.STATUS_LOGGING)
               {    cout << Date( ) << ": running second pass of CorrectPairs1, "
                         << "using " << Sum(to_edit)/2 << " pairs" << endl;    }
               REPORT_TIME( cp2_clock, "used in prep for CP2" );

               CorrectPairs1( tmp_mgr.dir(), 40, max_freq, creads, cquals, pairs, to_edit,
                    trace_ids, heur2, log_control, logc, corrected );
               } // end of heur.CP2

               double pclock = WallClockTime( );
               for ( int64_t id = 0; id < done.jsize( ); id++ )
               {    if ( done[id] ) corrected[id] = creads_done[id];    }
               REPORT_TIME( pclock, "used in pair correction copying" );    }

          else // Non-default!
          {    double oclock = WallClockTime( );
               Correct1( tmp_mgr.dir(), 40, max_freq, creads, cquals, pairs, to_edit, trim_to,
                    trace_ids, log_control, logc, heur );
               for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
                    if ( trim_to[id] == creads[id].isize( ) ) to_edit[id] = False;
               Correct1( tmp_mgr.dir(), 60, max_freq, creads, cquals, pairs, to_edit, trim_to,
                    trace_ids, log_control, logc, heur );
               for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
                    if ( trim_to[id] == creads[id].isize( ) ) to_edit[id] = False;
               Correct1( tmp_mgr.dir(), 80, max_freq, creads, cquals, pairs, to_edit, trim_to,
                    trace_ids, log_control, logc, heur );
               for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
                    if ( trim_to[id] == creads[id].isize( ) ) to_edit[id] = False;
               Correct1( tmp_mgr.dir(), 28, max_freq, creads, cquals, pairs, to_edit, trim_to,
                    trace_ids, log_control, logc, heur );
               for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
               {    creads[id].resize( Max( 1, trim_to[id] ) );
                    cquals[id].resize( Max( 1, trim_to[id] ) );    }

               // Convert to efasta.

               corrected.assign(creads.begin(),creads.end());
               REPORT_TIME( oclock, "used in old unpair correction" );    }

          // Write corrected pairs.

          double wclock = WallClockTime( );
          if(!bOrgReadsInMem){
               vecbasevector correctedb( creads.size( ) );
               Ofstream( out, tmp_mgr.dir() + "/frag_reads_mod.efasta" );
               for ( size_t i = 0; i < corrected.size( ); i++ )
               {    corrected[i].FlattenTo( correctedb[i] );
                    corrected[i].Print( out, ToString(i) );    }
               correctedb.WriteAll( tmp_mgr.dir() + "/frag_reads_mod.fastb" );
          }
          REPORT_TIME( wclock, "used writing corrected" );    }    }

// Define pairing info.  Note that for now we set all the library ids to 0.

void DefinePairingInfo( const LongProtoTmpDirManager& tmp_mgr, const vecbasevector& creads,
     const vec<Bool>& to_delete, vec<int>& cid, VecEFasta& corrected,
     vec<pairing_info>& cpartner, const long_logging& logc )
{    double clock = WallClockTime( );
     PairsManager const& pairs = tmp_mgr.get("frag_reads_orig").pairs();
//     pairs.Read( TMP + "/frag_reads_orig.pairs" );
     for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
          if ( !to_delete[id] ) cid.push_back(id);
     corrected.EraseIf(to_delete);
     cpartner.resize( creads.size( ) );
     for ( int64_t xid1 = 0; xid1 < (int64_t) corrected.size( ); xid1++ )
     {    int id1 = cid[xid1];
          int64_t pid = pairs.getPairID(id1);
          if ( pairs.isUnpaired(id1) ) cpartner[xid1] = pairing_info(0,-1,-1);
          else
          {    int id2 = pairs.getPartnerID(id1);
               int xid2 = BinPosition( cid, id2 );
               if ( xid2 < 0 ) cpartner[xid1] = pairing_info(0,-1,-1);
               else
               {    if ( pairs.ID1(pid) == id1 )
                         cpartner[xid1] = pairing_info(1,xid2,0);
                    else cpartner[xid1] = pairing_info(2,xid2,0);    }    }    }
     REPORT_TIME( clock, "used in load tail" );    }

