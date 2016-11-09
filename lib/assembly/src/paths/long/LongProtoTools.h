///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_PROTO_TOOLS_H
#define LONG_PROTO_TOOLS_H

// MakeDepend: library OMP

#include <omp.h>
#include <memory>
#include <unordered_map>

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "paths/BigMapTools.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/Friends.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "system/ParsedArgs.h"
#include "system/System.h"

// MACROS.

#define REPORT_TIME( CLOCK, MSG ) \
     if (logc.PRINT_TIME_USED) cout << TimeSince(CLOCK) << " " << MSG << endl;

#define REPORT_TIMEX( CLOCK, MSG, OUT )                              \
     if (logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 ) \
          OUT << TimeSince( CLOCK, 1.3 ) << " " << MSG << endl;

#define DATE_MSG(MSG) cout << Date( ) << ": " << MSG << endl;

#define FAIL_MSG(MSG) FatalErr(MSG)

// Class ref_loc: record position of something that's aligned to a reference
// sequence id from start to stop, with rc2 = True if it's reversed relative
// to the reference.

class ref_loc {

     public:

     ref_loc( ) : id(-1) { }
     ref_loc( const int id, const int start, const int stop, const Bool rc2 )
          : id(id), start(start), stop(stop), rc2(rc2) { }

     int id;
     int start;
     int stop;
     Bool rc2;

     Bool Defined( ) const { return id >= 0; }

};

class pairing_status {

     public:

     int partner; // or -1 if unpaired
     int lib;     // or -1 if unpaired

};

class simple_map {

     public:

     vec< pair<String,int> > x;

     int operator[]( const String& s ) const
     {    for ( int j = 0; j < x.isize( ); j++ )
               if ( x[j].first == s ) return x[j].second;
          return 0;    }

};

class long_logging_control {

     public:

// it's too dangerous to have a constructor that leaves a gazillion pointers
// uninitialized
//     long_logging_control( ) { }

     long_logging_control( const ref_data& ref, vec<ref_loc>* readlocs, 
          const String& OUT_INT_HEAD, const String& VERB )
          : G(&ref.G), G3(&ref.G3), G3plus(&ref.G3plus), Gplus_ext(&ref.Gplus_ext), 
          is_circular(&ref.is_circular), ploidy(&ref.ploidy), LG(ref.LG), 
          Glocs(&ref.Glocs), G3locs(&ref.G3locs), G3pluslocs(&ref.G3pluslocs), 
          readlocs(readlocs), OUT_INT_HEAD(OUT_INT_HEAD)
          {    vec<String> m;
               ParseStringSet( VERB, m );
               for ( int j = 0; j < m.isize( ); j++ )
               {    ForceAssert( m[j].Contains( ":" ) );
                    verb.x.push( m[j].Before( ":" ), 
                         m[j].After( ":" ).Int( ) );    }    }

     const vecbasevector* G;
     const vecbasevector* G3;
     const vecbasevector* G3plus;
     const vec<int>* Gplus_ext;
     const vec<bool>* is_circular;
     const vec<double>* ploidy;
     int LG;
     const VecIntPairVec* Glocs;
     const VecIntPairVec* G3locs;
     const VecIntPairVec* G3pluslocs;
     vec<ref_loc>* readlocs;
     String OUT_INT_HEAD;
     simple_map verb;

     vec<placementy> FindGenomicPlacements( const basevector& b ) const
     {    vec<placementy> p = FindGenomicPlacementsY( 0, b, LG, *G3plus, 
               *G3pluslocs );
          vec<Bool> to_delete( p.size( ), False );
          for ( int i = 0; i < p.isize( ); i++ )
          {    p[i].pos -= (*Gplus_ext)[ p[i].g ];
               p[i].Pos -= (*Gplus_ext)[ p[i].g ];
               int n = (*G)[ p[i].g ].size( );
               if ( p[i].pos >= n ) to_delete[i] = True;
               else if ( p[i].Pos > n ) p[i].Pos -= n;    }
          EraseIf( p, to_delete );
          return p;    }

};

// stores reads/quals/pairs and load from disk only if necessarily
class LongProtoReadsQualsPairs{
    //gut of the lazy load mechanism
    template<typename T>
    class path_or_data{
        template<typename U>
        struct has_ReadAll{
            typedef char yes; typedef char (&no)[2];
            template<class V> static yes test( decltype(&V::ReadAll) );
            template<class V> static no test( ... );
            static bool const value = sizeof(test<U>(0))==sizeof(yes);
        };
        template<typename U>
        static typename std::enable_if<has_ReadAll<T>::value,U&>::type
        load_file(U&data_,const String&path_){ data_.ReadAll(path_); return data_; }

        template<typename U>
        static typename std::enable_if<!has_ReadAll<T>::value,U&>::type
        load_file(U&data_,const String&path_){ data_.Read(path_); return data_; }

        template<typename U>
        static typename std::enable_if<has_ReadAll<T>::value,U&>::type
        write_file(U&data_,const String&path_){ data_.WriteAll(path_); return data_; }

        template<typename U>
        static typename std::enable_if<!has_ReadAll<T>::value,U&>::type
        write_file(U&data_,const String&path_){ data_.Write(path_); return data_; }

    public:
        path_or_data():sPath(),ptr(){};
        explicit path_or_data(String const&s)  :sPath(s), ptr(){}
        path_or_data(String const&s, T const&t):sPath(s), ptr(new T(t)){}
        path_or_data(path_or_data&& other):sPath(std::move(other.sPath)),ptr(other.ptr){}
        path_or_data& operator=(const path_or_data&& other){sPath=std::move(other.sPath); ptr=std::move(other.ptr); return *this;}

        const T & get(bool bResetFile=false)const {
            if(bResetFile && sPath.size()>0 && IsRegularFile(sPath)) {
                Remove(sPath);
            }
            if(!ptr){
                ptr.reset( new T() );
                if(sPath.size()>0 && IsRegularFile(sPath)) {
                    load_file(*ptr,sPath);
                }
            }
            return *ptr;
        }
        T& get(bool bResetFile=false){ return const_cast<T&>( static_cast<const path_or_data<T>&>(*this).get(bResetFile) ); }

        bool inMem()const{return static_cast<bool>(ptr);};

        void setPath(const String&s){ sPath=s; }
        void clearPath(){sPath.clear();};
        void clearMem()const{ptr.reset();};
        void write()const{ if(ptr && sPath.size()>0) write_file(*ptr,sPath); };
    private:
        path_or_data(path_or_data&);
        path_or_data& operator=(const path_or_data&);
        String sPath;
        mutable std::unique_ptr<T> ptr;
    };
public:
    explicit LongProtoReadsQualsPairs(const String& head)
        :reads_(head+".fastb") ,quals_(head+".qualb") ,pairs_(head+".pairs") {}
    LongProtoReadsQualsPairs(const String& head,vecbasevector const&a,vecqualvector const&b,PairsManager const&c)
        :reads_(head+".fastb",a) ,quals_(head+".qualb",b) ,pairs_(head+".pairs",c) {}
    vecbasevector const& reads(bool bResetFile=false)const{return reads_.get(bResetFile);}
    vecqualvector const& quals(bool bResetFile=false)const{return quals_.get(bResetFile);}
    PairsManager  const& pairs(bool bResetFile=false)const{return pairs_.get(bResetFile);}

    vecbasevector& reads(bool bResetFile=false){return reads_.get(bResetFile);}
    vecqualvector& quals(bool bResetFile=false){return quals_.get(bResetFile);}
    PairsManager&  pairs(bool bResetFile=false){return pairs_.get(bResetFile);}
    void write(){ reads_.write(); quals_.write(); pairs_.write(); }
    void clearMem()const{
        reads_.clearMem();
        quals_.clearMem();
        pairs_.clearMem();
    }
    bool inMem()const{return reads_.inMem() || quals_.inMem() || pairs_.inMem();};
private:
    path_or_data<vecbasevector> reads_;
    path_or_data<vecqualvector> quals_;
    path_or_data<PairsManager>  pairs_;
};

// deals with Tmp directory's file, intended to replace all "TMP" hard coded stuff
class LongProtoTmpDirManager{
public:
    typedef LongProtoReadsQualsPairs unit_t;
    typedef std::unique_ptr<unit_t> entry_t;
    typedef String key_t;
    struct key_hasher{ size_t operator()(const key_t& in)const{ return std::hash<std::string>()(in); } };
    typedef std::unordered_map<key_t,entry_t,key_hasher> dictionary_t;

    // no explicit here, this is deisgned to be a drop-in replacement of the legacy "String TMP" mechanism
    // passing in a String, instead of a manager instance, creates an unnamed manager which goes out of scope
    // right away and clear the memory
    LongProtoTmpDirManager(const String&base_directory):tmp_dir_(base_directory){};

    // return the unit_t object which gives access to the file tmp_dir_/key.[fastb/qualb/pairs]
    // initially, the object will load such files if the files exist
    const unit_t& operator[](const key_t& key)const{ return this->get(key);};
    unit_t& operator[](const key_t& key){ return this->get(key);};

    const unit_t& get(const key_t& key)const{
        auto itr=dict_.find(key);
        if( itr== dict_.end()){
            itr = dict_.insert( std::make_pair(key ,entry_t( new unit_t(tmp_dir_+"/"+String(key))))).first;
        }
        return *(*itr).second;
    }
    unit_t& get(const key_t&key){ return const_cast<unit_t&>( static_cast<const LongProtoTmpDirManager&>(*this).get(key) ); }
    const String& dir()const{return tmp_dir_;};

    void clearMem(const key_t&key)const{
        auto itr=dict_.find(key);
        if( itr != dict_.end() ){ (*itr).second->clearMem(); }
    }

private:
    LongProtoTmpDirManager(const LongProtoTmpDirManager&);
    LongProtoTmpDirManager(LongProtoTmpDirManager&&);
    LongProtoTmpDirManager& operator=(const LongProtoTmpDirManager&);
    LongProtoTmpDirManager& operator=(LongProtoTmpDirManager&&);
    mutable dictionary_t dict_;
    const String tmp_dir_;
};

template<int K> void BuildCorrectedReads( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid,
     VecEFasta& corrected, vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS );

void ReportPeakMem( const String msg = "" );

void Done( const double clock ) __attribute__((__noreturn__));

int SelectK2( const VecEFasta& corrected, const double K2frac,
     const long_logging& logc, const long_heuristics& heur );

void ReportExcessEdges( const HyperEfasta& he, const vecbasevector& G,
     const Bool PERF_STATS );

void DumpSimLocs( const vecbasevector& reads, const vec<ref_loc>& readlocs );

int StartStep( const String& IN_SHBV, const int START_STEP );

void TestSampleAndReads( const String& SAMPLE, const String& READS,
     long_logging& logc );

void PrintPerformanceStats( const double clock, const long_logging& logc );

void PrintCorrectedReadStats( const VecEFasta& corrected );

void PrintTail( const parsed_args& command, const String& X, const String& X_actual,
     const double clock, const long_logging& logc );

void LoadEfastaLong( const String& fn, VecEFasta& x,
     vec<int>& cid, vec<pairing_info>& cpartner, const long_logging& logc );

#endif
