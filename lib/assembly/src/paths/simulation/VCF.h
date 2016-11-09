///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//
// http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
//


#ifndef VCF_H_
#define VCF_H_

#include <string>
#include <vector>
#include <limits>
#include "Basevector.h"

using std::string;
using std::vector;


// VCFWriter use a streamlined data structure compared to VCF, which was hard coded to use VCFChromosome
class VCFWriter
{
public:
    enum{ VCF_CHROM, VCF_POS, VCF_ID, VCF_REF, VCF_ALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES};
    struct line_t{
        line_t(const String&chrom,size_t pos1):CHROM(chrom),POS1(pos1){};
        line_t():POS1(std::numeric_limits<size_t>::max()){};

        line_t(const String&chrom,size_t pos1,const String&id,const String&ref,const vec<String>&alts):CHROM(chrom),POS1(pos1),ID(id),REF(ref),ALTs(alts){};
        line_t(const String&chrom,size_t pos1,const String&id,const String&ref,const vec<String>&alts,const String&qual
              ,const vec<size_t>&filters
              ,const vec<std::pair<size_t,String>>&infos
              ,const vec<vec<String> >& gtfs):CHROM(chrom),POS1(pos1),ID(id),REF(ref),ALTs(alts),QUAL(qual),FILTERs(filters),INFOs(infos),GTFs(gtfs){};

        String CHROM;
        size_t POS1;
        String ID;
        String REF;
        vec<String> ALTs;
        String QUAL;
        vec<size_t> FILTERs;
        vec<std::pair<size_t,String>> INFOs;
        vec< vec<String> > GTFs;
        bool operator<(const VCFWriter::line_t&right)const {
            return (CHROM < right.CHROM) || ( CHROM==right.CHROM && POS1 < right.POS1);
        }
    };
    void sort(){ Sort(entries); }

    struct info_t{
        info_t(const String&a):ID(a){};
        info_t(const String&a,const String&b,const String&c,const String&d):ID(a),Number(b),Type(c),Description(d){};
        void set(const String&a,const String&b,const String&c,const String&d){ID=a;Number=b;Type=c;Description=d;};
        String ID;
        String Number;
        String Type;
        String Description;
    };

    VCFWriter();
    void AddEntry(const String&chrom,size_t pos1,const String&id,const String&ref,const vec<String>&alts){
        entries.push_back(line_t(chrom,pos1,id,ref,alts));
    };
    void AddEntry(const String&chrom                          //chromosome name
                 ,size_t pos1                                 //1-based position
                 ,const String&id                             //unique variant identifier
                 ,const String&ref                            //reference
                 ,const vec<String>&alts                      // list of alts
                 ,const String&qual                           // quality score, can be .
                 ,const vec<String>&filters                   // list of failed filters
                 ,const vec<std::pair<String,String>>&infos   // list of additional infos
                 ,const std::map<String,vec<std::pair<String,String>>> & gtfs); //sample, and a list of format-value
    void setSample(const String&name,const String&description);
    void setFormat(const String&name,const String& number, const String& type, const String&description);

    void AddMeta( const String& tag, const String& value ) {
	ostringstream s;
	s << "##" << tag << "=" << value;
	_svMETA.push_back( s.str() );
    }

    void AddMeta( const String& tag, std::vector< std::pair< const String, const String > > kv_pairs ) {
	ostringstream s;
	s << "##" << tag << "=<";

	for ( auto const& pair: kv_pairs ) {
	    s << pair.first << "=" << pair.second;
	    if ( &pair != &kv_pairs.back() ) s << ",";
	}

	s << ">";
	_svMETA.push_back(s.str());
    }

    void AddInfo( const info_t& info ) {
	_INFOs.push_back(info);
    }

private:
    vec<String> _svMETA; // line-by-line storage of the file META
    vec<String> _svHeader; //term-by-term of the data header, should be >8 in VCF v4.1
    vec<info_t> _SAMPLEs;

    vec<info_t> _FILTERs; //list of FILTER
    vec<info_t> _INFOs; //list of infos
    vec<info_t> _FORMATs; //list of genotype formats
    vec<line_t> entries;

    friend ostream& operator<<(ostream&,const VCFWriter&);

};
ostream& operator<<(ostream&,const VCFWriter&);


// class storing line-by-line data entries of a particular chromosome, behaves like a std::vector
class VCFChromosome{
public:
  //value of each line
  struct VCFLine {
    struct alt_t{ // in VCF spec, the alt strings
      alt_t(const std::string& sIn);
      std::string sRep;
      BaseVec     bRep;
      bool isBaseVec()const{return sRep.empty();}
    };
    size_t POS1;
    string ID;
    BaseVec REF;
    vector<alt_t> ALT;
    double QUAL;
    string sQUAL;
    string FILTER;
    string INFO;
    string GENOTYPE_FORMAT;				// e.g. GT:AD:DP:GQ:PL:TP
    vector<string> 		 GENOTYPE;		// e.g. 0|0:139,2:141:15:0,15,196:15
							// OR just ./., even if FORMAT has other fields
    vector<map<string, string> > GENOTYPE_FIELDS;	// e.g. GT->0|0, AD->139,2, DP->141, GQ->15, PL->0,15,196, TP->15
    vector<string> 		 GT;
    vector< vector< char> > 	 GT_IDX;
    vector<bool> 		 GT_PHASED;

//    string gt_format;
//    vector<string> gt;
    VCFLine() : POS1(0), QUAL(NAN) {};
    
    void Print( std::ostream& os ) const;
    void Write( const std::string& sChrom, std::ostream& os ) const;
  };

  // constructors: must provide a name of the chromosome
  VCFChromosome(const std::string& name):_sName(name){};
  VCFChromosome(const char* const cp):_sName(cp){};

  //name of the chromosome
  const std::string& name()  const{return _sName;};
  size_t size()const{return _lines.size();};

  //iterators for going through line-by-line data
  typedef typename std::vector<VCFLine>::iterator iterator;
  typedef typename std::vector<VCFLine>::const_iterator const_iterator;
  iterator begin(){return _lines.begin();};
  iterator end(){return _lines.end();};
  const_iterator begin()const{return _lines.begin();};
  const_iterator end()const{return _lines.end();};
  //c++11
  const_iterator cbegin() const noexcept{return _lines.cbegin();};
  const_iterator cend() const noexcept{return _lines.cend();};

  size_t idx( iterator i ) { return i - this->begin(); };
  size_t idx( const_iterator i ) const { return i - this->cbegin(); };

  //reference of the last line
  VCFLine& back(){return _lines.back();}

  //add an entry to the end
  void push_back(const VCFLine& val){_lines.push_back(val);};

private:
  vector<VCFLine> _lines;
  std::string _sName;
};

class VCF {
public:
  // VCF constructor   - read full VCF from a .vcf file.
  //          filename - input file name
  // sChromosome       - if != "", read only the data for tag==sChromosome
  // lBegin            - if sChromosome!="", read data line starting from 1<=lBegin for that chromosome
  // lEnd              - if sChromosome!="" and lEnd>lBegin>=1, read data lines indexed [lBegin,lEnd) for that chromosome
  // stDump	       - if sDump != "", then it's a filename to dump raw lines of header and matching locs to
  // (Obsolete) sColumnDelimiters - by default, VCF spec use \t, but it can be a list like " \t" here
  VCF( const std::string& filename
     , const std::string& sChromosome=""
     , size_t lBegin1=1
     , size_t lEnd1=numeric_limits<size_t>::max()
     , const string& sDump = "");
//     , const std::string& sColumnDelimiters="\t");

  //iterators for going through the chromosomes
  typedef typename std::vector<VCFChromosome>::iterator iterator;
  typedef typename std::vector<VCFChromosome>::const_iterator const_iterator;
  iterator begin() {return _vChromosomes.begin();};
  iterator end() {return _vChromosomes.end();};
  const_iterator begin() const {return _vChromosomes.begin(); };
  const_iterator end() const {return _vChromosomes.end(); };
  //c++11
  const_iterator cbegin() const noexcept{return _vChromosomes.cbegin();};
  const_iterator cend() const noexcept{return _vChromosomes.cend();};

  //number of chromosomes
  size_t size()const{return _vChromosomes.size();};

  //number of samples, i.e. number of tags after the FORMAT tag
  size_t getNSamples()const { return ((VCF_SAMPLES>_svHeader.size())?(0):(_svHeader.size()-VCF_SAMPLES)); };

  // get line-by-lie file meta
  const vector<std::string>& getMETA()const{return _svMETA;};

  // get tag-by-tag data header
  const vector<std::string>& getHeader()const{return _svHeader;};

  // get index of a sample with a given name, or -1 on error
  int getSampleIndex( const std::string& name ) {
      auto iter = std::find( _svSamples.begin(), _svSamples.end(), name ) ;

      if ( iter != _svSamples.end() )
	  return iter - _svSamples.begin();
      else
	  return -1;
  }

  void writeSimplifiedVCF(const std::string&filename);

private:
  vector<std::string> _svMETA; // line-by-line storage of the file META

  vector<std::string> _svHeader; //term-by-term of the data header, should be >8 in VCF v4.2

  vector<std::string> _svSamples; //sample headers, after the FORMAT tag

  vector< VCFChromosome > _vChromosomes; //list of chromosomes

  enum{ VCF_CHROM, VCF_POS, VCF_ID, VCF_REF, VCF_ALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES};

  //quick tokenizer which returns a std::pair such that str.substr(pair.first,pair.second) is the substring
  //starting at str[uStart] ending before the delimiters. uStart is then modified to be the first
  //position after delimiters, note that std::string:npos tricks are implied.
  std::pair<size_t,size_t> nextToken(const std::string&str,const std::string&delimiters,size_t&uStart);

  static void dumpLine( const string& filename, ofstream& stDump, const string& line ) {
      if ( stDump.is_open() ) {
	  stDump << line << endl;
	  if ( ! stDump.good() ) FatalErr("error writing VCF output to " << filename );
      }
  }
};


#endif /* VCF_H_ */
