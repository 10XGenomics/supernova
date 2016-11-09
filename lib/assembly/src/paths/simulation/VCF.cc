///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//

#include <cstdio>
#include "feudal/FeudalString.h"
#include <iostream>
#include <regex>
#include "paths/simulation/VCF.h"
#include "system/Assert.h"
#include <utility>

namespace {

// split - tokenize a string based on a separator and return a vector of strings representing the
// characters between each instance of a separator.  If there is no separator, then the vector will
// have one entry equal to the original string.
vector<string> split( const string& input, const char sep=':' )
{
    vector<string> tokens;

    for (auto i = input.begin(); i != input.end(); /* */ ) {
	tokens.push_back(string());

	while ( i != input.end() && *i != sep )
	    tokens.back().push_back(*i++);

	if ( i != input.end() ) i++;
    }

    return tokens;
}


// parse_genotype_fields -- take current line's GENOTYPE and break out into
// GENOTYPE_FIELDS dictionary.  GENOTYPE_FIELDS[i]["GT"] then gets broken out
// into GT_IDX[i] and GT_PHASED[i].
//
// e.g.:
// 	GT 0/1 0|2 indicating genotypes for two samples becomes:
//
// GENOTYPE_FIELDS[0]["GT"] = "0/1"		// first sample (index 0)
// GT_IDX[0][0] = 0, GT_IDX[0][1] = 1, GT_PHASED=false
//
// GENOTYPE_FIELDS[1]["GT"] = "0|2"		// second sample (index 1)
// GT_IDX[0][0] = 0, GT_IDX[0][1] = 2, GT_PHASED=true
//

void parse_genotype_fields( VCFChromosome::VCFLine& currentLine )
{
    /* NEED to parse GENOTYPE tags into GENOTYPE_FIELDS map, rather than
     * just the GT standalone.  So we'll do away with GT and just use
     * GENOTYPE_FIELDS[idx]["GT"] if we need it.  But GT_IDX and GT_PHASED should stay
     *
     */
    size_t nSamples = currentLine.GENOTYPE.size();
    ForceAssertEq( currentLine.GENOTYPE_FORMAT.find("GT"), 0U );	// first should be GT

    vector<string> field_names = split( currentLine.GENOTYPE_FORMAT );	// tokenize based on ':'
    ForceAssertGe(field_names.size(), 1U);


    //---- here we parse out all of the two-letter tags (e.g. "GT")
    for ( size_t i = 0; i < nSamples; ++i ) {		// for each sample
	// we should have either a single faux GT entry ./. or the entire thing from GENOTYPE_FORMAT
	vector<string> field_vals = split( currentLine.GENOTYPE[i] );

	currentLine.GENOTYPE_FIELDS.push_back(map<string,string>());

	if ( field_vals.size() == 1 && field_names.size() > 1 ) {	// case where only ./. is shown
	    ForceAssertEq(field_vals.back(), "./.");
	    for ( vector<string>::iterator j = field_names.begin(); j != field_names.end(); ++j ) {
		if ( j == field_names.begin() )
		    currentLine.GENOTYPE_FIELDS.back()[*j] = *field_vals.begin();
		else
		    currentLine.GENOTYPE_FIELDS.back()[*j] = "";
	    }
	} else { 							// normal case -- map vals to fields
	    ForceAssertEq( field_names.size(), field_vals.size() );
	    auto k = field_vals.begin();
	    for ( auto j = field_names.begin();
		    j != field_names.end() && k != field_vals.end();
		    ++j, ++k ) {
		currentLine.GENOTYPE_FIELDS.back()[*j] = *k;
	    }
	}
    }


    //---- "GT" tag processing -- gets translated out into GT_IDX and GT_PHASED
    ForceAssertEq(nSamples, currentLine.GENOTYPE_FIELDS.size() );
    currentLine.GT_IDX.resize( nSamples );
    currentLine.GT_PHASED.resize( nSamples, false );

    for ( size_t i = 0; i < nSamples; ++i ) {
	// grab genotype into string buffer
	istringstream s(currentLine.GENOTYPE_FIELDS[i]["GT"]);
	currentLine.GT_IDX[i].clear(); 	// clear index
	int idx;			// we use int here because we want a numeric conversion

	// idx gets a -1 if we have a missing field marker '.'
	// otherwise, the index in the field 0,1,...
	if ( s.peek() != '.' ) {
	    s >> idx;
	    ForceAssertLe(idx, numeric_limits<char>::max());
	}
	else {
	    s.ignore();
	    idx = -1;
	}
	currentLine.GT_IDX[i].push_back( idx );

	// still more left (i.e. diploid?)
	int ph;
	if ( ( ph = s.get() ) > 0 && !s.fail()  ) {
	    // grab the phasing character | or /
	    if (ph == '|')
		currentLine.GT_PHASED[i] = true;
	    else if ( ph != '/' ) {
		    cout << "about to die: currentLine sample i="  << i << " GT=[" << currentLine.GENOTYPE_FIELDS[i]["GT"] << "] genotype=[" << currentLine.GENOTYPE[i] << "]"<< endl;
		    FatalErr("received an unexpected phasing character in the genotype: int=" << ph << ", char=" << static_cast<char>(ph) );
		}

	    if ( s.peek() != '.' ) {
		s >> idx;
		ForceAssertLe(idx, numeric_limits<char>::max());
	    }
	    else
		idx = -1;
	    currentLine.GT_IDX[i].push_back( idx );
	}
    }
}

};
void
VCFChromosome::VCFLine::Write( const std::string&sChrom, ostream& os ) const
{
    os << sChrom << '\t' 
       << POS1 << '\t' 
       << ID << '\t' 
       << REF.ToString() << '\t';
    for ( size_t i = 0; i < ALT.size(); ++i ){
        if (i>0 ) os << ",";
        if ( ALT[i].isBaseVec() )
            os << ALT[i].bRep.ToString();
        else
            os << ALT[i].sRep;
    }
    os <<'\t' 
       << sQUAL << '\t'
       << FILTER << '\t'
       << INFO << '\t'
       << GENOTYPE_FORMAT;

    for(const auto& entry : GENOTYPE){
        os << '\t' << entry;
    }
    os << '\n';
}

void
VCFChromosome::VCFLine::Print( ostream& os ) const
{
    os << "POS1 " << POS1 << endl;
    os << "ID   " << ID << endl;
    os << "REF  " << REF.ToString() << endl;
    os << "ALT (" << ALT.size() << ")" << endl;
    for ( size_t i = 0; i < ALT.size(); ++i )
	if ( ALT[i].isBaseVec() )
	    cout << "\tALT " << i << " " << ALT[i].bRep.ToString() << endl;
	else
	    cout << "\tALT " << i << " " << ALT[i].sRep << endl;
    os << "FILTER " << FILTER << endl;
    os << "INFO   " << INFO << endl;
    os << "GENOTYPE_FORMAT " << GENOTYPE_FORMAT << endl;
    os << "GENOTYPE_FIELDS.size() " << GENOTYPE_FIELDS.size() << endl;
    for ( size_t i = 0; i < GENOTYPE_FIELDS.size(); ++i ) {

	for ( auto j = GENOTYPE_FIELDS[i].begin(); j != GENOTYPE_FIELDS[i].end(); ++j )
	    cout << "\tGENOTYPE_FIELDS " << i << " [" << (j->first) << "] "  << (j->second) << endl;

	cout << "\tGT " << i << " " << GT[i] << endl;
	for ( size_t j = 0; j < GT_IDX[i].size(); ++j )
	    cout << "\t\tGT_IDX " << i << " = " << static_cast<int>(GT_IDX[i][j]) << endl;
	cout << "\t\tGT_PHASED " << ( GT_PHASED[i] ? "true" : "false") << endl;

    }
}


VCFChromosome::VCFLine::alt_t::alt_t(const string& sIn)
  :sRep(sIn){
  if(sIn.find_first_of("<>.")==std::string::npos){
    //there could be space characters in some fields, and BaseVec doesn't like it
    sRep.erase( remove( sRep.begin(), sRep.end(), ' ' ), sRep.end() );
    bRep = BaseVec(sRep);
    sRep.clear();
  }
}

void VCF::writeSimplifiedVCF(const std::string&filename){
  std::ofstream outfile(filename, std::ios::out );
  for(const auto& entry: _svMETA){ outfile << entry << '\n'; }
  for(const auto& chrom: _vChromosomes){
      for(const auto& line:chrom){
          bool bSuccess=false;
          if( line.ALT.size()==1){
              const BaseVec& ref=line.REF;
              bool bFixable=true;
              if( !line.ALT[0].isBaseVec() ){
                  for(const auto& entry: line.ALT[0].sRep){
                      if(  entry!='a' && entry!='t' && entry!='g' && entry!='c'
                         &&entry!='A' && entry!='T' && entry!='G' && entry!='C'){
                        bFixable=false;
                        break;
                      }
                  }
              }
              BaseVec alt=( line.ALT[0].isBaseVec() || !bFixable)? (line.ALT[0].bRep):(BaseVec(String(line.ALT[0].sRep)));
              if (ref.size()>1 && alt.size()>1 && bFixable){
                  vec<int64_t> ref_to(ref.size()+1,-1);
                  ref_to[ref.size()]=alt.size();
                  if(ref.size()==alt.size()){
                      for(size_t ii=0;ii<ref.size();++ii){ ref_to[ii]=ii; }
                      bSuccess=true;
                  }
                  else if( ref.size()<alt.size()){
                      size_t ref_idx=0;
                      for( ; ref_idx<ref.size() && ref[ref_idx]==alt[ref_idx]; ++ref_idx){
                          ref_to[ref_idx]=ref_idx;
                      }
                      bSuccess=ref_idx>0;
                      if(bSuccess){
                          int64_t alt_lower=ref_idx;
                          for( ; ref_idx<ref.size() ;++ref_idx){
                              const int64_t upper=alt.size()+ref_idx-ref.size();
                              ref_to[ref_idx]=upper;
                              for( int64_t alt_idx=upper-1
                                 ; alt_idx >= alt_lower 
                                 ; --alt_idx){
                                  if( ref[ref_idx] == alt[alt_idx] ) {
                                      ref_to[ref_idx] = alt_idx;
                                  }
                              }
                              alt_lower=ref_to[ref_idx]+1;
                          }
                      }
                  }
                  else{ //alt is shorter
                      int ref_lower=0;
                      for( int64_t alt_idx=0 ; alt_idx<alt.size() ; ++alt_idx){
                          const int64_t upper= ref.size() + alt_idx - alt.size();
                          int64_t earliest_ref_idx=upper;
                          for( int64_t ref_idx=upper-1; ref_idx>=ref_lower;--ref_idx){
                              if( ref[ref_idx] == alt[alt_idx] ){
                                  earliest_ref_idx=ref_idx;
                              }
                          }
                          ref_to[earliest_ref_idx]=alt_idx;
                          ref_lower=earliest_ref_idx+1;
                      }
                      bSuccess=true;
                  }
                  if(bSuccess){
                      for( size_t ref_idx=0;ref_idx<ref.size();){
                          if( ref_to[ref_idx] < 0) continue;
                          auto new_line=line;
                          size_t ref_end=ref_idx+1;
                          for(;ref_end<ref.size() && ref_to[ref_end]<0 ; ++ref_end){}
                          new_line.REF.SetToSubOf(ref, ref_idx, ref_end-ref_idx);
                          new_line.ALT[0].sRep.clear();
                          new_line.ALT[0].bRep.SetToSubOf(alt, ref_to[ref_idx], ref_to[ref_end]-ref_to[ref_idx]);
                          new_line.POS1=line.POS1+ref_idx;
                          new_line.Write(chrom.name(), outfile);
                          ref_idx=ref_end;
                      }
                  }
              }
          }
          if(!bSuccess) line.Write(chrom.name(), outfile);
      }
  }
  outfile.close();
}

VCF::VCF( const std::string& filename
        , const std::string& sChromosome
        , const size_t lBegin1
        , const size_t lEnd1
        , const string& sDump /* = "" */
//        , const std::string& sColumnDelimiters
        )

{
  const std::string& sColumnDelimiters="\t"; // VCF spec said using tab, so let's use tabs

  ForceAssert( sChromosome == "" || lEnd1 > lBegin1 );

  //prepare to read file line-by-line
  std::ifstream infile(filename, std::ios::in );
  ForceAssert( infile.is_open() );

  string sBuffer; //for readline

  // output dump
  ofstream outfile;
  if ( sDump != "" ) OpenOfstream(outfile, sDump, sDump, ios::trunc);

  //log all the meta file info before the data header line
  bool bHeaderFound = false;
  while( std::getline(infile, sBuffer) ){
    _svMETA.push_back(sBuffer);
    dumpLine( sDump, outfile, sBuffer );
    if( std::regex_match(sBuffer,std::regex("#CHROM.*"))){
      bHeaderFound=true;
      break;
    }
  }
  ForceAssert(bHeaderFound);
  // std::cout << "Number of lines of file meta info:   " <<_svMETA.size() << std::endl;

  {//special treatment of data header -- make it "space" tolerent
    std::string sLocDelimiter="\t";

    size_t nTabs=std::count(sBuffer.begin(),sBuffer.end(),'\t');
    if(nTabs>0){
        if(nTabs < VCF_SAMPLES){
            std::cerr << "Data header is tab-delimited but does not have enough fields"<<std::endl;
            ForceAssert ( nTabs>=VCF_SAMPLES);
        }
    }
    else{
        std::cerr << "WARNING: Data header has no tabs, trying to delimit by space instead"<<std::endl;
        sLocDelimiter=" ";
    }

    //parse the data header the rest column-by-column
    for( size_t token_start = sBuffer.find_first_not_of(sLocDelimiter,0)
       ; token_start < sBuffer.size()
       ;
       ){
      std::pair<size_t,size_t> substr_params = nextToken(sBuffer,sLocDelimiter,token_start);
      if( substr_params.first>sBuffer.size() ) break;
      _svHeader.push_back(sBuffer.substr(substr_params.first,substr_params.second));
    }
  }


  _svSamples.assign( _svHeader.begin()+VCF_SAMPLES, _svHeader.end());

  /*
  std::cout << "Number of fields in the data header: " <<_svHeader.size() << std::endl;
  for(size_t ii=0;ii<_svHeader.size();++ii){
    std::cout << _svHeader[ii] << " ";
  }
  std::cout << std::endl;
  std::cout << "Number of samples in the data header: " <<_svSamples.size() << std::endl;
  for(size_t ii=0;ii<_svSamples.size();++ii){
    std::cout << _svSamples[ii] << " ";
  }
  std::cout << std::endl;
  */

  //now parse the rest line-by-line
  bool bFoundSelection=false;
  const bool bSelectRange = (sChromosome!="") && (lBegin1>0) && (lEnd1>lBegin1);
  // cout << "bSelectRange=" << (bSelectRange ? "true":"false") << ", lBegin1=" << lBegin1 << ", lEnd1=" << lEnd1 << endl;
  std::string sTmp;
  size_t nDataLine=0;

  // if we're matching on chromosome, push back an empty so that we have something,
  // even if we match no lines in the file.
  if ( sChromosome != "" )
      _vChromosomes.push_back( VCFChromosome( sChromosome ) );

  while( std::getline(infile, sBuffer) ){
    if( sBuffer.size() == 0) continue;

    size_t token_start = sBuffer.find_first_not_of(sColumnDelimiters,0);

    unsigned int column = 0;

    // Parse the CHROM tag
    ForceAssert ( token_start < sBuffer.size() );
    std::pair<size_t,size_t> substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    sTmp.assign(sBuffer,substr_params.first,substr_params.second);

    if(sChromosome!="" && sChromosome!=sTmp){
      if(!bFoundSelection) continue;
      else                 break;
    }
    else{
      bFoundSelection=true;
    }

    if( _vChromosomes.size()==0 || sTmp!= _vChromosomes.back().name()){
      _vChromosomes.push_back( VCFChromosome(sTmp));
      std::cout << "Reading chromosome " << sTmp << std::endl;
    }
    ++column;

    // Parse POS tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    size_t pos1 = std::atoll(sBuffer.substr(substr_params.first,substr_params.second).c_str());

    if ( bSelectRange && ( pos1 < lBegin1 || pos1 >= lEnd1 ) )
	 continue;

    dumpLine( sDump, outfile, sBuffer );

    // now we're committed to the line
    _vChromosomes.back().push_back(VCFChromosome::VCFLine());
    VCFChromosome::VCFLine & currentDataLine = _vChromosomes.back().back();

    // Store the POS tag
    currentDataLine.POS1 = pos1;
    ++column;

    // Parse ID tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    currentDataLine.ID.assign(sBuffer,substr_params.first,substr_params.second);
    ++column;

    // Parse REF tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    currentDataLine.REF = BaseVec(std::string(sBuffer,substr_params.first,substr_params.second));
    ++column;

    // Parse ALT tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    {
      std::string sBufferLoc(sBuffer,substr_params.first,substr_params.second);
      size_t start = 0;
      for( size_t found=sBufferLoc.find_first_of(",",start)
         ; found != std::string::npos
         ; found = sBufferLoc.find_first_of(",",start)
         ){
        currentDataLine.ALT.push_back(sBufferLoc.substr(start,found-start));
        start=found+1;
      }
      currentDataLine.ALT.push_back(sBufferLoc.substr(start,sBufferLoc.size()-start));
    }
    ++column;

    // Parse QUAL tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    if( sBuffer.substr(substr_params.first,substr_params.second) == ".") {
        currentDataLine.QUAL = -std::numeric_limits<double>::epsilon();
    }
    else{
        currentDataLine.QUAL = std::atof(sBuffer.substr(substr_params.first,substr_params.second).c_str());
    }
    currentDataLine.sQUAL = sBuffer.substr(substr_params.first,substr_params.second);
    ++column;

    // Parse FILTER tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    currentDataLine.FILTER.assign(sBuffer,substr_params.first,substr_params.second);
    ++column;

    // Parse INFO tag
    ForceAssert ( token_start < sBuffer.size() );
    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
    ForceAssert ( substr_params.first < sBuffer.size() );
    currentDataLine.INFO.assign(sBuffer,substr_params.first,substr_params.second);
    ++column;

    if( column < _svHeader.size() ){
        // Parse FORMAT tag, warn about  missing data
        if(token_start<sBuffer.size()){
            substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
            if ( substr_params.first < sBuffer.size() ){
                currentDataLine.GENOTYPE_FORMAT.assign(sBuffer,substr_params.first,substr_params.second);
                ++column;

                // Pasrse GENOTYPE_FORMAT entries, warn about missing data
                for( ; column < _svHeader.size() ; ++column){
                    if(token_start>=sBuffer.size()){
                        std::cerr << "WARNING: Missing entry for GENOTYPE" << std::endl;
                        break;
                    }
                    substr_params = nextToken(sBuffer,sColumnDelimiters,token_start);
                    if ( substr_params.first >= sBuffer.size() ){
                        std::cerr << "WARNING: Unable to split GENOTYPE entry" << std::endl;
                        break;
                    }
                    currentDataLine.GENOTYPE.push_back( sBuffer.substr(substr_params.first,substr_params.second));
                    currentDataLine.GT.push_back( currentDataLine.GENOTYPE.back().substr(0
                                                                                        ,currentDataLine.GENOTYPE.back().find_first_of(':')
                                                                                        )
                                                );
                }//for
               parse_genotype_fields( currentDataLine );
            } //if ( substr_params.first < sBuffer.size() ){
            else{
                std::cerr << "WARNING: Unable to split FORMAT entry" << std::endl;
            }
        }// if(token_start<sBuffer.size()){
        else{
            std::cerr << "WARNING: Missing entry for FORMAT" << std::endl;
        }
    }//if( column < _svHeader.size() ){
    ++nDataLine;
  } //while( std::getline(infile, sBuffer) )
  if ( sChromosome != "" && !bFoundSelection ) {
      std::ostringstream s;
      s << "You specified chromosome " << sChromosome <<
          ", which was not seen in the VCF file." << std::endl <<
           "If your genome specifies chromosome names like 1, 2, 3..., but "
           << std::endl <<
          "your VCF uses chr1, chr2, chr3..., then you can set VCF_CHR_PREFIX=chr" << endl;
      FatalErr( s.str() );
  }
  // std::cout << "Full processed " << nDataLine << " data lines."<<std::endl;


  infile.close();
  if ( sDump != "" ) outfile.close();

}

std::pair<size_t,size_t> VCF::nextToken(const std::string& str
                                       ,const std::string& delimiters
                                       ,size_t& uStart
                                       ){
  size_t delimiter_position = str.find_first_of(delimiters,uStart);
  std::pair<size_t,size_t> out(uStart,delimiter_position-uStart);
  uStart = str.find_first_not_of(delimiters,delimiter_position);
  return out;
};

VCFWriter::VCFWriter():_svMETA({"##fileformat=VCFv4.1"})
                  ,_svHeader({"CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"})
                  ,_FORMATs({info_t("GT","1","String","Genotype")})
{};
ostream& operator<<(ostream&os,const VCFWriter&in){
    for(const auto&entry:in._svMETA){ os<<entry<<std::endl; }
    for(const auto&entry:in._INFOs){ os<<"##INFO=<ID=" << entry.ID << ",Number="<< entry.Number<<",Type="<<entry.Type<<",Description=\""<<entry.Description<<"\">"<<std::endl; }
    for(const auto&entry:in._FILTERs){ os<<"##FILTER=<ID=" << entry.ID << ",Number="<< entry.Number<<",Type="<<entry.Type<<",Description=\""<<entry.Description<<"\">"<<std::endl; }
    for(const auto&entry:in._FORMATs){ os<<"##FORMAT=<ID=" << entry.ID << ",Number="<< entry.Number<<",Type="<<entry.Type<<",Description=\""<<entry.Description<<"\">"<<std::endl; }
    for(const auto&entry:in._SAMPLEs){ os<<"##SAMPLE=<ID=" << entry.ID << ",Description=\""<<entry.Description<<"\">"<<std::endl; }
    os<<"#";
    for(size_t ii=0;ii<in._svHeader.size();++ii){
        if(ii>0) os<<"\t";
        os<<in._svHeader[ii];
    }
    if(in._FORMATs.size()>0 || in._SAMPLEs.size()>0) os << "\tFORMAT";
    for(const auto&entry:in._SAMPLEs){ os<<"\t"<<entry.ID; }
    os<<std::endl;
    for(const auto&entry: in.entries){
        os << ((entry.CHROM.size()==0)?("."):(entry.CHROM)) << "\t"
           << entry.POS1 << "\t"
           << ((entry.ID.size()==0)?("."):(entry.ID)) << "\t"
           << ((entry.REF.size()==0)?("."):(entry.REF)) << "\t";
        if( entry.ALTs.size()==0){
            os<< ".";
        }
        else{
            os << ((entry.ALTs[0].size()==0)?("."):(entry.ALTs[0]));
            for(size_t ii=1;ii<entry.ALTs.size();++ii){
                os << "," << ((entry.ALTs[ii].size()==0)?("."):(entry.ALTs[ii]));
            }
        }
        os <<"\t";
        os << ((entry.QUAL.size()==0)?("."):(entry.QUAL)) << "\t";
        if( in._FILTERs.size()==0){
            os<< ".";
        }
        else{
            if( entry.FILTERs.size()==0){
                os<< "PASS";
            }
            else{
                os << in._FILTERs[entry.FILTERs[0]].ID;
                for(size_t ii=1;ii<entry.FILTERs.size();++ii){
                    os << "," << in._FILTERs[entry.FILTERs[ii]].ID;
                }
            }
        }
        os <<"\t";
        if( entry.INFOs.size()==0){
            os<< ".";
        }
        else{
            os << in._INFOs[entry.INFOs[0].first].ID;
            if( entry.INFOs[0].second.size()>0){
                os << "=" << entry.INFOs[0].second;
            }
            for(size_t ii=1;ii<entry.INFOs.size();++ii){
                os << in._INFOs[entry.INFOs[ii].first].ID;
                if( entry.INFOs[ii].second.size()>0){
                    os << "=" << entry.INFOs[ii].second;
                }
            }
        }
        if(in._FORMATs.size()>0 || in._SAMPLEs.size()>0){
            os <<"\t";
            if( in._FORMATs.size() == 0){
                os << "."; // this shouldn't happen
            }
            else{
                os << in._FORMATs[0].ID;
                for(size_t ii=1;ii<in._FORMATs.size();++ii){
                    os << ":" << in._FORMATs[ii].ID;
                }
            }
            for(size_t ss=0;ss<in._SAMPLEs.size();++ss){
                os <<"\t";
                if(ss < entry.GTFs.size()){
                    const auto& fields=entry.GTFs[ss];
                    for(size_t ff=0;ff<in._FORMATs.size();++ff){
                        if(ff>0) os<<":";
                        if( ff < fields.size()){
                            os<<fields[ff];
                        }
                        else{
                            os<<".";
                        }
                    }
                }
                else{
                    for(size_t ff=0;ff<in._FORMATs.size();++ff){
                        if(ff>0) os<<":";
                        os<<".";
                    }
                }
            }
        }
        os <<"\n";
    }
    return os;
}
void VCFWriter::AddEntry(const String&chrom,size_t pos1,const String&id,const String&ref,const vec<String>&alts,const String&qual
                 ,const vec<String>&filters
                 ,const vec<std::pair<String,String>>&infos
                 ,const std::map<String,vec<std::pair<String,String>>> & gtfs){
    vec<size_t> filter_indices;
    for(const auto&entry:filters){
        size_t ii=0;
        for(ii=0; ii<_FILTERs.size()&&_FILTERs[ii].ID!=entry;++ii){ }
        if( ii == _FILTERs.size()){
            _FILTERs.push_back(info_t(entry));
        }
        filter_indices.push_back(ii);
    }
    vec<std::pair<size_t,String>> indexed_infos;
    for(const auto&entry:infos){
        size_t ii=0;
        for(ii=0; ii<_INFOs.size()&&_INFOs[ii].ID!=entry.first;++ii){ }
        if( ii == _INFOs.size()){
            _INFOs.push_back(info_t(entry.first));
        }
        indexed_infos.push_back(std::make_pair(ii,entry.second));
    }
    for(const auto&entry:gtfs){
        size_t ii=0;
        for(ii=0; ii<_SAMPLEs.size()&&_SAMPLEs[ii].ID!=entry.first;++ii){ }
        if( ii == _SAMPLEs.size()){
            _SAMPLEs.push_back(info_t(entry.first));
        }
        for(const auto&field_entry:entry.second){
            size_t jj=0;
            for(jj=0; jj<_FORMATs.size()&&_FORMATs[jj].ID!=field_entry.first;++jj){ }
            if( jj == _FORMATs.size()){
                _FORMATs.push_back(info_t(field_entry.first));
            }
        }
    }
    vec<vec<String>> samples(_SAMPLEs.size());
    for(const auto&entry:gtfs){
        size_t ii=0;
        for(ii=0; ii<_SAMPLEs.size()&&_SAMPLEs[ii].ID!=entry.first;++ii){ }
        samples[ii].clear();
        samples[ii].resize(_FORMATs.size());
        for(const auto&field_entry:entry.second){
            size_t jj=0;
            for(jj=0; jj<_FORMATs.size()&&_FORMATs[jj].ID!=field_entry.first;++jj){ }
            samples[ii][jj]=field_entry.second;
        }
    }
    entries.push_back(line_t(chrom,pos1,id,ref,alts,qual,filter_indices,indexed_infos,samples));
};
void VCFWriter::setSample(const String&name,const String&description){
    size_t ii=0;
    for(ii=0; ii<_SAMPLEs.size()&&_SAMPLEs[ii].ID!=name;++ii){ }
    if( ii == _SAMPLEs.size()){
        _SAMPLEs.push_back(info_t(name,"","",description));
    }
    else{
        _SAMPLEs[ii].set(name,"","",description);
    }
}
void VCFWriter::setFormat(const String&name,const String& number, const String&type, const String&description){
    size_t ii=0;
    for(ii=0; ii<_FORMATs.size()&&_FORMATs[ii].ID!=name;++ii){ }
    if( ii == _FORMATs.size()){
        _FORMATs.push_back(info_t(name,number,type,description));
    }
    else{
        _FORMATs[ii].set(name,number,type,description);
    }
}
