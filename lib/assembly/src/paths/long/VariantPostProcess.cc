///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * VariantPostProcess.cc
 *
 *  Created on: Sep 19, 2013
 *      Author: blau
 */


#include "paths/long/VariantPostProcess.h"
#include <cmath>
#include <map>
#include <fstream>
#include "TokenizeString.h"
#include <cstdio>
#include <functional>


//below is a line-by-line translation of translation of paths/long/scripts/filter_vcfs.py at r47584

//# filter_vcfs - filter variant lines where phred score in either REFP or
//# ALTP exceeds a threshold.
//#
//# The following rules apply:
//#
//# - already filtered lines are left alone (lines not marked as "PASS" or
//# ".")
//#
//# - lines that are "PASS" or "." with a phred score for any element in
//# ALTP or REFP which is >0 AND <threshold are marked as filtered (e.g.
//# "PhredFilter10" for lines filtered with a threshold of 10.
//#
//# - lines that are "." and pass the phred threshold are marked as
//# "PASS"
//#

//import sys
//import argparse
//import math


//def filter_vcf( input, output ):
bool filter_vcf(const String&input,const String& output){
    bool bError=false;
    
//  phred90 = int(-10.0*math.log10( 1.0-0.9 ))
    const int phred90 = int(-10.0*log10( 1.0-0.9 ));
//  phred995 = int(-10.0*math.log10( 1.0-0.995 ))
    const int phred995 = int(-10.0*log10( 1.0-0.995 ));
//  print "input={}, output={}, phred90={}, phred995={}".format(input, output, phred90, phred995)
    std::cout << "input=" << input << ", output=" << output << ", phred90=" << phred90 << ", phred995=" << phred995 << std::endl;

//  # new filters must be added below as symbol=("name","description")
//  # this is just to enforce consistency between what's in the head and
//  # in the body of the vcf
/*
  filters = dict(
      Filter090=( "PhredFilter{}".format(phred90), "Some allele was below prob 0.90 = phred {}".format(phred90) ),
      NoAllelesPass0995=( \
          "NoAllelesPassPhred{}".format(phred995), "Not one allele was above prob 0.995 = phred {}".format(phred995) ),
      RefCallOnly=( "RefCallOnly", "Non-variant line; only the reference was called" ),
      TooManyCalls=( "TooManyCalls", "Multi-allelic site (greater than two alleles)" )
  )
*/
    std::map<String,std::pair<String,String>> filters;
    filters["SFilter090"]=std::make_pair("PhredFilter"+ToString(phred90)
                                       ,"Some allele was below prob 0.90 = phred "+ToString(phred90));
    filters["NoAllelesPass0995"]=std::make_pair("NoAllelesPassPhred"+ToString(phred995)
                                               ,"Not one allele was above prob 0.995 = phred "+ToString(phred995));
    filters["RefCallOnly"]=std::make_pair("RefCallOnly","Non-variant line; only the reference was called");
    filters["TooManyCalls"]=std::make_pair("TooManyCalls","Multi-allelic site (greater than two alleles)");



//  with open(output, 'w') as fd:
    ofstream os(output,ios::out);
    ifstream is(input,ios::in);
    if( os.is_open() && is.is_open()){
//      for line in open(input, 'r'):
        String line;
        while(getline(is,line)){
//          if line[0] == '#':
            if( line[0] == '#'){
//              if line[:len("#CHROM")] == "#CHROM":
                if (line.substr(0,String("#CHROM").size()) == "#CHROM"){
//                  # dump out the ##FILTER lines before the #CHROM
//                  for id, descr in filters.viewvalues():
                    for( const auto& entry: filters){
//                      fd.write("##FILTER=<ID={id},Description=\"{descr}\">\n".format( id=id, descr=descr ))
                        os << "##FILTER=<ID="<<entry.second.first<<",Description=\""<<entry.second.second<<"\">"<<std::endl;
                    }
                }
//              fd.write(line)
                os<<line<<std::endl;
            }
//          else:
            else{
//              fields=line.split()
                vec<String> fields; Tokenize(line,fields);
//              (chrom,pos,id,ref,alt,qual,filter,info,format)=fields[:9]
                String chrom=fields[0], pos=fields[1], id=fields[2], ref=fields[3], alt=fields[4], qual=fields[5], filter=fields[6], info=fields[7], format=fields[8];
//              samples = fields[9:]
                vec<String> samples(fields.begin()+9,fields.end());
//              if len(samples) > 1:
//                  raise Exception("Does not yet work for multi-sample data.")
                if ( samples.size() > 1 ){
                    std::cout << "Does not yet work for multi-sample data." << std::endl;
                    bError=true;
                    break;
                }


//              # don't filter unless we're already passing or not yet filtered
//              if filter == '.' or filter == 'PASS':
                vec<String> new_samples;
                if (filter=="." or filter == "PASS"){
//                    names = format.split(':')
                    vec<String> names; Tokenize(format,':',names);
//                  if (not "REFP" in names) or (not "ALTP" in names):
//                      raise Exception("missing REFP and ALTP tags in line {}".format(line.strip() ) )
                    if (!Member(names,String("REFP"))  || !Member(names,String("ALTP"))){
                        std::cout << "missing REFP and ALTP tags in line " << line << std::endl;;
                        bError=true;
                        break;
                    }

//                  new_samples=[]
                    new_samples.clear();
//                  for sample in samples:
                    for( const auto& sample:samples){
//                      vals = sample.split(':')
                        vec<String> vals; Tokenize(sample,':',vals);

//                      if len(vals) != len(names):
//                          raise Exception("sample values {} doesn't match format {} for line {}".format(
//                              sample, format, line.strip()))
                        if(vals.size()!=names.size()){
                            std::cout << "sample values " << sample << " doesn't match format " << format << "for line " << line << std::endl;;
                            bError=true;
                            break;
                        }

//                      sample_info = dict( zip( names, vals ) )
                        std::map<String,String> sample_info;
                        for(size_t ii=0;ii<names.size();++ii) sample_info[names[ii]]=vals[ii];

//                      probs = [ float( sample_info['REFP'] ) ]
                        vec<double> probs(1,sample_info["REFP"].Double());
//                      probs.extend( map(float, sample_info['ALTP'].split(',') ) )
                        {
                            vec<String> tmp; Tokenize(sample_info["ALTP"],',',tmp);
                            std::transform(tmp.begin(),tmp.end(),back_inserter(probs),std::mem_fn(&String::Double));
                        }

//                      if True in map( lambda x: x>0 and x<phred90, probs ):
                        if ( any_of(probs.begin(),probs.end(),[&phred90](double x){return x>0 && x<phred90;} )){
//                          new_samples = samples
                            new_samples = samples;
//                          filter=filters['Filter090'][0]
                            filter=filters["SFilter090"].first;
//                          break
                            break;
                        }


//                      calls = map( lambda x: x > phred995, probs )
                        vec<unsigned int> calls;
                        std::transform(probs.begin(),probs.end(),back_inserter(calls),[&phred995](double x){return (unsigned int)(x>phred995);});

//                      if sum(calls) == 0:
                        if ( std::accumulate(calls.begin(),calls.end(),0) == 0 ){
//                          new_samples = samples
                            new_samples = samples;
//                          filter=filters['NoAllelesPass0995'][0]
                            filter=filters["NoAllelesPass0995"].first;
//                          break
                            break;
                        }

//                      if sum(calls) == 1 and calls[0] == True:
                        if ( std::accumulate(calls.begin(),calls.end(),0) == 1 && calls[0] > 0){
//                          new_samples = samples
                            new_samples = samples;
//                          filter=filters['RefCallOnly'][0]
                            filter=filters["RefCallOnly"].first;
//                          break
                            break;
                        }

//                      if sum(calls) > 2:
                        if ( std::accumulate(calls.begin(),calls.end(),0) > 2 ){
//                          new_samples = samples
                            new_samples = samples;
//                          filter=filters['TooManyCalls'][0]
                            filter=filters["TooManyCalls"].first;
//                          break
                            break;
                        }

//                      allele1 = calls.index(True)
                        size_t allele1 = find_if(calls.begin(),calls.end(),[](unsigned int x){return x!=0;})-calls.begin();
                        size_t allele2;

//                      if sum(calls) == 1:
                        if ( std::accumulate(calls.begin(),calls.end(),0) == 1 ){
//                          allele2 = allele1
                            allele2 = allele1;
                        }
//                      else:
                        else{
//                          allele2 = calls.index(True, allele1+1)
                            allele2 = find_if(calls.begin()+allele1+1,calls.end(),[](unsigned int x){return x!=0;})-calls.begin();
                        }


//                      gt_idx = names.index('GT')
                        size_t gt_idx = find_if(names.begin(),names.end(),[](String x){return x=="GT";})-names.begin();
//                      vals[gt_idx] = '{}/{}'.format(allele1, allele2)
                        vals[gt_idx] = ToString(allele1)+"/"+ToString(allele2);

//                      new_samples.append( ':'.join(vals) )
                        new_samples.push_back(vals[0]);
                        for(size_t ii=1;ii<vals.size();++ii){ new_samples.back().append(":"+vals[ii]); }
                    }
                    if(bError) break;
                }


//              # make passing "." lines a "PASS"
//              if filter == '.': filter = 'PASS'
                if (filter == ".") filter = "PASS";

//              fd.write( "\t".join( [chrom, pos, id, ref, alt, qual, filter, info, format] ) )
                os << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
                   << qual << "\t" << filter << "\t" << info << "\t" << format;
//              fd.write( "\t" )
//              fd.write( "\t".join( new_samples ) )
                for(const auto& entry:new_samples){ os << "\t" << entry;}
//              fd.write( "\n" )
                os<<std::endl;
            }
        }
        os.close();
        is.close();
        if(bError){
            remove(output.c_str());
        }
    }
    else{
        if( !is.is_open()){
            std::cout << "failed to open " << input << std::endl;
        }
        if( !os.is_open()){
            std::cout << "failed to open " << output << std::endl;
        }
        else{
            remove(output.c_str());
        }
        bError=true;
    }
    return bError;
}

/*
def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='filter Discovar-generated VCF based on REFP/ALTP')
    parser.add_argument( '-o', '--output', help='VCF output file', required=True)
    parser.add_argument( 'input', help='VCF input file' )
    args = parser.parse_args(argv[1:])

    return(filter_vcf( args.input, args.output ) )

if __name__ == "__main__":
    sys.exit(main(sys.argv))
*/
