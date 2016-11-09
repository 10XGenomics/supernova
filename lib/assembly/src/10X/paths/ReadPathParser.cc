// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef READPATHPARSER_CC_
#define READPATHPARSER_CC_

#include "10X/paths/ReadPathParser.h"


    /* P S E U D O   C O N S T R U C T O R S */

/**
 * @brief zips a readpath into a readpathX
 *
 * @param bitseq: vec<unsigned char>&
 * @param rp: const ReadPath&
 * @param hb: const HyperBasevectorX&
 */
    void RPParser::LLzip(uCharArray& bitseq, const ReadPath& rp, const HyperBasevectorX& hb){

        int array_size = (rp.size()>0)? (rp.size()-1+3)/4+sizeof(unsigned char)+sizeof(int16_t)+sizeof(uint32_t) : sizeof(unsigned char); 

        int64_t idx = bitseq.size();
        bitseq.resize(bitseq.size()+array_size,0);

        // numEdges
        bitseq[idx] = static_cast<unsigned char>(rp.size()); idx++;

        // append offset
        if(bitseq[idx-1]==0)
            return;
        reinterpret_cast<int16_t*>(&bitseq[idx])[0] = static_cast<int16_t>(rp.getOffset());
        idx += sizeof(int16_t);

        // edge1
        reinterpret_cast<uint32_t*>(&bitseq[idx])[0] = static_cast<uint32_t>(rp[0]);
        idx += sizeof(uint32_t);

        // access edgelist and convert to branchlist
        unsigned int sub_idx = 0;
        for(unsigned int e = 0; e < rp.size()-1; e++){
            int64_t w = hb.ToRight(rp[e]);
            for(unsigned int j = 0; j < hb.From(w).size(); j++){
                if(hb.IFrom(w,j) == rp[e+1]){

                    // encode scheme, insert to bitseq[idx] using current sub_idx
                    LLencodeBranchId(bitseq,j,sub_idx,idx); //FIFO
                    break;
                }
            }
        }
    }
 
/**
 * @brief zips a edge list and an offset into a readpathX
 *
 * @param bitseq is vec<unsigned char> &
 * @param rp: const vec<int> &
 * @param offset: const int
 * @param hb: const HyperBasevectorX&
 */
    void RPParser::LLzip(uCharArray& bitseq, const vec<int>& rp, const int offset, const HyperBasevectorX& hb){
         
        int array_size = (rp.size()>0)? (rp.size()-1+3)/4+sizeof(unsigned char)+sizeof(int16_t)+sizeof(uint32_t) : sizeof(unsigned char);

        int64_t idx = bitseq.size();
        bitseq.resize(bitseq.size()+array_size,0);

        // numEdges
        bitseq[idx] = static_cast<unsigned char>(rp.size()); idx++;

        // append offset
        if(bitseq[idx-1]==0)
            return;
        reinterpret_cast<int16_t*>(&bitseq[idx])[0] = static_cast<int16_t>(offset);
        idx += sizeof(int16_t);

        // edge1
        reinterpret_cast<uint32_t*>(&bitseq[idx])[0] = static_cast<uint32_t>(rp[0]);
        idx += sizeof(uint32_t);

        // access edgelist and convert to branchlist
        unsigned int sub_idx = 0;
        for(unsigned int e = 0; e < rp.size()-1; e++){
            int64_t w = hb.ToRight(rp[e]);
            for(unsigned int j = 0; j < hb.From(w).size(); j++){
                if(hb.IFrom(w,j) == rp[e+1]){

                    // encode scheme, insert to bitseq[idx] using current sub_idx
                    LLencodeBranchId(bitseq,j,sub_idx,idx); //FIFO
                    break;
                }
            }
        }
    }

    /* A C C E S S O R S: U N C O M P R E S S O R S */

/**
 * @brief unzips a readpath at index in bitseq
 *
 * @param bitseq: const vec<unsigned char>&
 * @param rp: ReadPath&
 * @param hb: const HyperBasevectorX&
 * @param index: const int64_t
 */
    void RPParser::LLunzip(const uCharArray& bitseq, ReadPath& rp, const HyperBasevectorX& hb, const int64_t index)const{
        if (bitseq.size()==0)
            exit(1);

        // set numEdges
        rp.resize(bitseq[index],0); 

        // set offset
        if(rp.size()==0){
            rp.setOffset(0);
            return;
        }
        rp.setOffset(LLgetOffset(bitseq,index));

        // set edge1
        rp[0] = LLgetFirstEdge(bitseq,index);

        int64_t idx = index+sizeof(unsigned char) + sizeof(int16_t) + sizeof(uint32_t);
        // access branchid and convert to edgelist
        unsigned int sub_idx = 0;
        for(unsigned int e = 0; e < rp.size()-1; e++){
            int64_t w = hb.ToRight(rp[e]);

            // decode scheme, extract from bitseq[idx] using current sub_idx
            rp[e+1] = hb.IFrom(w,LLdecodeBranchId(bitseq,sub_idx,idx));
        }
    }
    
/**
 * @brief unzips an edge list and offset at index in bitseq
 *
 * @param bitseq: const vec<unsigned char>&
 * @param rp: vec<int>&
 * @param offset: int&
 * @param hb: const HyperBasevectorX&
 * @param index: const int64_t
 */
    void RPParser::LLunzip(const uCharArray& bitseq, vec<int>& rp, int& offset, const HyperBasevectorX& hb, const int64_t index) const{
        if (bitseq.size()==0)
            exit(1);
        
        // set numEdges
        rp.resize(bitseq[index],0); 

        // set offset
        if(rp.size()==0){
            offset = 0;
            return;
        }
        offset = LLgetOffset(bitseq,index);

        // set edge1
        rp[0] = LLgetFirstEdge(bitseq,index);

        int64_t idx = index+sizeof(unsigned char) + sizeof(int16_t) + sizeof(uint32_t);

        // access branchid and convert to edgelist
        unsigned int sub_idx = 0;
        for(unsigned int e = 0; e < rp.size()-1; e++){
            int64_t w = hb.ToRight(rp[e]);

            // decode scheme, extract from bitseq[idx] using current sub_idx
            rp[e+1] = hb.IFrom(w,LLdecodeBranchId(bitseq,sub_idx,idx));
        }
    }


    /* E N C R Y P T I O N   A L G O R I T H M S */

/**
 * @brief encodes a branch id using a fixed width encoding
 *
 * @param bitseq: vec<unsigned char>&
 * @param id: const unsigned int&
 * @param sub_idx: unsigned int&
 * @param id: int64_t&
 * @param sch: CODING_SCHEME
 */
     void RPParser::LLencodeBranchId(uCharArray& bitseq, const unsigned int& id, unsigned int& sub_idx, int64_t& idx, CODING_SCHEME sch){
        /* // This piece of code causes seg-fault via branch prediction- processor does not clean up afterwards!! */
        /* if(idx==bitseq.size()) // resize for rare event in Huffman scheme */
        /*     bitseq.resize(idx+1,0); */
        if(sch == FIXED_WIDTH){
            bitseq[idx] +=((static_cast<unsigned char>(id))<<sub_idx);
            sub_idx += 2;
            // reset indices upon overflow
            if(sub_idx>7){
                idx++;
                sub_idx = 0;
            }
        }
        else{}
    }

/**
 * @brief decodes a branch id using a fixed width encoding
 *
 * @param bitseq: const vec<unsigned char>&
 * @param sub_idx: unsigned int&
 * @param idx: int64_t&
 * @param sch: CODING_SCHEME
 *
 * @return 
 */
    int RPParser::LLdecodeBranchId(const uCharArray& bitseq, unsigned int& sub_idx, int64_t& idx, CODING_SCHEME sch) const{
        int ret = 0;
        unsigned char word = bitseq[idx];
        if(sch == FIXED_WIDTH){
            ret = ( (word >> sub_idx ) & 3 );
            sub_idx += 2;
            // reset indices upon overflow
            if(sub_idx>7){
                idx++;
                sub_idx = 0;
            }
        }
        else{}
        return ret;
    }


    /* A C C E S S O R S  &  M O D I F I E R S */

/**
 * @brief gets an offset at a read index in the bitseq
 *
 * @param bitseq: const vec<unsigned char>&
 * @param idx: const int64_t
 *
 * @return 
 */
    int RPParser::LLgetOffset(const uCharArray& bitseq, const int64_t idx) const{
        if(LLgetNumEdges(bitseq,idx)>0){
            return (static_cast<int>((reinterpret_cast<const int16_t*>(&bitseq[idx+1])[0])));
        }
        else{
            // WARNING: number of aligned edges is 0, so default offset is 0
            return 0;
        }
    }

/**
 * @brief gets the first edge for a read index in the bitseq
 *
 * @param bitseq:const vec<unsigned char>&
 * @param idx: const int64_t
 *
 * @return 
 */
    int RPParser::LLgetFirstEdge(const uCharArray& bitseq, const int64_t idx) const{
        return (static_cast<int>((reinterpret_cast<const uint32_t*>(&bitseq[1+sizeof(int16_t)+idx])[0])));
    }

    int RPParser::LLgetNumEdges(const uCharArray& bitseq, const int64_t idx) const {
        return static_cast<int>(bitseq[idx]);
    }

    int RPParser::LLgetNumBytes(const uCharArray& bitseq) const {
        return bitseq.size();
    }

    /* /1* S T D   O/P *1/ */

/**
 * @brief well formatted output, use only with readpathX
 *
 * @param bitseq: const vec<unsigned char>&
 * @param hb: const HyperBasevectorX&
 * @param messg: std::string
 * @param index: const int64_t
 */
    void RPParser::LLprintMe(const uCharArray& bitseq,const HyperBasevectorX& hb, const std::string messg,const int64_t index) const {
        cout<<messg<<"[";
        int size = LLgetNumEdges(bitseq,index);
        cout<<"offset: "<<LLgetOffset(bitseq,index)<<" | numEdges: "<<size<<" | Edges: ";

        if(size==0){
            cout<<"<empty> ]"<<endl;
            return;
        }

        // set edge1
        int edge = LLgetFirstEdge(bitseq,index);
        cout<<edge;

        int64_t idx = index+sizeof(unsigned char)+sizeof(int16_t)+sizeof(uint32_t);

        if(size==1){
            cout<<"]"<<endl;
            return;
        }

        unsigned int sub_idx = 0;
        for(int e = 0; e < size-1; e++){
            int64_t w = hb.ToRight(edge);

            // decode scheme, extract from bitseq[idx] using current sub_idx
            edge = hb.IFrom(w,LLdecodeBranchId(bitseq,sub_idx,idx));
            cout<<","<<edge;

        }
        cout<<"]"<<endl;
        return;
    }

#endif /* READPATHPARSER_CC_ */
