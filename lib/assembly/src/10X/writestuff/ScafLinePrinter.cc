// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#include "10X/writestuff/ScafLinePrinter.h"

void ScafLinePrinter::FindIncorrectSeqInsertions( HyperBasevectorX const& hbd2, vec<bool>& baddies )
{
     ForceAssertGt(_to_left.size(), 0 );       // ensure initialized
     ForceAssertGt(_to_right.size(), 0 );      // ensure initialized

     int count = 0;
     baddies.resize(_D.E(), false);
     for ( int d = 0; d < _D.E(); ++d ) {
          if ( IsSequence( _D.O(d) ) && hbd2.O(d).size() == 0 ) {
               baddies[d] = true;
               count++;
          }
     }

}

void ScafLinePrinter::ExpandDGraphCellEdge( FastaEdgeWriter& fw, int edge_id, cell const& this_cell )
{
     // we start with original vertex ids and switch to 
     // D2 ids if/when we have a valid path through
     int left = _to_left[edge_id];
     int right = _to_right[edge_id];

     vec<int> path;
     this_cell.FindPath(path);

     if ( path.size() == 0 ) {
          fwLog(fw,"\t\t\t\tcell has NO path through -> GAP" << endl);
          fw.AddGapEdge( left, right, edge_id, 10 );
          return;
     }


     fwLog(fw, "*** tail is: " << fw.Tail() );

     fwLog(fw, endl << "*** first edge of path: ");
     int cell_edge0_id = path[0];
     vec<int> cell_edge0 = this_cell.G().O(cell_edge0_id);

     fwLog(fw, "cell edge0 base graph edge sequence: " << printSeq(cell_edge0) << endl);
     _hb.O( cell_edge0[0] ).Print(fw.Log());

     fwLog(fw, endl);

     vec<int> base_path;
     for ( int base_edge_id : path ) {
          vec<int> const& base_edge = this_cell.G().O(base_edge_id);
          base_path.append( base_edge );
     }

     try { 
          for ( size_t i = 0; i < path.size(); ++i ) {
               path[i] += _dcell_map.at(edge_id);
          }
     } catch ( std::out_of_range& ) {
          FatalErr("BUG: edge_id missing from cell translation map");    // TODO
     }

     // we pull vertices from the reinserted graph if we have a valid path
     left=_to_left2[path.front()];
     right=_to_right2[path.back()];

     if ( base_path.size() == 1 && base_path[0] < 0 ) {
          fw.AddGapEdge(left, right, path[0]);
          return;
     }

     std::string seq;

     auto si = base_path.begin();
     if ( *si < 0 ) fw.PushTail();      // if we start with a gap, push last K-1 first
     do {
          // find runs prior to internal gaps
          auto ei = std::find( si, base_path.end(), -1 );
          if ( ei != si ) {
               vec<int> vtmp(ei-si);
               std::copy( si, ei, vtmp.begin() );
               auto bseq = _hb.Cat(vtmp);
               std::transform( bseq.begin(), bseq.end(), 
                         std::back_inserter(seq), BaseToCharMapper() );
          }
          // add 100 N's if an internal gap was found
          if ( ei != base_path.end() ) {
               seq.append(100, 'N');
               ei++;
          }
          si = ei;
     } while ( si != base_path.end() );

     fwLog(fw,"base path sequence: " << seq << endl );


     fw.AddEdgeSeq(left, right, path, seq);
     fw.PushTail();
     if ( seq.back() != 'N') {
          fw.AppendAnonGap(10,fw.K()-1);
     }

}

void ScafLinePrinter::ExpandDGraphGapEdge( FastaEdgeWriter& fw, int edge_id  )
{
     const int bc_gap_repr = 3000;
     auto edge_path = _D.O(edge_id);
     int left = _to_left2[edge_id];
     int right = _to_right2[edge_id];
     ForceAssertGt(edge_path.size(), 0u);
     if ( IsPairGap( edge_path ) )
          fw.AddGapEdge( left, right, edge_id );
     else if ( IsBarcodeOnlyGap( edge_path ) )
     {    int gap = bc_gap_repr;
          if ( edge_path.size( ) >= 2 ) gap = edge_path[1];
          fw.AddGapEdge( left, right, edge_id, gap );    }
     else if ( IsSequence( edge_path ) ) {
          if ( _baddies[edge_id] ) {
               fwLog(fw, "\t\t\tedge_id=" << edge_id << " marked as a baddie.. adding gap" << endl);
               fw.AddGapEdge(left, right, edge_id);
          } else {
               int ltrim, rtrim;
               basevector seq;
               GapToSeq( edge_path, ltrim, rtrim, seq );
               fwLog(fw, endl << "\t\t\t ltrim=" << ltrim << ", rtrim=" << rtrim << ", seq=" );
               seq.Print(fw.Log());
               fwLog(fw, "\t\t\t _seq right=");
               auto tmp_seq = fw.Seq();
               if ( fw.IsLog() ) {
                    for ( size_t i = tmp_seq.size() - ltrim - (_hbd.K()-1); i < tmp_seq.size() ; ++i )
                         fw.Log() << tmp_seq[i];
                    fw.Log() << ", tail=" << fw.Tail();
                    fw.Log() << "\t\t\t";
               }
               // check tail and add a gap if we don't match
               auto tail = fw.Tail();
               ForceAssertGe( seq.size(), tail.size() );
               string head_seq;
               std::transform( seq.begin(), seq.begin() + tail.size(), 
                         std::back_inserter(head_seq), BaseToCharMapper() );
               if ( tail == head_seq ) {
                    fwLog(fw, "GOOD sequence insert part of record " << fw.Count()+1 << endl << "\t\t\t");
                    fw.AddEdge( left, right, edge_id, seq );
               } else {
                    fwLog(fw, "BAD mismatched sequence insert.... adding gap as record " << fw.Count()+1 << endl << "\t\t\t");
                    fw.AddGapEdge(left, right, edge_id);
               }
               fwLog(fw, endl);
          }
     } else if ( IsCell( edge_path ) ) {
          // fwLog(fw, "\t\t\tCoding cell as gap" << endl);
          // fw.AddGapEdge( left, right, edge_id, 10 );
          fwLog(fw, "\t\t\tExpanding Cell" << endl);
          cell this_cell;
          this_cell.CellDecode(edge_path);
          ExpandDGraphCellEdge( fw, edge_id, this_cell );
          fwLog(fw, "\n");
     } else
          FatalErr( "BUG: unrecognized gap type in graph" );

}

void ScafLinePrinter::ExpandDLinePath( FastaEdgeWriter& fw, vec<int> const& path ) 
{
     fwLog(fw,"\t\t");
     for ( int edgei = 0; edgei < path.isize(); ++edgei ) {
          auto edge = path[edgei];
          int left = _to_left2[edge];
          int right = _to_right2[edge];
          if ( _hbd.O(edge).size() > 0 ) {
               fw.AddEdge( left, right, edge, _hbd.O(edge) );
               fwLog(fw," {" << edge << "} len so far=" << fw.Seq().size());
          } else {
               fwLog(fw, " {GAP: " << _D.O(edge)[0] << "} ");
               ExpandDGraphGapEdge(fw, edge);
          }
     }
     fwLog(fw, endl);
}

void ScafLinePrinter::ExpandDLineGapCell( FastaEdgeWriter& fw, vec<vec<vec<int>>> const& line, int const celli )
{
     // recover gap edge from D graph
     ForceAssertGt(celli, 0);
     ForceAssertLt(celli, line.size()-1 );
     ForceAssertEq(celli%2, 1);

     // vertex left of gap should be right of last edge of previous (even-numbered) cell
     auto vleft = _to_right[ line[celli-1][0].back() ];
     auto vright = _to_left[ line[celli+1][0].front() ];

     ForceAssertEq( _D.FromSize(vleft), 1 );
     auto gap_edge_id = _D.EdgeObjectIndexByIndexFrom( vleft, 0 );
     auto gap_edge = _D.O(gap_edge_id);

     if ( gap_edge[0] < 0 ) {
          fwLog(fw, "\t\t{*GAP: " << gap_edge[0] << " -> " << gap_edge_id << "} " << endl);
     } else {
          FatalErr("MakeFasta BUG: pulled back a gap edge that was not a gap");
     }

     ExpandDGraphGapEdge( fw, gap_edge_id );
}

int ScafLinePrinter::ExpandDLineCell( FastaEdgeWriter& fw, vec<vec<int>> const& cell )
{
     ForceAssertGt( cell.size(), 0u );
     vec<int> idx( cell.size(), vec<int>::IDENTITY );

     if ( cell.size() > 1 && _breakBubbles ) {
          int left = _to_left[ cell[0][0] ];
          for ( size_t i = 1; i < cell.size(); ++i )
               ForceAssertEq( left, _to_left[ cell[i][0] ] );
          fw.Break();
          for ( int pathi = 0; pathi < cell.isize(); ++pathi ) {
               ExpandDLinePath( fw, cell[pathi] );
               fw.Break();
          }

          return _to_right[cell.back().back()];
     } 
     
     
     if ( cell.size() > 1 && !_breakBubbles ) {
          // we have a coverage-based choice to make & try to avoid gap edges
          vec<pair<int,int>> gapcov( cell.size( ), make_pair(1,0) );
          for ( int pathi = 0; pathi < cell.isize(); ++pathi ) {
               for ( auto edge : cell[pathi] )
               {    const vec<int>& x = _D.O(edge);
                    if ( IsPairGap(x) || IsBarcodeOnlyGap(x) || IsCell(x) )
                    {    gapcov[pathi].first = 0;    }
                    gapcov[pathi].second += _dedge_counts[edge];
               }
          }
          ReverseSortSync(gapcov, idx);
          // fwLog(fw, "\t\texpanding cell, coverages=" << printSeq(cov) << endl);
     }

     ExpandDLinePath( fw, cell[ idx[0] ] );
     return _to_right[cell[idx[0]].back()];
}

void ScafLinePrinter::ExpandDLine( FastaEdgeWriter& fw, int linei )
{
     auto const& line = _dlines[linei];
     fwLog(fw, "\texpanding L" << linei  << " of size " << line.size() << " cells " << endl);

     if ( line.size() == 1 && line[0].size() == 1 && line[0][0].size() == 0 ) {
          FatalErr("BUG: {special empty cell type 1}");
     }

     for ( size_t celli = 0; celli < line.size(); ++celli ) {
          fwLog(fw, "\tcell " << celli << ", size=" << line[celli].size() << endl);
          if ( line[celli].size() == 1 && line[celli][0].size() == 0 ) {
               if ( celli == line.size()-1) FatalErr("MakeFasta BUG: gap cell at end of non-empty line");
               if ( celli == 0 ) FatalErr("MakeFasta BUG: gap cell at start of line");
               fwLog(fw, "\t\t(expanding gap cell)" << endl);
               ExpandDLineGapCell(fw, line, celli );
          } else {
               ExpandDLineCell( fw, line[celli] );
          }
     }

}

void ScafLinePrinter::ExpandMegabubbleArm( FastaEdgeWriter& fw,  vec<int> const& arm)
{
     // expand a path of lines in a megabubble (i.e. a path in a cell of a lines-of-lines)
     // so this is a scaffolding of lines from dlines
     for ( auto const linei : arm ) {
          // pull out line and walk 
          ExpandDLine( fw, linei );
     }
     if ( !_mashMegaBubbles) fw.Break();
}

void ScafLinePrinter::BustMegabubble( FastaEdgeWriter& fw, vec<vec<int>> const& llcell ) {
     // dump any previous sequence
     int vleft = _to_left[ llcell.front().front() ];
     fw.BreakIfSeq();
     // just output each line in the cell as separate records
     vec<int> lines;
     for ( int i = 0; i < llcell.isize(); ++i )
          for ( int j = 0; j < llcell[i].isize(); ++j ) 
               lines.push_back(llcell[i][j]);
     UniqueSort(lines);

     for ( auto const linei: lines )  {
          fw.Index();
          ExpandDLine( fw, linei );
          fw.Break();  // bustable megabubbles are odd cells, so can't be terminal
     }
}


void ScafLinePrinter::WalkScaffoldLines( FastaEdgeWriter& fw, size_t choose /* = 0 */ ) {

     if ( _mashMegaBubbles && _breakBubbles )
          FatalErr("can't mash megabubbles but break bubbles");

     // walk the lines of lines graph
     int ndots=0, stopped=0;
     cout << "writing 100 dots in 10 groups of 10:" << endl;
     for ( int linei = 0; linei < _dlines2.isize(); ++linei ) {
          MakeDots(stopped, ndots, _dlines2.isize() );

          // skip RC lines if we're not keeping them
          if ( _keepRc == false && _linv2[linei] < linei ) {
               fwLog(fw, "Scaffold Line " << linei << " RC of line " << _linv2[linei] << " - skipping" << endl);
               continue;
          } else if ( _linv2[linei] == linei ) {
               fwLog(fw, "Scaffold Line " << linei << " RC of itself" << endl);
          } else {
               fwLog(fw, "Scaffold Line " << linei << " RC of line " << _linv2[linei] << " - keeping" << endl);
          }

          auto& line=_dlines2[linei];
          ForceAssertGt(line.size(), 0u);
          for ( size_t celli = 0; celli < line.size(); ++celli ) {
               auto& cell=line[celli];
               ForceAssertGt( cell.size(), 0 );
               ForceAssert( celli % 2 != 0 || cell.size() == 1 );
               if ( cell.size() > 2 ) {
                    fwLog(fw, "Scaffold Line " << linei << ", busting many-arm (" << cell.size() << ") cell " << celli << endl);
                    BustMegabubble(fw, cell);
               } else {
                    size_t ncell = _mashMegaBubbles ? 1 : cell.size();
                    size_t start_cell = _mashMegaBubbles ? std::min(choose,cell.size()-1) : 0;
                    size_t end_cell = start_cell + ncell;
                    fwLog(fw, "Scaffold Line " << linei << ", expanding cell " << celli << endl);
                    for ( size_t armi = start_cell; armi < end_cell; ++armi )  {
                         fwLog(fw, "Scaffold Line " << linei << ", expanding cell " << celli <<
                             ", megabubble cell arm " << armi+1 << " of " << cell.size() << endl);
                         fw.Index();
                         ExpandMegabubbleArm(fw, cell[armi]);
                    }
               }
          }
          if ( _mashMegaBubbles ) fw.Break();
     }
}

