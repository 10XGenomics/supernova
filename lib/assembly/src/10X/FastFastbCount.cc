// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#include "MainTools.h"
#include "feudal/ObjectManager.h"
#include "feudal/VirtualMasterVec.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"


int main(int argc, char *argv[])
{
     RunTime( );

     Bool quiet=False;
     for ( int i = 1; i < argc; ++i ) {
          if ( String(argv[i]) == "QUIET=True" ) {
               quiet = True;
               break;
          }
     }

     parsed_args command(argc, argv, quiet);       // essentially BeginCommandArguments
     CommandArgument_String_OrDefault(FASTB, "");
     CommandArgument_String_OrDefault(QUALP, "");
     CommandArgument_String_OrDefault(QUALB, "");
     CommandArgument_Bool_OrDefault_Doc(QUIET, False, "don't print any decorations, just numbers" );
     CommandArgument_Bool_OrDefault_Doc(MAXLEN, False, "print max read length instead" );
     EndCommandArguments;

     Bool did = False;

     if ( FASTB != "" ) {
          VirtualMasterVec<basevector> bases( FASTB );
          if ( MAXLEN ) {
               auto maxl = bases[0].size();
               for ( size_t i = 0; i < 1000 && i < bases.size(); ++i )
                    maxl = std::max(maxl, bases[i].size());
               cout << (QUIET?"":"#maxlen=") << maxl << endl;
          } else {
               cout << (QUIET?"":"#fastb=") << bases.size() << endl;
          }
          did = True;
     }

     if ( QUALP != "" ) {
          VirtualMasterVec<PQVec> quals( QUALP );
          cout << (QUIET?"":"#qualp=") << quals.size() << endl;
          did = True;
     }

     if ( QUALB != "" ) {
          VirtualMasterVec<qualvector> quals( QUALB );
          cout << (QUIET?"":"#qualb=") << quals.size() << endl;
          did = True;
     }

     if ( !did ) FatalErr("you didn't specify anything to do");

     return 0;
}


