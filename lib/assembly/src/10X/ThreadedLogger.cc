#include "10X/ThreadedLogger.h"
#include "MainTools.h"

int main()
{
     ThreadedLogger log("neil.txt", 4);

     log.get(3) << "this is three" << endl << endl;
     log.get(0) << "this is zero" << endl;
     log.get(1) << "this is one" << endl << endl;
     log.join();
     log.get(2) << "this is two" << endl << endl;
}
