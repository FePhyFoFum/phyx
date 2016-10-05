// log program calls. useful for documenting and replication

#include <fstream>

using namespace std;

#include "log.h"

void log_call (int argc, char * argv[]) {
    ofstream phyxlog;
    phyxlog.open ("phyx.logfile", ios::out | ios::app);
    for (int i = 0; i < argc; i++) {
        phyxlog << argv[i];
        if (i < (argc - 1)) {
            phyxlog << " ";
        }
    }
    phyxlog << endl;
    phyxlog.close();
}
