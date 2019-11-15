// log program calls. useful for documenting and replication

#include <fstream>

#include "log.h"

void log_call (int argc, char * argv[]) {
    std::ofstream phyxlog;
    phyxlog.open ("phyx.logfile", std::ios::out | std::ios::app);
    for (int i = 0; i < argc; i++) {
        phyxlog << argv[i];
        if (i < (argc - 1)) {
            phyxlog << " ";
        }
    }
    phyxlog << std::endl;
    phyxlog.close();
}
