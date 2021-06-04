#include <fstream>
#include <vector>

#include "log.h"

void log_call (int argc, char * argv[]) {
    std::ofstream phyxlog;
    phyxlog.open ("phyx.logfile", std::ios::out | std::ios::app);
    int count = 0;
    for (auto& s : std::vector<char*>(argv, argv+argc)) {
        phyxlog << s;
        if (count < (argc - 1)) {
            phyxlog << " ";
        }
        count++;
    } 
    phyxlog << std::endl;
    phyxlog.close();
}
