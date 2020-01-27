#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>

class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time_point_;
    unsigned long start_;
    unsigned long end_;
    unsigned long duration_;
    
public:
    Timer ();
    void compute_duration ();
    ~Timer ();  
};

#endif /* _TIMER_H_ */
