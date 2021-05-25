#ifndef PX_TIMER_H
#define PX_TIMER_H

#include <chrono>

class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time_point_;
    long start_;
    long end_;
    long duration_;
    
public:
    Timer ();
    void compute_duration ();
    ~Timer ();  
};

#endif /* PX_TIMER_H */
