#ifndef PX__TIMER_H
#define PX__TIMER_H

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

#endif /* PX__TIMER_H */
