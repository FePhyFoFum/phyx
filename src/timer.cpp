#include <chrono>
#include <iostream>

#include "timer.h"

// a barebones RAII timer


Timer::Timer ():start_(0), end_(0), duration_(0) {
    start_time_point_ = std::chrono::high_resolution_clock::now();
}

void Timer::compute_duration () {
    end_time_point_ = std::chrono::high_resolution_clock::now();
    start_ = std::chrono::time_point_cast<std::chrono::microseconds>(start_time_point_).time_since_epoch().count();
    end_ = std::chrono::time_point_cast<std::chrono::microseconds>(end_time_point_).time_since_epoch().count();
    duration_ = end_ - start_;
    double millisecs = duration_ * 0.001;

    std::cout << duration_ << " µs (" << millisecs << " ms)" << std::endl;
}

Timer::~Timer () {
    compute_duration();
}
