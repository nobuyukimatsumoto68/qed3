#pragma once

#include <chrono>

struct Timer
{
  std::chrono::time_point<std::chrono::system_clock> t_start;
  std::chrono::time_point<std::chrono::system_clock> t_stop;

  void start(){ t_start = std::chrono::system_clock::now(); }
  void restart(){ t_start = std::chrono::system_clock::now(); }
  void stop(){ t_stop = std::chrono::system_clock::now(); }

  Timer()
  {
    start();
  }

  inline double elapsedMilliseconds() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t_stop - t_start).count();
  }
  inline double elapsedSeconds() const { return elapsedMilliseconds() / 1000.0; }

  inline double currentMilliseconds() const {
    const std::chrono::time_point<std::chrono::system_clock> t_current = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t_current - t_start).count();
  }
  inline double currentSeconds() const { return currentMilliseconds() / 1000.0; }

};
