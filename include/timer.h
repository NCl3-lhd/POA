#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

class Timer {
public:
  static Timer& instance() {
    static Timer inst;
    return inst;
  }

  void start(const std::string& name) {
    starts_[name] = std::chrono::high_resolution_clock::now();
  }

  void stop(const std::string& name) {
    auto it = starts_.find(name);
    if (it == starts_.end()) return;
    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - it->second).count();
    records_[name].push_back(ms);
    starts_.erase(it);
  }

  void print() const {
    if (records_.empty()) return;
    std::cerr << "\n=== Timer ===\n";
    std::vector<std::pair<std::string, Record>> sorted(records_.begin(), records_.end());
    std::sort(sorted.begin(), sorted.end(), [](const std::pair<std::string, Record> &a, const std::pair<std::string, Record> &b) {
      return total(b.second) > total(a.second);
    });
    double grand = 0;
    for (auto& [name, rec] : sorted) {
      double t = total(rec);
      grand += t;
      std::cerr << std::left << std::setw(24) << name
                << std::right << std::fixed << std::setprecision(2)
                << std::setw(10) << t << " ms"
                << "  (cnt=" << rec.size()
                << "  avg=" << (rec.empty() ? 0 : t / rec.size())
                << ")\n";
    }
    std::cerr << std::string(50, '-') << "\n"
              << std::left << std::setw(24) << "TOTAL"
              << std::right << std::fixed << std::setprecision(2)
              << std::setw(10) << grand << " ms\n";
  }

  void reset() {
    records_.clear();
    starts_.clear();
  }

private:
  Timer() = default;
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;
  using Record = std::vector<double>;

  static double total(const Record& r) {
    double s = 0;
    for (double v : r) s += v;
    return s;
  }

  std::unordered_map<std::string, TimePoint> starts_;
  std::unordered_map<std::string, Record> records_;
};

#endif
