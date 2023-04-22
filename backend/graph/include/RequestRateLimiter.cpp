#include <deque>
#include <string>
#include <chrono>
#include <unordered_map>
#include <mutex>
#include <iostream>

#include "RequestRateLimiter.h"

using namespace std;

typedef chrono::time_point<std::chrono::high_resolution_clock> time_point;

int duration(time_point start, time_point end){
    return (int) chrono::duration_cast<chrono::milliseconds>(end - start).count();
}


TimeInterval::TimeInterval(time_point start, time_point end){
    this->start = start;
    this->end = end;
    this->intervalInMilliseconds = duration(start, end);
}

RequestRateLimiter::RequestRateLimiter(int lockCount){
    this->lockCount = lockCount;
}

void RequestRateLimiter::cleanRequests(){
    time_point now = chrono::high_resolution_clock::now();
    for (auto it=requests.begin(); it!=requests.end();){
        deque<TimeInterval>& timeIntervals = it->second;
        while (!timeIntervals.empty() && duration(timeIntervals.front().end, now) > MEASURE_TIME_1){
            timeIntervals.pop_front();
        }
        if (timeIntervals.empty()){
            it = requests.erase(it);
        } else {
            it++;
        }
    }
}

bool RequestRateLimiter::isAllowed(const string& ipAddress, time_point& now){
    auto it = requests.find(ipAddress);
    if (it == requests.end()){
        return true;
    }
    deque<TimeInterval>& timeIntervals = it->second;
    cout << ipAddress << " " << timeIntervals.size() << endl;


    // remove all expired records
    while (!timeIntervals.empty() && duration(timeIntervals.front().end, now) > MEASURE_TIME_1){
        timeIntervals.pop_front();
    }

    int timeUsed = 0;
    for (TimeInterval& timeInterval: timeIntervals){
        cout << duration(timeInterval.end, now) << " " << duration(timeInterval.start, now) << endl;
        if (duration(timeInterval.start, now) > MEASURE_TIME_1){
            timeUsed += MEASURE_TIME_1 - duration(timeInterval.end, now);
        } else {
            timeUsed += timeInterval.intervalInMilliseconds;
        }
    }

    return timeUsed < ALLOWED_TIME_1;
}

void RequestRateLimiter::addRecord(const string& ipAddress, time_point& start, time_point& end){
    requests[ipAddress].emplace_back(start, end);
}

bool RequestRateLimiter::acquireLock(){
    this->lock.lock();
    cout << this->lockCount << endl;
    if (this->lockCount == 0){
        this->lock.unlock();
        return false;
    } else {
        this->lockCount -= 1;
        this->lock.unlock();
        return true;
    }

}

void RequestRateLimiter::releaseLock(){
    this->lock.lock();
    this->lockCount += 1;
    cout << this->lockCount << endl;
    this->lock.unlock();
}