#include <deque>
#include <string>
#include <chrono>
#include <unordered_map>
#include <mutex>
#include <iostream>
#include <thread>
#include <limits.h>

#include "RequestRateLimiter.h"

using namespace std;

typedef chrono::time_point<std::chrono::high_resolution_clock> time_point;

long duration(time_point start, time_point end){
    return (int) chrono::duration_cast<chrono::microseconds>(end - start).count();
}



Record::Record(const time_point start, const int numLock){
    this->start = start;
    this->numLock = numLock;
    this->interval = 0;
}

int Record::getRecentUsage(time_point now, long time) {
    if (this->interval != 0){
        int timeSinceEnd = duration(this->end, now);
        return max(min(this->interval, time-timeSinceEnd), 0L);
    } else {
        return min(duration(this->start, now), time);
    }
}

void Record::setEndTime(time_point end){
    this->end = end;
    this->interval = duration(this->start, this->end);
}

time_point Record::getEnd(){
    return this->end;
}

int Record::getNumLock() {
    return this->numLock;
}

int Record::getInterval() {
    return this->interval;
}

RequestRateLimiter::RequestRateLimiter(const int lockCount, const int maxRequests){
    this->initialLockCount = lockCount;
    this->lockCount = lockCount;
    this->maxRequests = maxRequests;
    this->nextRequestId = 0;

}

void RequestRateLimiter::repeatCleanRequest(int seconds){

    while (true){
        this_thread::sleep_for(std::chrono::seconds(seconds));

        this->lock.lock();
        int newLockCount = this->initialLockCount;
        time_point now = chrono::high_resolution_clock::now();
        for (auto it=records.begin(); it != records.end();){
            unordered_map<int, Record*>& recordsPerIp = it->second;
            for (auto it2=recordsPerIp.begin(); it2 != recordsPerIp.end();){
                Record* record = it2->second;
                if (record->getInterval() == 0){
                    newLockCount--;
                    it2++;
                } else if (duration(record->getEnd(), now) > EXPIRE_TIME){
                    it2 = recordsPerIp.erase(it2);
                } else {
                    it2++;
                }
            }
            if (recordsPerIp.empty()){
                it = records.erase(it);
            } else {
                it++;
            }
        }
        cout << "clean request (lockCount: " << this->lockCount << " -> " << newLockCount << ")" << endl;
        this->lockCount = newLockCount;
        this->lock.unlock();
    }
}

bool RequestRateLimiter::checkTimeCriteria(const string& ipAddress, const time_point& now){

    auto it = records.find(ipAddress);
    if (it == records.end()){
        return true;
    }
    unordered_map<int, Record*>& recordsPerIp = it->second;

    for (int i=0; i<NUM_TIME_CONSTRAINTS; i++) {
        int timeUsed = 0;
        for (pair<int, Record*> pair1: recordsPerIp){
            timeUsed += pair1.second->getRecentUsage(now, MEASURE_TIME[i]);
        }
//        cout << (double) timeUsed / 1000000 << " " << (double) MEASURE_TIME[i] / 1000000 << " " << (double) ALLOWED_TIME[i] / 1000000 << endl;
        if (timeUsed > ALLOWED_TIME[i]){
            return false;
        }
    }

    return true;

}

int RequestRateLimiter::acquireLock(const std::string &ipAddress, int numLock) {

    // check if num is valid
    if (numLock > this->maxRequests){
        cout << "Number of lock being acquired is bigger than the number of max requests." << endl;
        return -3;
    }

    time_point now = chrono::high_resolution_clock::now();

    // check for exceeding computation time per IP address
    if (!checkTimeCriteria(ipAddress, now)){
        return -2;
    }

    // Check if lock is available
    this->lock.lock();
    if (this->lockCount < this->maxRequests){
        this->lock.unlock();
        return -1;
    }
    this->lockCount -= numLock;
    int requestId = this->nextRequestId;
    if (this->nextRequestId == INT_MAX){
        this->nextRequestId = 0;
    } else {
        this->nextRequestId++;
    }
    this->lock.unlock();
    cout << "lock acquired (requestId: " << requestId << ", IP address: " << ipAddress << ", num: "  << numLock << ") lockCount: " << this->lockCount << endl;
    this->records[ipAddress][requestId] = new Record(now, numLock);

    return requestId;
}


void RequestRateLimiter::releaseLock(const string& ipAddress, int requestId){
    time_point now = chrono::high_resolution_clock::now();
    this->lock.lock();
    Record* record = records[ipAddress][requestId];
    record->setEndTime(now);
    this->lockCount += record->getNumLock();
    this->lock.unlock();
    cout << "lock released (requestId: " << requestId << ", IP address: " << ipAddress << ", num: "  << record->getNumLock() << ") lockCount: " << this->lockCount;
    cout << " interval: " << (double) record->getInterval() / 1000000 << " seconds" << endl;
}
