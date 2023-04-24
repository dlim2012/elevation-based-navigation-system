#include <deque>
#include <string>
#include <chrono>
#include <unordered_map>

using namespace std;


//const int MEASURE_TIME = 120000; // 2 minutes;
//const int ALLOWED_TIME = (int) (MEASURE_TIME * 3.0);
const int NUM_TIME_CONSTRAINTS = 4;

const int MEASURE_TIME [NUM_TIME_CONSTRAINTS] {
        3000000, // 3 seconds
        10000000, // 10 seconds
        30000000, // 30 seconds
        60000000, // 60 seconds
//        120000000 // 120 seconds
};

const int ALLOWED_TIME [NUM_TIME_CONSTRAINTS] {
    3000000 * 6,
    10000000 * 6,
    30000000 * 5,
    60000000 * 4,
//    120000000 * 6
};

const int EXPIRE_TIME = MEASURE_TIME[NUM_TIME_CONSTRAINTS-1]; // 4 minutes

typedef chrono::time_point<std::chrono::high_resolution_clock> time_point;

int duration(time_point start, time_point end);

class Record{
private:
    time_point start;
    time_point end;
    int interval;
    int numLock;
public:
    Record(time_point start, int numLock);
    int getRecentUsage(time_point now, int time);
    void setEndTime(time_point end);
    time_point getEnd();
    int getNumLock();
    int getInterval();
};

class RequestRateLimiter {
private:
    int initialLockCount;
    int lockCount;
    int maxRequests;
    int nextRequestId;

    // A lock for acquiring permissions
    mutex lock;

    // records per IP Address
    // { IP Address: { request ID : Record } }
    unordered_map<string, unordered_map<int, Record*>> records;

public:
    RequestRateLimiter(int lockCount, int maxRequests);

    // clean-up expired records or locks and reset for any potential inconsistencies
    void repeatCleanRequest(int seconds);

    // check time criteria for an IP address
    bool checkTimeCriteria(const string& ipAddress, const time_point& now);

    // acquire permission: return request ID
    int acquireLock(const string& ipAddress, int numLock);

    // release lock and add to record
    void releaseLock(const string& ipAddress, int requestId);

};