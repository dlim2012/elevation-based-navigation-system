#include <deque>
#include <string>
#include <chrono>
#include <unordered_map>

using namespace std;


//const int MEASURE_TIME = 120000; // 2 minutes;
//const int ALLOWED_TIME = (int) (MEASURE_TIME * 3.0);
const int NUM_TIME_CONSTRAINTS = 0;

const long MEASURE_TIME [NUM_TIME_CONSTRAINTS] {
//        3000000, // 3 seconds
//        10000000, // 10 seconds
//        30000000, // 30 seconds
//        60000000, // 60 seconds
//        120000000 // 120 seconds,
//        3600000000 // 1 hour
};

const long ALLOWED_TIME [NUM_TIME_CONSTRAINTS] {
//    3000000 * 4,
//    10000000 * 3,
//    30000000 * 5 / 2,
//    60000000 * 2,
//    120000000 * 9 / 5
};

const int EXPIRE_TIME = NUM_TIME_CONSTRAINTS > 0 ? MEASURE_TIME[NUM_TIME_CONSTRAINTS-1] : 0;

typedef chrono::time_point<std::chrono::high_resolution_clock> time_point;

long duration(time_point start, time_point end);

class Record{
private:
    time_point start;
    time_point end;
    long interval;
    int numLock;
public:
    Record(time_point start, int numLock);
    int getRecentUsage(time_point now, long time);
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