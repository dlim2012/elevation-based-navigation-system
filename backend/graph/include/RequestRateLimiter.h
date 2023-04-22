#include <deque>
#include <string>
#include <chrono>
#include <unordered_map>

using namespace std;

//const int MEASURE_TIME = 120000; // 2 minutes;
//const int ALLOWED_TIME = (int) (MEASURE_TIME * 3.0);

const int MEASURE_TIME_1 = 10000; // 10 seconds
const int ALLOWED_TIME_1 = (int) (MEASURE_TIME_1 * 12.0);

const int MEASURE_TIME_2 = 30000; // 30 seconds
const int ALLOWED_TIME_2 = (int) (MEASURE_TIME_2 * 10.0);

const int MEASURE_TIME_3 = 60000; // 30 seconds
const int ALLOWED_TIME_3 = (int) (MEASURE_TIME_3 * 8.0);

const int MEASURE_TIME_4 = 120000; // 120 seconds
const int ALLOWED_TIME_4 = (int) (MEASURE_TIME_4 * 6.0);

typedef chrono::time_point<std::chrono::high_resolution_clock> time_point;

int duration(time_point start, time_point end);

struct TimeInterval{
    time_point start;
    time_point end;
    int intervalInMilliseconds;
    TimeInterval(time_point start, time_point end);
};

class RequestRateLimiter {
public:
    unordered_map<string, deque<TimeInterval>> requests;
    int lockCount;
    mutex lock;

    RequestRateLimiter(int lockCount);
    void cleanRequests();
    bool isAllowed(const string& ipAddress, time_point& now);
    void addRecord(const string& ipAddress, time_point& start, time_point& end);
    bool acquireLock();
    void releaseLock();

};