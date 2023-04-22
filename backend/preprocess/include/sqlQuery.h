#include <string>
#include <pqxx/pqxx>

using namespace std;
using namespace pqxx;

void timedExecution(connection& C, string query, bool skipFail);
long getNumberOfRows(connection& C, long lastID,  string tableName);
long getNumberOfPoints(connection& C, long lastID,  string tableName);
long getLastId(connection& C, string tableName);