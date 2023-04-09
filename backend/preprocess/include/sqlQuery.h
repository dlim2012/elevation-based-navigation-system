#include <string>
#include <pqxx/pqxx>

using namespace std;
using namespace pqxx;

void timedExecution(connection& C, string query);
int getNumberOfRows(connection& C, string tableName);
int getNumberOfPoints(connection& C, string tableName);