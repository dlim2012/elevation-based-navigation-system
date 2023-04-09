#include <string>
#include <iostream>
#include <chrono>
#include <pqxx/pqxx>

using namespace std;
using namespace pqxx;

void timedExecution(connection& C, string query, bool skipFail){
    auto t0 = chrono::high_resolution_clock::now();
    cout << query << endl;
    work w(C);
    try {
        w.exec(query);
        w.commit();
    } catch(const std::exception &e) {
        cout << "failed: " << e.what() << endl;
        if (!skipFail)
            throw e;
    };
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
}

int getNumberOfRows(connection& C, string tableName){
    auto t0 = chrono::high_resolution_clock::now();
    const string numRowsQuery = "SELECT count(*) FROM " + tableName + ";";
    cout << numRowsQuery << endl;
    work queryCount(C);
    result r = queryCount.exec(numRowsQuery);
    int numRows = stoi(r[0][0].c_str());
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return numRows;
}

int getNumberOfPoints(connection& C, string tableName){
    auto t0 = chrono::high_resolution_clock::now();
    const string numRowsQuery = "SELECT sum(ST_NPOINTS(geom)) FROM " + tableName + ";";
    cout << numRowsQuery << endl;
    work queryCount(C);
    result r = queryCount.exec(numRowsQuery);
    int numRows = stoi(r[0][0].c_str());
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return numRows;
}
