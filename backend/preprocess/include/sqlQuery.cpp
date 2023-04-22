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
        cout << "failed: " << e.what();
        if (!skipFail) {
            throw e;
        }
        cout << "(skipping)" << endl;
    };
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
}


int timedExecutionReturnsInt(connection& C, string query, bool skipFail){
    auto t0 = chrono::high_resolution_clock::now();
    cout << query << endl;
    work w(C);
    result r;
    try {
        r = w.exec(query);
        w.commit();
    } catch(const std::exception &e) {
        cout << "failed: " << e.what() << endl;
        if (!skipFail)
            throw e;
    };
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return stoi(r[0][0].c_str());
}

long getNumberOfRows(connection& C, long lastID, string tableName){
    auto t0 = chrono::high_resolution_clock::now();
    work queryCount(C);
    string numRowsQuery;
    if (lastID > 0){
        numRowsQuery = "SELECT count(*) FROM " + tableName + " WHERE id > " + queryCount.esc(to_string(lastID)) + ";";
    } else {
        numRowsQuery = "SELECT count(*) FROM " + tableName + ";";
    }
    cout << numRowsQuery << endl;
    result r = queryCount.exec(numRowsQuery);
    int numRows = stol(r[0][0].c_str());
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return numRows;
}

long getNumberOfPoints(connection& C, long lastID, string tableName){
    auto t0 = chrono::high_resolution_clock::now();
    work queryCount(C);
    string numRowsQuery;
    if (lastID > 0){
        numRowsQuery = "SELECT COALESCE(sum(ST_NPOINTS(geom)), 0) FROM " + tableName + " WHERE id > " + queryCount.esc(to_string(lastID)) + ";";
    } else {
        numRowsQuery = "SELECT COALESCE(sum(ST_NPOINTS(geom)), 0) FROM " + tableName + ";";
    }
    cout << numRowsQuery << endl;
    result r = queryCount.exec(numRowsQuery);
    cout << r[0][0].c_str() << endl;
    int numRows = stol(r[0][0].c_str());
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return numRows;
}

long getLastId(connection& C, string tableName){
    auto t0 = chrono::high_resolution_clock::now();
    work queryCount(C);
    string numRowsQuery = "SELECT COALESCE(MAX(id), -1) FROM " + tableName + ";";
    cout << numRowsQuery << endl;
    result r = queryCount.exec(numRowsQuery);
    long numRows = stol(r[0][0].c_str());
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "time: " << t << " seconds" << endl << endl;
    return numRows;
}