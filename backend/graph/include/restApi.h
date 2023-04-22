#include <string>

#ifndef PATHFINDING_H
#define PATHFINDING_H

#include "pathfinding.h"

#endif

using namespace std;


crow::json::wvalue pathToGeoJson(
        Path *path, unordered_map<long, Geometry*> &umap);

crow::json::wvalue jsonFromNode(Node* node);




void runServer(string& connectionString, int port);