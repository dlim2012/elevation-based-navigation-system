
## REST API
```
[GET] Random Node: /api/v1/elena/random-node

[GET] Shortest Path: /api/v1/elena/shortest-path/random

[GET] Minimize elevation: /api/v1/elena/elena-minimize/random/max-length-ratio=1.5
(max-length-ratio configurable between 1 and 10)

[GET] Maximize elevation: /api/v1/elena/elena-maximize/random/max-length-ratio=1.5
(max-length-ratio configurable between 1 and 10)

[POST] Nearest Node: /api/v1/elena/nearest-node
Inputs: “lon” (longitude), “lat” (latitude)
Sample Input:
{
"lon": -71.0505,
"lat": 42.1017
}

[POST] Shortest Path:  /api/v1/elena/shortest-path
Inputs: “lon1” (longitude), “lat1” (latitude), “lon2” (longitude), “lat2” (latitude)

Sample Input:
{
"lon1": -72.518914,
"lat1": 42.375643,
"lon2": -72.520979,
"lat2": 42.377473
}

[GET] Minimize elevation: /api/v1/elena/elena-minimize
Inputs: “lon1” (longitude), “lat1” (latitude), “lon2” (longitude), “lat2” (latitude), “maxLengthRatio” (max Length Ratio)
Sample Input:
{
"lon1": -72.518914,
"lat1": 42.375643,
"lon2": -72.520979,
"lat2": 42.377473,
"maxLengthRatio": 1.5
}

[POST] Nearest Node: /api/v1/elena/elena-maximize
Inputs: “lon1” (longitude), “lat1” (latitude), “lon2” (longitude), “lat2” (latitude), “maxLengthRatio” (max Length Ratio), “numProduce”: number of child to attempt to produce in each generation, “numMaxSelect”: maximum number of instances surviving at the end of each generation, “numEpoch”: number of generations
Sample Input:
{
"lon1": -72.518914,
"lat1": 42.375643,
"lon2": -72.520979,
"lat2": 42.377473,
"maxLengthRatio": 1.5,
"numProduce": 30,
"numMaxSelect": 10,
"numEpoch": 10
}

```
## Outputs
```
* Node outputs:
“id” (node id), “lon” (longitude), “lat” (latitude), “elevation” (elevation), “restriction” (restrictions, if exists)

* Path outputs:
  “elevation” (path elevation), “length” (path length), “path” (list of “lat”, “lon”)
```