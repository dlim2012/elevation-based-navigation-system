


/* ===========================================================
            REST API - clients
============================================================== */

var host = "http://76.23.247.67:8080";
// var host = "http://127.0.0.1:18080";

// brasilia
// var center = [-15.8264, -47.9245]
// var zoom = 12

// mass
// var center = [42.434343, -71.957276]
// var zoom = 8

// northeast
// var center = [40.879454, -74.917818]

// us
var center = [38.237592, -95.209850]
var zoom = 10


const checkStatus = response => {
    console.log(response)
    if (response.ok) {
        return response;
    } else if (response.status === 503){
        if (markerNodeList.alertAvailable){
            alert("Server Busy.");
            markerNodeList.alertAvailable = false;
        }
    } else if (response.status === 429){
        if (markerNodeList.alertAvailable){
            alert("Too many requests.");
            markerNodeList.alertAvailable = false;
        }
    }

    const error = new Error(response.statusText);
    error.response = response;
    return Promise.reject(error);
}

export const getRandomNode = () => {
    const requestOptions = {
        method: 'GET'
    }
    return fetch(host + "/api/v1/elena/random-node", requestOptions).then(checkStatus);
}

export const getTwoNearNodes = (json) => {
    const requestOptions = {
        method: 'POST',
        body: JSON.stringify(json)
    }
    return fetch(host + "/api/v1/elena/two-near-nodes", requestOptions).then(checkStatus);
}

export const getShortestPath = (json) => {
    const requestOptions = {
        method: 'POST',
        body: JSON.stringify(json)
    }
    return fetch(host + "/api/v1/elena/shortest-path", requestOptions).then(checkStatus);
}

export const getElenaMinimizePath = (json) => {
    const requestOptions = {
        method: 'POST',
        body: JSON.stringify(json)
    }
    return fetch(host + "/api/v1/elena/elena-minimize", requestOptions).then(checkStatus);
}


export const getElenaMaximizePath = (json) => {
    const requestOptions = {
        method: 'POST',
        body: JSON.stringify(json)
    }
    return fetch(host + "/api/v1/elena/elena-maximize", requestOptions).then(checkStatus);
}



/* ===========================================================
            Global functions - map
============================================================== */


function showCoordinates (e) {
	alert(e.latlng);
}

function centerMap (e) {
	map.panTo(e.latlng);
}

function zoomIn (e) {
	map.zoomIn();
}

function zoomOut (e) {
	map.zoomOut();
}

var notAny = arr => arr.every(v => v === false);

function getDistance(latlng1, latlng2){
    var lat1_rad = latlng1.lat * Math.PI / 180;
    var lon1_rad = latlng1.lng * Math.PI / 180;
    var lat2_rad = latlng2.lat * Math.PI / 180;
    var lon2_rad = latlng2.lng * Math.PI / 180;
    return Math.acos(Math.sin(lat1_rad) * Math.sin(lat2_rad) +
        Math.cos(lat1_rad) * Math.cos(lat2_rad) * Math.cos(lon2_rad - lon1_rad)) * 6371
}

function onMarkerDragend(event) {
    var marker = event.target;
    var prevNode = marker.node.prevNode;

    if ((!useComputationLock || markerNodeList.lockCount === 0) && getAutoSearchOnMarkerDrag()){
        var graphType;
        var edgeBased;
        [graphType, edgeBased] = getGraphType();
        if (getComputationTime() <= 60
            && (
               (
                   marker.node.nextNode != null &&
                   !checkDistanceLimit(
                       getDistance(marker.getLatLng(), marker.node.nextNode.getLatLng()),
                       edgeBased
                   )
               )|| (
                    prevNode != null &&
                    !checkDistanceLimit(
                        getDistance(prevNode.getLatLng(), marker.getLatLng(),
                            edgeBased
                    ))
            ))){
            ckptView.updateCkptRows();
            alert ("The maximum distance exceeded. Please change the computation time criteria from the sidebar to increase maximum distance.")
            return;
        }
        marker.node.updatePathsNoLock(graphType, edgeBased, getPathTypes());
        if (prevNode != null) {
            prevNode.updatePathsNoLock(graphType, edgeBased, getPathTypes());
        }
    } else {
        if (prevNode != null) {
            prevNode.addToMapIfUpdated(getPathTypes())
        }
        marker.node.addToMapIfUpdated(getPathTypes())
    }
    updateAllView();
}

function checkDistanceLimit(distance, edgeBased){
    if (getComputationTime() <= 60 && getComputationTime() * 17 * (edgeBased + 1) < distance){
        return false;
    }
    return true;
}


/* ===========================================================
            Global functions - getters
============================================================== */

function getGraphType(){
    var checked = parseInt(document.querySelector('input[name="transportation_type"]:checked').value, 10);
    var graphType;
    var edgeBased;
    if (checked === 1){
        graphType = "hiking";
        edgeBased = 0;
    } else if (checked === 2){
        graphType = "cycling";
        edgeBased = 0;
    } else if (checked === 3){
        graphType = "cycling";
        edgeBased = 1;
    }
    // else if (checked === 4){
    //     graphType = "main_motorway";
    //     edgeBased = 0;
    // } else if (checked === 5){
    //     graphType = "main_motorway";
    //     edgeBased = 1;
    // }
    return [graphType, edgeBased];
}

function getPathTypes(){
    var pathTypes = [];
    if (document.getElementById('shortest-length').checked){
        pathTypes.push(0);
    }
    if (document.getElementById('minimize-elevation').checked){
        pathTypes.push(1);
    }
    if (document.getElementById('maximize-elevation').checked){
        pathTypes.push(2);
    }
    return pathTypes;
}

function getDuplicateEdgePolicy(){

    var checked = parseInt(
        document.querySelector('input[name="duplicate_edge_policy"]:checked').value, 10);
    var res;
    if (checked === 1){
        res = 0;
    } else if (checked === 2){
        res = 2;
    } else {
        res = 4;
    }
    return res;
}

function getComputationTime(){
    var checked = document.querySelector('input[name="computation_time"]:checked').value;
    return parseInt(checked);
}

function getAutoSearchOnMarkerDrag(){
    return document.getElementById('auto-search-on-marker-drag').checked;
}

function getAutoSearchOnMarkerRemove(){
    return document.getElementById('auto-search-on-marker-remove').checked;
}

function getUnifiedMaxLengthRatio(){
    return document.getElementById('unified-max-length-ratio').checked;
}

function getSearchPathAgain() {
    return document.getElementById('search-again').checked;
}

function copyDict(json){
    var newJson = {};
    for (var key in dict){
        newJson[key] = dict[key];
    }
    return newJson
}



/* ===========================================================
            Global variables
============================================================== */
var maxNumPaths = 4;
var numPathTypes = 3;
var allPathTypes = [0, 1, 2];
var colors = ['red', 'blue', 'green']
var tooltipMessages = [
    "Shortest path",
    "Minimize elevation gain",
    "Maximize elevation gain"
]
var markerPopup = {
    0: "Start",
    end: "End"
}
for (var index=1; index<maxNumPaths; index++){
    markerPopup[index] = `Checkpoint ${index}`
}
var decimalPlaces = 5;
var maxLengthRatioDecimalPlaces = 2;

var markerNodeList;
var ckptView;

var useComputationLock = false;

/* ===========================================================
            Map components - Model
============================================================== */


class Path {
    polyline;
    length;
    elevation;
    searchInfo;
    constructor(polyline, length, elevation, searchInfo){
        this.polyline = polyline;
        this.length = length;
        this.elevation = elevation;
        this.searchInfo = searchInfo
    }
    static pathFromCoordinates(coordinates, length, elevation, searchInfo,
                               color='', tooltipMessage=''){
        var polyline = L.polyline(
            coordinates,
            {color: color}
        )
        polyline.bindTooltip(L.tooltip()
            .setContent(
                tooltipMessage += " (length: " + length.toFixed(2) +
                " km, elevation gain: " +
                (elevation / 1000).toFixed(2) +
                " km)"
            )
        );
        return new Path(polyline, length, elevation, searchInfo);

    }
    static emptyPath(){
        var emptyJson = {
            "graph": "none",
            "lat1": 0.0,
            "lng1": 0.0,
            "lat2": 0.0,
            "lng2": 0.0,
            "maxLengthRatio": -1.0,
            "edgeBased": -1
        }
        return new Path(L.polyline([]), 0.0, 0.0, emptyJson)
    }
    addTo(map){
        this.polyline.addTo(map);
    }
    removeFrom(map){
        map.removeLayer(this.polyline);
    }
}

class MarkerNode {
    markerNodeList;
    id; // unique id
    marker; // marker added to map
    prevNode; // current prev node, null if not exists
    nextNode; // current next node, null if not exists
    maxLengthRatio;

    paths; // list of Path objects
    constructor(markerNodeList, id, marker, prevNode, nextNode, maxLengthRatio, numPathTypes){
        this.markerNodeList = markerNodeList;
        this.id = id;
        this.marker = marker;
        this.prevNode = prevNode;
        this.nextNode = nextNode;
        this.maxLengthRatio = maxLengthRatio;

        this.paths = [];
        for (let i=0; i<=numPathTypes; i++){
            this.paths.push(Path.emptyPath());
        }
    }

    compareNumber(val1, val2){
        return val1.toFixed(decimalPlaces) === val2.toFixed(decimalPlaces);
    }

    getMaxLengthRatio(){
        if (getUnifiedMaxLengthRatio()){
            return this.markerNodeList.startMarkerNode.maxLengthRatio;
        } else {
            return this.maxLengthRatio;
        }
    }

    compareJsonBase(json1, json2){
        for (var key in json1){
            if (typeof(json1[key]) === "string" && json1[key] !== json2[key]){
                return false;
            }
            if (typeof(json1[key]) === "number" && !this.compareNumber(json1[key], json2[key])){
                return false;
            }
        }
        return true;
    }

    compareLocation(json){
        var latlng1 = this.marker.getLatLng();
        if (this.nextNode == null){
            return false;
        }
        var latlng2 = this.nextNode.marker.getLatLng();
        if (!this.compareNumber(latlng1.lat, json["lat1"])){
            return false;
        }
        if (!this.compareNumber(latlng1.lng, json["lon1"])){
            return false;
        }
        if (!this.compareNumber(latlng2.lat, json["lat2"])){
            return false;
        }
        if (!this.compareNumber(latlng2.lng, json["lon2"])){
            return false;
        }
        // todo add criteria for edgeBased
        return true;
    }

    isUpdated(pathType, pathInfo){
        var graphType;
        var edgeBased;
        [graphType, edgeBased] = getGraphType();

        if (!this.compareLocation(pathInfo)
            || pathInfo["edgeBased"] !== edgeBased
            || pathInfo["graph"] !== graphType) {
            return false;
        }

        if (pathType === 0) {
            return true;
        }

        if (!this.compareNumber(this.getMaxLengthRatio(), pathInfo["maxLengthRatio"])) {
            return false;
        }

        if (pathType === 1) {
            return true;
        }

        if (getDuplicateEdgePolicy() !== pathInfo["duplicateEdges"]) {
            return false;
        }

        if (pathType === 2) {
            return true;
        }
    }

    getUpdated(pathTypes){
        var res = []
        for (var i=0; i<pathTypes.length; i++){
            var pathType = pathTypes[i];
            var pathInfo = this.paths[pathType].searchInfo
            res.push(this.isUpdated(pathType, pathInfo))
        }
        return res;
    }
    addToMapIfUpdated(pathTypes){
        var updated = this.getUpdated(pathTypes);
        var checkedPathTypes = getPathTypes();
        for (var i=0; i<pathTypes.length; i++){
            if (!checkedPathTypes.includes(pathTypes[i])){
                continue;
            }
            if (updated[i]) {
                this.paths[pathTypes[i]].addTo(this.markerNodeList.map);
            } else {
                this.paths[pathTypes[i]].removeFrom(this.markerNodeList.map);
            }
        }
    }

    saveAndAddToMapIfUpdated(pathType, json, data){
        if (!this.isUpdated(pathType, json)){
            return;
        }

        // save the calculated path
        this.paths[pathType] = Path.pathFromCoordinates(
            data.GeoJSON.coordinates, data.length, data.elevation, json,
            colors[pathType], tooltipMessages[pathType]
        )

        if (getPathTypes().includes(pathType)){
            return;
        }
        this.paths[pathType].addTo(this.markerNodeList.map);
    }

    shortestPathRequest(graphType, pathType, json){
        if (this.compareJsonBase(json, this.paths[pathType].searchInfo)){
            this.paths[pathType].addTo(this.markerNodeList.map);
            return;
        }
        this.paths[pathType].removeFrom(this.markerNodeList.map)
        this.lock();
        getShortestPath(json)
            .then(response => response.json())
            .then(data => {
                this.paths[pathType].removeFrom(this.markerNodeList.map)
                // Check if the node still exists in the list
                if (this.markerNodeList.getMarkerIndex(this) === -1){
                    return;
                }

                this.saveAndAddToMapIfUpdated(pathType, json, data)
                //
                // if (!this.isUpdated(pathType, json)){
                //     return;
                // }
                //
                // // save the calculated path
                // this.paths[pathType] = Path.pathFromCoordinates(
                //     data.GeoJSON.coordinates, data.length, data.elevation, json,
                //     colors[pathType], tooltipMessages[pathType]
                // )
                //
                //
                //
                // this.addToMapIfUpdated([pathType])
                console.log('shortest path: ', (new Date() - this.markerNodeList.startTime) /1000, " seconds");
            }).catch(error => {
            console.log(error);
        }).finally(() =>{
            this.unlock();
        });
    }
    elenaMinimizePathRequest(graphType, pathType, json){
        json["maxLengthRatio"] = this.getMaxLengthRatio();
        if (this.compareJsonBase(json, this.paths[pathType].searchInfo)){
            this.paths[pathType].addTo(this.markerNodeList.map);
            return;
        }
        this.paths[pathType].removeFrom(this.markerNodeList.map)
        this.lock();
        json["maxLength"] = -1;

        getElenaMinimizePath(json)
            .then(response => response.json())
            .then(data => {
                this.paths[pathType].removeFrom(this.markerNodeList.map)
                // Check if the node still exists in the list
                if (this.markerNodeList.getMarkerIndex(this) === -1){
                    return;
                }

                if (!this.isUpdated(pathType, json)){
                    return;
                }

                // save the calculated path
                this.paths[pathType] = Path.pathFromCoordinates(
                    data.GeoJSON.coordinates, data.length, data.elevation, json,
                    colors[pathType], tooltipMessages[pathType]
                )

                this.addToMapIfUpdated([pathType])
                console.log('elena-minimize: ', (new Date() - this.markerNodeList.startTime) /1000, " seconds");
            }).catch(error => {
                console.log(error);
            }).finally(() =>{
                this.unlock();
            });
    }
    elenaMaximizeRequest(graphType, pathType, json){
        json["maxLengthRatio"] = this.getMaxLengthRatio();
        json["duplicateEdges"] = getDuplicateEdgePolicy();
        json["seconds"] = getComputationTime();
        if (!getSearchPathAgain() && this.compareJsonBase(json, this.paths[pathType].searchInfo)){
            this.paths[pathType].addTo(this.markerNodeList.map);
            return;
        }
        this.paths[pathType].removeFrom(this.markerNodeList.map)
        this.lock();
        json["numProduce"] = 100;
        json["numMaxSelect"] = 30;
        json["numEpoch"] = 1000;
        json["maxLength"] = -1;
        getElenaMaximizePath(json)
            .then(response => response.json())
            .then(data => {
                this.paths[pathType].removeFrom(this.markerNodeList.map)
                // Check if the node still exists in the list
                if (this.markerNodeList.getMarkerIndex(this) === -1){
                    return;
                }

                if (!this.isUpdated(pathType, json)){
                    return;
                }

                // save the calculated path
                this.paths[pathType] = Path.pathFromCoordinates(
                    data.GeoJSON.coordinates, data.length, data.elevation, json,
                    colors[pathType], tooltipMessages[pathType]
                )

                this.addToMapIfUpdated([pathType])
                console.log('elena-maximize: ', (new Date() - this.markerNodeList.startTime) /1000, " seconds");
            }).catch(error => {
            console.log(error);
        }).finally(() =>{
            this.unlock();
        });

    }

    lock(){
        this.markerNodeList.lockCount += 1;
        ckptView.setLoading(this);
        searchButton.state('loading');
    }

    unlock(){
        this.markerNodeList.lockCount -= 1;
        ckptView.unsetLoading(this);
        if (this.markerNodeList.lockCount === 0){
            searchButton.state('search-route');
            this.markerNodeList.alertAvailable = true;
        }
        updateAllView()
    }

    updatePathsNoLock(graphType, edgeBased, pathTypes, recalculate=0){
        if (this.size < 2){
            return;
        }
        this.markerNodeList.startTime = new Date()

        var latlng1 = this.marker.getLatLng();
        if (this.nextNode == null){
            return;
        }
        var latlng2 = this.nextNode.marker.getLatLng();
        console.log("distance:", getDistance(this.marker.getLatLng(), this.nextNode.marker.getLatLng()), "km")

        for (var pathType of pathTypes){
            var json = {
                graph: graphType,
                lat1 : latlng1.lat,
                lon1 : latlng1.lng,
                lat2 : latlng2.lat,
                lon2 : latlng2.lng,
                edgeBased: edgeBased
            }

            if (pathType === 0){
                this.shortestPathRequest(graphType, pathType, json)
            } else if (pathType === 1){
                this.elenaMinimizePathRequest(graphType, pathType, json);
            } else if (pathType === 2){
                this.elenaMaximizeRequest(graphType, pathType, json, recalculate > 0);
            }
        }
    }

    updatePaths(graphType, edgeBased, pathTypes, recalculate=0){
        if (useComputationLock && this.markerNodeList.lockCount !== 0){
            return;
        }
        this.updatePathsNoLock(graphType, edgeBased, pathTypes, recalculate);
    }

    updateRelevantPathsNoLock(){
        var graphType;
        var edgeBased;
        [graphType, edgeBased] = getGraphType();
        var pathTypes = getPathTypes();
        if (this.nextNode != null){
            this.updatePathsNoLock(graphType, edgeBased, pathTypes);
        }
        if (this.prevNode != null){
            this.prevNode.updatePathsNoLock(graphType, edgeBased, pathTypes)
        }
    }
    updateRelevantPaths(){
        if (useComputationLock && this.markerNodeList.lockCount !== 0){
            return;
        }
        this.updateRelevantPathsNoLock();
    }
    togglePathType(pathType, show, update=false){
        if (show) {
            var graphType;
            var edgeBased;
            [graphType, edgeBased] = getGraphType();
            if (update) {
                this.updatePathsNoLock(graphType, edgeBased, [pathType])
            } else {
                this.paths[pathType].addTo(this.markerNodeList.map);
            }
        } else {
            this.paths[pathType].removeFrom(this.markerNodeList.map);
        }
    }
    removeAllPaths(){
        for (let i=0; i<numPathTypes; i++){
            this.paths[i].removeFrom(this.markerNodeList.map);
        }
    }
    getLatLng(){
        return this.marker.getLatLng();
    }
    setIcon(icon){
        this.marker.setIcon(icon);
    }
    setPrevNode(prevNode){
        this.prevNode = prevNode;
    }
    setNextNode(nextNode){
        this.nextNode = nextNode;
    }
    setLatLng(latlng){
        this.marker.setLatLng(latlng);
    }
}

class MarkerNodeList{
    map;
    size; // number of MarkerNode including the endNode
    startMarkerNode; // start marker
    endMarkerNode;
    numPathTypes;

    startTime;

    lastMarkerNodeId;
    lockCount;
    alertAvailable;
    constructor(map, startLatLng, endLatLng, numPathTypes){
        this.map = map;
        this.lastMarkerNodeId = 0;
        this.numPathTypes = numPathTypes;
        this.endMarkerNode = null;
        this.lockCount = 0;
        this.alertAvailable = true;

        this.startTime = new Date();

        // make end Marker Node
        this.size = 2;
        this.endMarkerNode = new MarkerNode(
            this,
            this.newMarkerNodeId(),
            this.newMarker(1, endLatLng).addTo(this.map),
            null,
            null,
            1.5,
            numPathTypes
        )

        // add start node
        // this.startMarkerNode = this.addNewMarkerNodeSecondLast(startLatLng, false);
        this.startMarkerNode = new MarkerNode(
            this,
            this.newMarkerNodeId(),
            this.newMarker(0, startLatLng),
            null,
            this.endMarkerNode,
            1.5,
            this.numPathTypes

        )
        this.endMarkerNode.setPrevNode(this.startMarkerNode);
        this.startMarkerNode.marker.addTo(this.map);
        this.endMarkerNode.marker.addTo(this.map);

        this.startMarkerNode.marker.node = this.startMarkerNode;
        this.endMarkerNode.marker.node = this.endMarkerNode;

        this.startMarkerNode.marker.bindPopup(markerPopup[0]);
        this.endMarkerNode.marker.bindPopup(markerPopup['end']);
    }

    addNewMarkerNodeSecondLast(latlng){
        if (this.size === 0 || this.endMarkerNode == null){
            console.error("EndMarkerNode does not exist.")
        }
        this.size += 1;
        var markerId = this.newMarkerNodeId();
        var prevNode;
        if (this.size === 1){
            prevNode = null;
        } else {
            prevNode = this.endMarkerNode.prevNode;
        }
        var index = this.size - 2;
        var marker = this.newMarker(index, latlng);

        var markerNode = new MarkerNode(
            this,
            markerId,
            marker,
            prevNode,
            this.endMarkerNode,
            this.startMarkerNode.maxLengthRatio,
            this.numPathTypes
        );
        this.endMarkerNode.prevNode = markerNode;

        if (prevNode != null) {
            prevNode.setNextNode(markerNode);
            prevNode.addToMapIfUpdated(getPathTypes())
        }

        marker.addTo(this.map);
        marker.bindPopup(markerPopup[index])
        marker.node = markerNode;

        updateAllView()
        // ckptView.rows[index + 1].children[1].children[0].value = null;
        // ckptView.rows[index + 1].children[2].children[0].value = null;
        // ckptView.rows[index + 1].children[3].children[0].value = null;

        return markerNode;
    }
    newMarkerNodeId(){
        if (this.lastMarkerNodeId === Number.MAX_VALUE){
            console.error("too big node Id");
        }
        this.lastMarkerNodeId += 1;
        return this.lastMarkerNodeId;
    }
    newMarkerIcon(index){
        var icon;
        if (index === this.size - 1){
            // end icon
            icon =  L.AwesomeMarkers.icon({
                icon: 'flag',
                prefix: 'fa',
                markerColor: 'blue'
              });
        } else if (index === 0){
            // starter icon
            icon = L.AwesomeMarkers.icon({
                icon: 'location-arrow',
                prefix: 'fa',
                markerColor: 'blue'
              });
        } else {
            // numbered icon
            icon =  L.AwesomeMarkers.icon({
                icon: '',
                prefix: 'fa',
                markerColor: 'blue',
                html: index
              });
        }
        return icon;
    }
    newMarker(index, latlng){
        var marker = L.marker(
            latlng,
            {
                icon: this.newMarkerIcon(index),
                draggable: 'true',
                contextmenu: true,
                contextmenuItems: [{
                    text: 'Remove this marker',
                    callback: (e) => {
                        removeMarkerNode(markerNodeList.getMarkerIndex(e.relatedTarget.node))
                    }
                }]
            }
        );
        marker.on('dragend', onMarkerDragend)
        return marker;
    }
    extractMarkerNode(index, remove=true){
        if (index > this.size){
            console.error("error")
        }
        var markerNode = this.getMarkerNode(index);
        markerNode.removeAllPaths();
        if (markerNode.prevNode != null){
            markerNode.prevNode.nextNode = markerNode.nextNode;
            markerNode.prevNode.removeAllPaths();
        } else {
            this.startMarkerNode = markerNode.nextNode;
        }
        if (markerNode.nextNode != null){
            markerNode.nextNode.prevNode = markerNode.prevNode;
        } else {
            this.endMarkerNode = markerNode.prevNode;
        }

        if (remove) {
            this.map.removeLayer(markerNode.marker);
        }
        this.size -= 1;
        return markerNode;

    }
    insertMarkerNode(markerNode, index, add){
        if (index > this.size){
            console.error("error")
        }
        if (this.size === 0){
            markerNode.prevNode = null;
            markerNode.nextNode = null;
            this.startMarkerNode = markerNode;
            this.endMarkerNode = markerNode;
            return;
        }

        if (index === 0) {
            markerNode.prevNode = null;
            markerNode.nextNode = this.startMarkerNode;
            this.startMarkerNode.prevNode = markerNode;
            this.startMarkerNode = markerNode;
        } else if (index === this.size){
            markerNode.prevNode = this.endMarkerNode;
            markerNode.nextNode = null;
            this.endMarkerNode.nextNode = markerNode;
            this.endMarkerNode = markerNode;
        } else {
            var nextNode = this.getMarkerNode(index);
            var prevNode = nextNode.prevNode
            markerNode.nextNode = nextNode;
            markerNode.prevNode = prevNode;
            nextNode.prevNode = markerNode;
            prevNode.nextNode = markerNode;
        }

        if (add){
            markerNode.marker.addTo(this.map)
        }
        this.size += 1;
    }
    removeMarkerNode(index){
        this.extractMarkerNode(index, true);
    }
    moveMarkerNode(fromIndex, toIndex){
        if (fromIndex === toIndex){
            return;
        }
        var markerNode = this.extractMarkerNode(fromIndex, false);
        this.insertMarkerNode(markerNode, toIndex, false);
    }
    switchMarkerNodes(index1, index2){
        if (index1 === index2){
            return;
        }
        if (index1 > index2){
            [index1, index2]  = [index2, index1]
        }

        var markerNode1 = this.getMarkerNode(index1);
        var markerNode2 = this.getMarkerNode(index2);

        markerNode1.removeAllPaths();
        markerNode2.removeAllPaths();

        [markerNode1.prevNode, markerNode2.prevNode] = [markerNode2.prevNode, markerNode1.prevNode]
        if (markerNode1.prevNode == null){
            this.startMarkerNode = markerNode1;
        } else {
            markerNode1.prevNode.removeAllPaths();
            markerNode1.prevNode.nextNode = markerNode1;
        }
        if (markerNode2.prevNode == null){
            this.startMarkerNode = markerNode2;
        } else {
            markerNode2.prevNode.removeAllPaths();
            markerNode2.prevNode.nextNode = markerNode2;
        }

        [markerNode1.nextNode, markerNode2.nextNode] = [markerNode2.nextNode, markerNode1.nextNode]
        if (markerNode1.nextNode == null){
            this.endMarkerNode = markerNode1;
        } else {
            markerNode1.nextNode.prevNode = markerNode1;
        }
        if (markerNode2.nextNode == null){
            this.endMarkerNode = markerNode2;
        } else {
            markerNode2.nextNode.prevNode = markerNode2;
        }
    }

    printAll(){
        var node = this.startMarkerNode;
        while (node != null){
            console.log('print ', node.id, node.marker.getLatLng())
            node = node.nextNode;
        }
        console.log('print done')
    }

    getSize(){
        return this.size;
    }
    getMarkerNode(index){
        if (index >= this.size){
            console.error("index error")
        }
        var markerNode = this.startMarkerNode;
        for(let i=0; i<index; i++){
            markerNode = markerNode.nextNode;
        }
        return markerNode;
    }
    getMarkerIndex(markerNode){
        var node = this.startMarkerNode;
        let index = 0;
        while (node != null){
            if (node == markerNode){
                return index;
            }
            index += 1;
            node = node.nextNode;
        }
        return -1;
    }
    updateMarkerIconsPopupsAndContextMenu(){
        var markerNode = this.startMarkerNode;
        var index = 0;
        while (markerNode != null){
            markerNode.setIcon(this.newMarkerIcon(index));
            markerNode.marker.bindPopup(markerPopup[index]);
            markerNode.marker.options.contextmenuItems[0].disabled = (this.getSize() <= 2);
            markerNode = markerNode.nextNode;
            index += 1;
        }
        this.endMarkerNode.marker.bindPopup(markerPopup["end"])
    }
    updateMaxLengthRatio(maxLengthRatio){
        var markerNode = this.startMarkerNode;
        while (markerNode != null) {
            markerNode.maxLengthRatio = maxLengthRatio;
            markerNode = markerNode.nextNode;
        }
    }
    updateAllPaths(){
        var graphType;
        var edgeBased;
        [graphType, edgeBased] = getGraphType();

        // Check max distance
        var maxDistance = 0
        var markerNode = this.startMarkerNode;
        for(let i=0; i<this.size-1; i++){
            var distance = getDistance(markerNode.getLatLng(),
                markerNode.nextNode.getLatLng());
            if (maxDistance < distance){
                maxDistance = distance;
            }
            markerNode = markerNode.nextNode;
        }
        if (!checkDistanceLimit(maxDistance, edgeBased)){
            alert ("The maximum distance exceeded. Please change the computation time criteria from the sidebar to increase maximum distance.")
            return;
        }

        if (!useComputationLock || this.lockCount !== 0) {
            var markerNode = this.startMarkerNode;
            for (let i = 0; i < this.size - 1; i++) {
                markerNode.updatePathsNoLock(graphType, edgeBased, getPathTypes());
                markerNode = markerNode.nextNode;
            }
        }
    }
    togglePathType(pathType, show, update=false){
        var node = this.startMarkerNode;
        while (node != null){
            node.togglePathType(pathType, show, update)
            node = node.nextNode
        }
    }
    updateAllPathView(){
        var pathTypes = getPathTypes();
        var node = this.startMarkerNode;
        while (node != null){
            node.addToMapIfUpdated(pathTypes)
            node = node.nextNode
        }
    }
    getTotalLengthAndElevationGain(){
        var pathTypes = getPathTypes();
        var pathTypesUpdated = [];
        var lengths = [];
        var elevationGains = [];
        for (var i=0; i<pathTypes.length; i++){
            pathTypesUpdated.push(0);
            lengths.push(0);
            elevationGains.push(0);
        }
        var node = this.startMarkerNode;
        while (node != null){
            var updated = node.getUpdated(pathTypes);
            for (var i=0; i<pathTypes.length; i++){
                if (updated[i]){
                    console.log(node.paths[pathTypes[i]])
                    pathTypesUpdated[i] += 1;
                    lengths[i] += node.paths[pathTypes[i]].length;
                    elevationGains[i] += node.paths[pathTypes[i]].elevation;
                }
            }
            node = node.nextNode;
        }
        return [pathTypes, lengths, elevationGains, pathTypesUpdated];
    }
}


/* ===========================================================
            Map components - View - Checkpoint Table
============================================================== */

class ckptTable{
    markerNodeList;
    rows;
    rowIdToMarkerNodeIndex;

    constructor(markerNodeList) {
        this.markerNodeList = markerNodeList;
        this.rows = document.getElementById("ckptTable").rows;
        this.rowIdToMarkerNodeIndex = {};


        this.updateCkptRows();

        for (let rowIndex=1; rowIndex<maxNumPaths+3; rowIndex++){
            if (rowIndex === maxNumPaths + 1){
                continue;
            }
            // listener: locate
            this.rows[rowIndex].children[0].children[0].addEventListener(
                "click",
                onRowLocateButton
            )

            // listener: on input changes
            this.rows[rowIndex].addEventListener("input", onCkptTableRowInputChange)
            this.rows[rowIndex].addEventListener("propertychange", onCkptTableRowInputChange)

            // listener: search path
            this.rows[rowIndex].children[4].children[0].addEventListener(
                "click",
                (e) => {
                    var markerNode = this.markerNodeList
                        .getMarkerNode(this.rowIdToMarkerNodeIndex[e.target.parentNode.parentNode.parentNode.id])

                    var graphType;
                    var edgeBased;
                    [graphType, edgeBased] = getGraphType();
                    markerNode.updatePaths(graphType, edgeBased, getPathTypes());
                }
            )

            // listener: dragging
            this.rows[rowIndex].addEventListener("dragstart", onCkptTableDragStart);
            this.rows[rowIndex].addEventListener("dragover", onCkptTableDragover);
            this.rows[rowIndex].addEventListener("dragend", onCkptTableDragend);


        }
        // listener: search all
        document.getElementById("search-all")
            .addEventListener(
                "click",
                (e) => {markerNodeList.updateAllPaths()}
            );

        // listener: remove nodes
        for (let i=0; i<maxNumPaths+1; i++){
            document.getElementById(`ckpt-${i}-remove`)
                .addEventListener(
                    "click",
                    onRemoveMarkerNode
                );
        }

        // add checkpoint
        document.getElementById('add-ckpt-row')
            .addEventListener(
                "click",
                (e) => this.addCheckpoint(map.getCenter(), false));


        // set initial display of checkpoints to 'none'
        for (let i=1; i<maxNumPaths; i++){
            this.rows[i+1].style.display = 'none';
        }
    }

    setCkptTableRow(markerNode, index){
        var row;
        if (index === 0){
            row = this.rows[1];
            row.children[0].children[0].innerHTML = 'Start'
            row.children[3].children[0].disabled = false;
        }
        else if (index === this.markerNodeList.getSize()-1){
            row = this.rows[maxNumPaths + 2];
            row.children[0].children[0].innerHTML = 'End'
            row.children[3].children[0].disabled = true;
        } else {
            row = this.rows[index+1];
            row.children[0].children[0].innerHTML = index
            if (getUnifiedMaxLengthRatio()){
                row.children[3].children[0].disabled = true;
            } else {
                row.children[3].children[0].disabled = false;
            }
        }


        row.children[4].children[0].disabled = false;
        row.children[4].children[1].disabled = false;

        row.children[1].children[0].value = markerNode.marker.getLatLng().lat.toFixed(decimalPlaces);
        row.children[2].children[0].value = markerNode.marker.getLatLng().lng.toFixed(decimalPlaces);

        if (getUnifiedMaxLengthRatio()) {
            if (index === 0){
                row.children[3].children[0].value = parseFloat(this.markerNodeList.startMarkerNode.maxLengthRatio).toFixed(maxLengthRatioDecimalPlaces);
            } else {
                row.children[3].children[0].value = null;
            }
        } else {
            row.children[3].children[0].value = parseFloat(markerNode.maxLengthRatio).toFixed(maxLengthRatioDecimalPlaces);
        }

        row.style.display = '';
    }

    updateCkptRows (){
        var markerNode = this.markerNodeList.startMarkerNode;
        var index = 0;
        while (index<maxNumPaths){
            this.rowIdToMarkerNodeIndex[this.rows[index+1].id] = index;
            this.setCkptTableRow(markerNode, index)
            markerNode = markerNode.nextNode;
            this.rows[index+1].style.display = '';
            index += 1;
            if (markerNode == null){
                break;
            }

            // set last node
            if (markerNode.nextNode == null){
                this.rowIdToMarkerNodeIndex[this.rows[maxNumPaths+2].id] = this.markerNodeList.getSize()-1;
                this.setCkptTableRow(markerNode, this.markerNodeList.getSize()-1)
                break;
            }
        }

        // disable display for removed rows
        if (index < maxNumPaths) {
            while (index < maxNumPaths) {
                this.rows[index + 1].style.display = 'none';
                index += 1;
            }
            this.rows[maxNumPaths+1].style.display = '';
        }

        // disable remove button if number of checkpoints is less than 3
        if (this.markerNodeList.getSize() < 3) {
            for (var rowIndex = 1; rowIndex < maxNumPaths + 3; rowIndex++) {
                if (rowIndex === maxNumPaths + 1) {
                    continue;
                }
                this.rows[rowIndex].children[4].children[1].disabled = true;
            }
        }
    }


    addCheckpoint(latlng, calculate=false){
        if (this.markerNodeList.getSize() === maxNumPaths + 1){
            return;
        }
        this.rows[this.markerNodeList.getSize()].style.display = '';
        var markerNode = this.markerNodeList.addNewMarkerNodeSecondLast(latlng);
        if (this.markerNodeList.getSize() === maxNumPaths + 1){
            this.rows[maxNumPaths+1].style.display = 'none';
        }
        if (calculate){
            markerNode.updateRelevantPaths()
        }
    }

    setLoading(markerNode){
        // todo
        // var index = this.markerNodeList.getMarkerIndex(markerNode);
        // if (index === -1 || index === this.markerNodeList.getSize() - 1){
        //     return;
        // }


        // this.rows[index+1].children[4].children[0].innerHTML = '<i className="fa-solid fa-spinner"></i>'

    }

    unsetLoading(markerNode){
        // todo
        // var index = this.markerNodeList.getMarkerIndex(markerNode);
        // if (index === -1 || index === this.markerNodeList.getSize() - 1){
        //     return;
        // }
        // this.rows[index+1].children[4].children[0].innerHTML = '<i className="fa-solid fa-route"></i>'
    }
}

function onRowLocateButton(e){
    var markerNode = markerNodeList
        .getMarkerNode(ckptView.rowIdToMarkerNodeIndex[e.target.parentNode.parentNode.id])
    markerNode.marker.openPopup();
    map.panTo(markerNode.getLatLng());
}

function onCkptTableRowInputChange(e){
    var markerIndex = ckptView.rowIdToMarkerNodeIndex[e.target.parentNode.parentNode.id]
    var markerNode = ckptView.markerNodeList.getMarkerNode(markerIndex);

    var key = e.target.id.slice(e.target.id.lastIndexOf('-')+1);
    var value = e.target.value;

    // Update latitude and longitude
    if (key === 'lat'){
        if (value <= -90 || value >= 90){
            return;
        }
        markerNode.setLatLng(L.latLng(value, markerNode.getLatLng().lng))
    } else if (key === 'lng'){
        if (value < -180 || value > 180){
            return;
        }
        markerNode.setLatLng(L.latLng(markerNode.getLatLng().lat, value))
    } else if (key === 'r'){
        if (getUnifiedMaxLengthRatio()){
            markerNodeList.updateMaxLengthRatio(parseFloat(value));
        } else {
            markerNode.maxLengthRatio = parseFloat(value);
        }
    }
    updateAllView()
}


var fromIndex;
var toIndex;
var dragSource;
var dragTo;

function rowSwap(i1, i2){
    if (i1 !== -1 && i2 !== -1 && i1 !== i2) {
        if (i1 === maxNumPaths || i2 === maxNumPaths) {
            return;
        }

        if (i1 < maxNumPaths && i2 < maxNumPaths) {
            if (i2 > i1)
                dragTo.after(dragSource);
            else
                dragTo.before(dragSource);
        } else if (i2 === maxNumPaths + 1) {
            dragTo.after(dragSource)
            ckptView.rows[markerNodeList.getSize() - 1].before(dragTo);
        } else if (i1 === maxNumPaths + 1) {
            dragSource.after(dragTo);
            ckptView.rows[markerNodeList.getSize() - 1].before(dragSource);

        }
    }
}
function onCkptTableDragStart (e){
    dragSource = e.target;
    fromIndex = ckptView.rowIdToMarkerNodeIndex[dragSource.id]
}

function onCkptTableDragover(e){
    e.preventDefault();
    if (e.target.nodeName === 'TD'){
        dragTo = e.target.parentNode;
    } else if (e.target.nodeName === 'INPUT'){
        dragTo = e.target.parentNode.parentNode;
    }

    var children=Array.from(dragTo.parentNode.children);

    var i1 = children.indexOf(dragSource);
    var i2 = children.indexOf(dragTo);

    var setting = document.querySelector('input[name="ckpt-table-drag-policy"]:checked').value;
    if (setting === "2"){
        // record the swap and swap at dragend
        fromIndex = i1;
        toIndex = i2;
    } else if (setting === "1") {
        // swap now
        rowSwap(i1, i2);
    }
}

function onCkptTableDragend(e){
    var setting = document.querySelector('input[name="ckpt-table-drag-policy"]:checked').value;
    if (setting === "2"){
        rowSwap(fromIndex, toIndex)
        if (fromIndex === maxNumPaths + 1){
            fromIndex = markerNodeList.getSize() - 1;
        }
        if (toIndex === maxNumPaths + 1){
            toIndex = markerNodeList.getSize() - 1;
        }
        markerNodeList.switchMarkerNodes(fromIndex, toIndex);

    } else if (setting === "1") {
        if (fromIndex === maxNumPaths){
            fromIndex = markerNodeList.getSize() - 1;
        }
        toIndex = -1;
        for (let i=0; i<maxNumPaths; i++){
            if (dragSource.id === ckptView.rows[i+1].id){
                toIndex = i;
                break;
            }
        }
        if (toIndex === -1){
            toIndex = markerNodeList.getSize() - 1;
        }
        markerNodeList.moveMarkerNode(fromIndex, toIndex);
    }
    console.log(fromIndex, toIndex)
    if (getAutoSearchOnMarkerDrag()){
        if (!useComputationLock || markerNodeList.lockCount === 0) {
            markerNodeList.getMarkerNode(fromIndex).updateRelevantPathsNoLock();
            markerNodeList.getMarkerNode(toIndex).updateRelevantPathsNoLock()
        }
    }
    markerNodeList.updateMarkerIconsPopupsAndContextMenu();
    ckptView.updateCkptRows(markerNodeList);

}


function removeMarkerNode(index){
    console.log(index)
    markerNodeList.removeMarkerNode(index);
    markerNodeList.updateMarkerIconsPopupsAndContextMenu();
    if (getAutoSearchOnMarkerRemove()){
        if (index >= 1){
            var graphType;
            var edgeBased;
            [graphType, edgeBased] = getGraphType()
            markerNodeList.getMarkerNode(index-1).updatePaths(
                graphType, edgeBased, getPathTypes()
            )
        }
    }
    ckptView.updateCkptRows();

}
function onRemoveMarkerNode(e){
    if (markerNodeList.getSize() < 3){
        return;
    }
    var index = ckptView.rowIdToMarkerNodeIndex[e.target.parentNode.parentNode.parentNode.id]
    removeMarkerNode(index)
}

function updateTotalLengthAndElevationGain(){
    var pathTypes;
    var lengths;
    var elevationGains;
    var pathTypesUpdated;
    [pathTypes, lengths, elevationGains, pathTypesUpdated] = markerNodeList.getTotalLengthAndElevationGain()

    var output = "<p>"
    for (var i=0; i<pathTypes.length; i++){
        var pathType = pathTypes[i];
        output += "*" + tooltipMessages[pathType] + "<br>";
        output += "&emsp; Number of paths:  &emsp;" + pathTypesUpdated[i] + "<br>"
        output += "&emsp; Total length: &nbsp; &emsp;&emsp; &emsp;" + lengths[i].toFixed(2) +
            " km <br>";
        output += "&emsp; Elevation Gain: &emsp; &emsp;" + (elevationGains[i]/1000).toFixed(2) +
            " km<br>";
        if (i < pathTypes.length - 1){
            output += "<br>"
        }
    }
    output += "</p>";
    document.getElementById("results").innerHTML = output;
}

function updateContextMenu(){
    // Add checkpoint option
    if (markerNodeList.getSize() === maxNumPaths + 1){
        map.contextmenu.setDisabled(4, true)
    } else {
        map.contextmenu.setDisabled(4, false)
    }
    
    // Remove checkpoint option
    if (markerNodeList.getSize() === maxNumPaths + 1){
        
    }
}

function updateAllView(){
    ckptView.updateCkptRows();
    updateTotalLengthAndElevationGain();
    updateContextMenu();

    markerNodeList.updateMarkerIconsPopupsAndContextMenu();
    markerNodeList.updateAllPathView();
}

/* ===========================================================
            Listener: settings
============================================================== */

document.getElementById("shortest-length").addEventListener(
    "click",
    () => {
        var update = !useComputationLock || markerNodeList.lockCount === 0;
        markerNodeList.togglePathType(
            0, document.getElementById('shortest-length').checked, update)
        updateAllView();
    }, false);

document.getElementById("minimize-elevation").addEventListener(
    "click",
     () => {
         var update = !useComputationLock || markerNodeList.lockCount === 0;
         markerNodeList.togglePathType(
             1, document.getElementById('minimize-elevation').checked, update)
         updateAllView();
    },
     false);
document.getElementById("maximize-elevation").addEventListener(
    "click",
     () => {
         var update = !useComputationLock || markerNodeList.lockCount === 0;
         markerNodeList.togglePathType(
             2, document.getElementById('maximize-elevation').checked, update)
         updateAllView();
    },
    false);
document.getElementById("unified-max-length-ratio").addEventListener(
    "click",
    ()=> {
        ckptView.updateCkptRows();
    }
)

/* ===========================================================
            Initialize
============================================================== */


/* ===========================================================
            Functions
============================================================== */


function getArrows(arrLatlngs, color, arrowCount, mapObj) {

    if (typeof arrLatlngs === undefined || arrLatlngs == null ||
(!arrLatlngs.length) || arrLatlngs.length < 2)
    return [];

    if (typeof arrowCount === 'undefined' || arrowCount == null)
        arrowCount = 1;

    if (typeof color === 'undefined' || color == null)
        color = '';
    else
        color = 'color:' + color;

    var result = [];
    var interval = 10.0;
    var remain = interval;
    var curDist, curRemain;
    var prev, cur = arrLatlngs[0];
    for (var i = 1; i < arrLatlngs.length; i++) {
        var icon = L.divIcon({ className: 'arrow-icon', bgPos: [5, 5], html: '<div style="' + color + ';transform: rotate(' + getAngle(arrLatlngs[i - 1], arrLatlngs[i], -1).toString() + 'deg)">▶</div>' });
        // for (var c = 1; c <= arrowCount; c++) {
        //     result.push(L.marker(myMidPoint(arrLatlngs[i], arrLatlngs[i - 1], (c / (arrowCount + 1)), mapObj), { icon: icon }));
        // }
        prev = cur;
        cur = arrLatlngs[i];
        curDist = getDistance(prev, cur);
        if (curDist < remain){
            remain -= curDist;
            continue;
        }
        curRemain = curDist - remain;
        result.push(L.marker(pointBetween(prev, cur, curRemain, curDist), {icon: icon}));
        while (curRemain > interval){
            curRemain -= interval;
            result.push(L.marker(pointBetween(prev, cur, curRemain, curDist), {icon: icon}));
        }
        remain = interval - curRemain;
    }
    return result;
}

function getAngle(latLng1, latlng2, coef) {
    var dy = latlng2[0] - latLng1[0];
    var dx = Math.cos(Math.PI / 180 * latLng1[0]) * (latlng2[1] - latLng1[1]);
    var ang = ((Math.atan2(dy, dx) / Math.PI) * 180 * coef);
    return (ang).toFixed(2);
}

function pointBetween(prev, cur, curDist, dist){
    return [cur[0] - (cur[0] - prev[0]) * (curDist / dist), cur[1] - (cur[1] - prev[1]) * (curDist / dist)];
}




/* ===========================================================
            Map components - main
============================================================== */


var map = L.map('map', {
    contextmenu: true,
    contextmenuWidth: 140,
    zoomControl: false,
    contextmenuItems: [
    {
	    text: 'Center map here',
	    callback: centerMap
	},
    // '-',
    {
	    text: 'Zoom in',
	    // icon: 'images/zoom-in.png',
	    callback: zoomIn
	},
    {
	    text: 'Zoom out',
	    // icon: 'fas fa-location-dot',
	    callback: zoomOut
	},
    '-',
    {
        text: 'Add checkpoint here',
        callback: (event) => {
            ckptView.addCheckpoint(event.latlng, getAutoSearchOnMarkerDrag());
            updateAllView();
        },
    },
    {
        text: 'Set start marker here',
        // icon: 'fas fa-location-pen',
        callback: (event) => {
            markerNodeList.startMarkerNode.setLatLng(event.latlng)
            updateAllView();
            if (getAutoSearchOnMarkerDrag()){
                var graphType;
                var edgeBased;
                [graphType, edgeBased] = getGraphType();
                markerNodeList.startMarkerNode.updatePaths(graphType, edgeBased, getPathTypes())
            }
        }
    },
    {
        text: 'Set end marker here',
        callback: (event) => {
            markerNodeList.endMarkerNode.setLatLng(event.latlng)
            updateAllView();
            if (getAutoSearchOnMarkerDrag()){
                var graphType;
                var edgeBased;
                [graphType, edgeBased] = getGraphType();
                markerNodeList.endMarkerNode.prevNode.updatePaths(graphType, edgeBased, getPathTypes());
            }
        }
    },
    // {
    //     text: 'Search paths',
    //     callback: calculatePaths
    // }
    ]
}).setView(center, zoom);



var osm = L.tileLayer(
    'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
    {
        maxZoom: 18,
        attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
    }
);
osm.addTo(map);



// var OpenTopoMap = L.tileLayer(
//     'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png',
//     {
//         maxZoom: 17,
//         attribution: 'Map data: &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, <a href="http://viewfinderpanoramas.org">SRTM</a> | Map style: &copy; <a href="https://opentopomap.org">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>)'
//
//     }
// );
// OpenTopoMap.addTo(map);

var CartoDB_DarkMatter = L.tileLayer(
    'https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
    {
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
        subdomains: 'abcd',
        maxZoom: 20
    }
);
CartoDB_DarkMatter.addTo(map);

var CartoDB_Voyager = L.tileLayer(
    'https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png',
    {
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
        subdomains: 'abcd',
        maxZoom: 20
    }
);
CartoDB_Voyager.addTo(map)

var googleStreets = L.tileLayer(
    'http://{s}.google.com/vt?lyrs=m&x={x}&y={y}&z={z}',
    {
        maxZoom: 20,
        subdomains:['mt0','mt1','mt2','mt3']
    }
);
googleStreets.addTo(map);

var googleSat = L.tileLayer(
    'http://{s}.google.com/vt?lyrs=s&x={x}&y={y}&z={z}',
    {
        maxZoom: 20,
        subdomains:['mt0','mt1','mt2','mt3']
    }
);
googleSat.addTo(map);

var googleTerrain = L.tileLayer(
    'http://{s}.google.com/vt?lyrs=p&x={x}&y={y}&z={z}',
    {
        maxZoom: 20,
        subdomains:['mt0','mt1','mt2','mt3']
    }
);
googleTerrain.addTo(map)


// Layer controller
var baseMaps = {
    "OSM": osm,
    // "Topology": OpenTopoMap,
    "Dark": CartoDB_DarkMatter,
    "Voyager": CartoDB_Voyager,
    "Streets": googleStreets,
    // "Hybrid": googleHybrid,
    "Satellite": googleSat,
    "Terrain": googleTerrain
}

var overlayMaps = {
    // "Start Marker": startMarker
};

L.control.layers(baseMaps, overlayMaps, {'collapsed': false}).addTo(map);

new L.cascadeButtons([
    {icon: 'fas fa-location-dot', ignoreActiveState:true , command: () =>{
        markerNodeList.endMarkerNode.marker.openPopup();
        map.panTo(markerNodeList.endMarkerNode.getLatLng());
     }},
    {icon: 'fas fa-location-dot', ignoreActiveState:true , command: () =>{
        markerNodeList.startMarkerNode.marker.openPopup();
        map.panTo(markerNodeList.startMarkerNode.getLatLng())
     }},
], {position:'topright', direction:'horizontal'}).addTo(map);

L.control.zoom({
    position: 'topright'
}).addTo(map);

var geocoder = L.Control.geocoder({
    defaultMarkGeocode: false
  })
    .on('markgeocode', function(e) {
      var bbox = e.geocode.bbox;
      var poly = L.polygon([
        bbox.getSouthEast(),
        bbox.getNorthEast(),
        bbox.getNorthWest(),
        bbox.getSouthWest()
      ]).addTo(map);
      map.fitBounds(poly.getBounds());
    });
geocoder.addTo(map);

var searchButton = L.easyButton('fas fa-route',
    {
        position: 'topright',
        states: [{
            stateName: 'search-route',
            icon: 'fas fa-route',
            onClick: function (button, map) {
                button.state('loading')
                markerNodeList.updateAllPaths()
                ckptView.updateCkptRows()

                if (markerNodeList.lockCount === 0){
                    button.state('search-route');
                }
            }
        },
        {
            stateName: 'loading',
            icon: 'fas fa-spinner',
            onClick: function (button, map){
                // button.state('search-route')
        }
        }]
}).addTo( map );

var sidebar = L.control.sidebar({ container: 'sidebar' })
            .addTo(map);

sidebar.open('checkpoints')

sidebar.on('content', function (ev) {
    switch (ev.id) {
        case 'autopan':
        sidebar.options.autopan = true;
        break;
        default:
        sidebar.options.autopan = false;
    }
});


map.on('mousemove', function(e) {
    document.getElementsByClassName('coordinate')[0].innerHTML = 'lat: ' + e.latlng.lat.toFixed(decimalPlaces)
        + '\tlon: ' + e.latlng.lng.toFixed(decimalPlaces);
})

map.on('click', function(e) {
    sidebar.close()
});


getTwoNearNodes({
    graph: getGraphType()[0],
    edgeBased: getGraphType()[1],
    maxDistance: 33_300_00
})
.then(response=>response.json())
.then(data => {

    markerNodeList = new MarkerNodeList(
        map,
        [data['lat1'], data['lon1']],
        [data['lat2'], data['lon2']],
        numPathTypes
    )
    ckptView = new ckptTable(markerNodeList);


    markerNodeList.updateAllPaths();
    ckptView.updateCkptRows(markerNodeList)
    map.panTo(new L.LatLng((data['lat1'] + data['lat2']) / 2, (data['lon1'] + data['lon2']) / 2))


    // setCkptTableContents(0);
    // setCkptTableContents(-1);
    // calculateShortestPath(0, getGraphType())
})


/* ===========================================================
            Others (unclassified)
============================================================== */
