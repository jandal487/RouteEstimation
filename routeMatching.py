# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:33:22 2016

@author: Amir Abu Jandal
"""
import math
import folium
from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from math import radians, cos, sin, asin, sqrt
import timeit
from IPython.display import HTML, display

#### mainMethod(): The main method to execute the whole algorithm
def mainMethod():    
    #==>[1]: Read the trip points data into a dataframe    
    tripPoints_file = "/testAlgo.csv"
    sample_df = pd.read_csv(tripPoints_file, sep=";")
    #header_file = "/trippoint_headerInfo.txt"    
    #sample_df.columns = open(header_file, 'r').readlines()[0].split(',')
    
    #==>[2]: Get a connection to PostgreSQL
    dbParameters = "dbname='trieste' user='postgres' host='localhost' password='123456'"
    con = getDBconnection(dbParameters)
    
    #==>[3]: Initialize a map and the final required dataframe
    map_osm = folium.Map(location=[45.6, 13.7], zoom_start= 11)
    map_osm.save(outfile='/osm.html')
    
    tripMatchDF = pd.DataFrame()
    
    #==>[4]: Loop through every trip and process it for matching and routing
    unique_ids = sample_df["trip_id"].unique()
    for i in range(len(unique_ids)):
        trip_id = unique_ids[i]
        dfSubset = sample_df[sample_df['trip_id'] == trip_id]
        
        #Getting appropriate index values for the dfSubset dataframe
        dfSubset.index = range(len(dfSubset.index.values))
        
        # Sorting the trip from start to the end
        dfSubset = dfSubset.sort(["corrected_timestamp"], ascending=[True])
        dfSubset = dfSubset.sort(["point_timestamp"], ascending=[True])
        
        # Getting the road segments and matching points for the provided GPS trippoints
        start_time = timeit.default_timer()
        returnedList = processTrip(dfSubset, con, map_osm)
        elapsed = timeit.default_timer() - start_time
        print("Time taken to process 1 trip is " + str(elapsed) + " seconds.")
        
        tripMatchDF.append(returnedList[0])
        map_osm = returnedList[1]
        map_osm.save(outfile='/osm_result.html')
    #END of FOR Loop
#END of mainMethod()  

#### getDBconnection(): get connection to spatial database
def getDBconnection(dbParameters):
    try:
        engine = create_engine('postgresql://postgres:123456@localhost:5432/trieste', encoding='latin1')
        con = engine.connect()
    except:
        print("Couldn't connect to the database!")
    
    return(con)
#END of getDBconnection()
    
#### processTrip(): process each trip to find matching points and road segments
def processTrip(dfSubset, con, map_osm):
    ways_df = pd.DataFrame(np.nan, index=range(len(dfSubset)), columns=list(["trip_id", "osm_id", "name", "mylinedistance", "way"]))
    ways_df.iloc[0] = [0, 0, "", 0, ""]
    
    matchedNode_df = pd.DataFrame(np.nan, index=range(len(dfSubset)), columns=list(["trip_id", "matchedLat", "matchedLon"]))
    matchedNode_df.iloc[0] = [0, 0, 0]
    
    trips_df_columns = dfSubset.columns.union(ways_df.columns.union(matchedNode_df.columns))
    trips_df = pd.DataFrame(np.nan, index=range(len(dfSubset)), columns=trips_df_columns)

    # For Snapping the road segments
    roadSnap_df = pd.DataFrame([list([""])], index=range(len(dfSubset)-1), columns=["roadSnap"])

    #ID of the trip being processed
    trip_id = dfSubset['trip_id'][0]
    
    #Road name of the previous trippoint
    previousRoadName = u''
    
    for j in range(len(dfSubset)):
        tripPoint = dfSubset.iloc[[j]]
        gpsLat = tripPoint['latitude'][j]
        gpsLon = tripPoint['longitude'][j]
            
        ways_df = nearestNeighbourWay(tripPoint, con, trip_id)
        
        #For the starting trippoint j==0
        if j == 0:
            selectedWay_df = ways_df.iloc[[0]]
            
            road = selectedWay_df['way'][0][11:len(selectedWay_df['way'][0])-1].split(',')
            matchingPointData = getMatchingPointData(road, gpsLat, gpsLon)
            matchedCoordinates = matchingPointData[0]
            matchedNode_df.iloc[0] = [trip_id, matchedCoordinates[0], matchedCoordinates[1]]
            
            # After processing trippoint, setting it as previous road name for processing the upcoming trippoint                
            previousRoadName = ways_df['name'][0]
        elif len(ways_df) < 1:
            matchedNode_df.iloc[0] = [trip_id, 0, 0]
            
            ways_df = pd.DataFrame(np.nan, index=[0], columns=list(["trip_id", "osm_id", "name", "mylinedistance", "way"]))
            ways_df.iloc[0] = [trip_id, 0, "", 0, ""]
            
        elif len(ways_df) > 0:
            shortestPathList = [1000]*len(ways_df)
            
            if sum(ways_df['name'].isin([previousRoadName])) > 0:
                bestIndex = list(ways_df['name']).index(previousRoadName)

                roadString = ways_df['way'][bestIndex][11:len(ways_df['way'][bestIndex])-1]
                road = roadString.split(',')
                matchingPointData = getMatchingPointData(road, gpsLat, gpsLon)
                matchedCoordinates = matchingPointData[0]
                matchedNode_df.loc[j] = [trip_id, matchedCoordinates[0], matchedCoordinates[1]]                
                
                shortestPathList[0] = getShortestPathResult(con, j, matchingPointData, matchedNode_df)
                
                roadSnap_df.append(shortestPathList[0])
                roadSnap_df.loc[j-1] = roadString
                
                # After processing trippoint, setting it as previous road name for processing the upcoming trippoint                
                previousRoadName = ways_df['name'][bestIndex]
            else:
                #SwitchArray is used to classify btw tricky and non-tricky points        
                switchArray = getSwitchArray(ways_df)
                
                #Variables to decide actual road for tricky trippoints
                distanceCost = [1000]*len(ways_df)
                pathCost = [1000]*len(ways_df)
                nearestWay = [1000]*len(ways_df)
                wayEstimationValue = [1000]*len(ways_df)
        
                for q in range(sum(switchArray)):
                    selectedWay_df = ways_df.iloc[[q]]
                    
                    road = selectedWay_df['way'][q][11:len(selectedWay_df['way'][q])-1].split(',')
                    matchingPointData = getMatchingPointData(road, gpsLat, gpsLon)

                    shortestPathList[q] = getShortestPathResult(con, j, matchingPointData, matchedNode_df)
                    
                    matchingPointDistance = matchingPointData[1]
                    shortestPathCost = shortestPathList[q]['cost'].sum()
                    nearestWayDistance = ways_df.iloc[q]['mylinedistance']
                    
                    # FINALLY, the Decision Value can be stated as w = x + y + z
                    distanceCost[q] = matchingPointDistance
                    pathCost[q] = shortestPathCost
                    nearestWay[q] = nearestWayDistance
                    wayEstimationValue[q] = matchingPointDistance + shortestPathCost + nearestWayDistance
                #END of FOR Loop
                
                # Getting the index of the way with least cost
                bestIndex = wayEstimationValue.index(min(wayEstimationValue))
                
                # Extracting the road segment
                roadString = ways_df['way'][bestIndex][11:len(ways_df['way'][bestIndex])-1]
                road = roadString.split(',')
                
                # Getting the matched points on the road segment
                matchingPointData = getMatchingPointData(road, gpsLat, gpsLon)
                matchedCoordinates = matchingPointData[0]
                matchedNode_df.loc[j] = [trip_id, matchedCoordinates[0], matchedCoordinates[1]]
                
                # Snapping the road segments 
                roadSnap_df.loc[j-1] = shortestPathList[bestIndex]['st_astext']

                # After processing trippoint, setting it as previous road name for processing the upcoming trippoint                
                previousRoadName = ways_df['name'][bestIndex]
            #END of inner ELSE
        #END of last ELSE
        
        trips_df.loc[[j]] = pd.merge(pd.merge(tripPoint, selectedWay_df, how='outer', on=['trip_id']), matchedNode_df, how='outer', on=['trip_id'])
    #END of for loop        
    
    map_osm = viewMatchedPointsOnMap(map_osm, trips_df)
    
    map_osm = viewRoadSnapOnMap(map_osm, roadSnap_df, trips_df)
    
    return([trips_df, map_osm])
#END of processTrip()

#### viewMatchedPointsOnLeaflet(): View gps points and matched points on leaflet
def viewMatchedPointsOnMap(m, trips_df):
    gpxData = [trips_df['latitude'], trips_df['longitude']]
    matchedData = [trips_df['matchedLat'], trips_df['matchedLon']]
    
    for k in range(len(trips_df)):
        folium.CircleMarker(location=[gpxData[0][k], gpxData[1][k]], radius=4, color='#31a354', fill_color='#31a354').add_to(m)
        folium.CircleMarker(location=[matchedData[0][k], matchedData[1][k]], radius=5, color='red', fill_color='red').add_to(m)
    #END of for loop
    
    return(m)
#END of viewMatchedPointsOnLeaflet()
    
#### viewRoadSnapOnMap(): View road snap for a trip on leaflet
def viewRoadSnapOnMap(map_osm, roadSnap_df, trips_df):
    if(len(roadSnap_df) > 0):
        for w in range(len(roadSnap_df)):
            roadSegment = roadSnap_df['roadSnap'][w][11:len(roadSnap_df['roadSnap'][w])-1].split(',')

            lonVector = [0]*len(roadSegment)
            latVector = [0]*len(roadSegment)            
            
            for u in range(len(roadSegment)):
                latlon = roadSegment[u].split(" ")
                lonVector[u] = latlon[0]
                latVector[u] = latlon[1]
            #END of FOR loop
            
            for k in range(len(latVector)):
                folium.CircleMarker(location=[latVector[k], lonVector[k]], radius=2.5, color='#blue', fill_color='#blue').add_to(map_osm)
                #folium.PolyLine([latVector[k], lonVector[k]], color="black", weight=2, opacity=1).add_to(map_osm)
            #END of for loop            
        # END of outer FOR LOOP
            
    return(map_osm)
#END of viewRoadSnapOnMap()

#### getShortestPathResult(): Execute Dijkstra to find shortest path between 2 osm nodes
def getShortestPathResult(con, j, matchingPointData, matchedNode_df):
    matchedCoordinates = matchingPointData[0]
    
    sourceLon = round(matchedNode_df.iloc[j-1][2], 5)
    sourceLat = round(matchedNode_df.iloc[j-1][1], 5)
    targetLon = round(matchedCoordinates[1], 5)
    targetLat = round(matchedCoordinates[0], 5)
            
    sourceNodeQuery = ("SELECT source FROM hh_2po_4pgr ORDER BY geom_way <-> "+
                        "ST_SetSRID(ST_Makepoint(%f,%f), 4326) LIMIT 1;") % (sourceLon, sourceLat)
    targetNodeQuery = ("SELECT target FROM hh_2po_4pgr ORDER BY geom_way <-> "+
                        "ST_SetSRID(ST_Makepoint(%f,%f), 4326) LIMIT 1;") % (targetLon, targetLat)
                                    
    sourceNode = pd.read_sql_query(sourceNodeQuery,con=con).iloc[0][0]
    targetNode = pd.read_sql_query(targetNodeQuery,con=con).iloc[0][0]
    
    #boundingBox = getBoundingBox(targetLat, targetLon, 5.0)
    
    #mnLon = boundingBox[0]
    #mnLat = boundingBox[1] 
    #mxLon = boundingBox[2]
    #mxLat = boundingBox[3]
    
    if sourceNode != targetNode :            
        shortestPathQuery = ("SELECT seq, id1 AS node, id2 AS edge, di.cost, ST_AsText(ST_Transform(geom_way, 4326))"+ 
                            " FROM pgr_dijkstra('SELECT id, source, target, st_length(geom_way) as cost FROM hh_2po_4pgr WHERE "+
                            " public.hh_2po_4pgr.geom_way && ST_MakeEnvelope(13.724327,45.606592,13.848267,45.70534, 4326)',"
                            "%d, %d, false, false ) as di JOIN hh_2po_4pgr ON di.id2 = hh_2po_4pgr.id") % (sourceNode,targetNode)
        #print(shortestPathQuery)
        return(pd.read_sql_query(shortestPathQuery,con=con).iloc[0])
    else:
        return(pd.DataFrame([-1000], index=[0], columns=["cost"]))
#END of getShortestPathResult()

#### nearestNeighbourWay(): Get nearest way from PostGIS using the getBoundingBox() function
def nearestNeighbourWay(tripPoint, con, trip_id):    
    radLat = tripPoint['latitude']
    radLon = tripPoint['longitude']
    boundingBox = getBoundingBox(radLat, radLon, 1.2)
    
    minLon = boundingBox[0]
    minLat = boundingBox[1] 
    maxLon = boundingBox[2]
    maxLat = boundingBox[3]
    
    spatialQuery = ("with distance_query as(SELECT DISTINCT ON (way)"
    +" osm_id, name, ST_Distance(ST_Transform(ST_GeomFromText('POINT(%f %f)',4326),3857),"
    +" ST_Transform(way, 3857)) AS myLineDistance, ST_AsText(ST_Transform(way, 4326)) as way"
    +" FROM public.planet_osm_line"
    +" WHERE public.planet_osm_line.way && ST_Transform("
    +" ST_MakeEnvelope(%f,%f,%f,%f, 4326), 3857)"
    +" AND ST_Distance(ST_Transform(ST_GeomFromText('POINT(%f %f)',4326),3857),"
    +" ST_Transform(way, 3857)) < 10.0 AND osm_id > 0)"
    +" SELECT *FROM distance_query" 
    +" ORDER BY myLineDistance LIMIT 3;"
    ) % (radLon, radLat, minLon, minLat, maxLon, maxLat, radLon, radLat) 

    #way: Execute the query and get the results in a simple df with geometry as text
    ways_table = pd.read_sql_query(spatialQuery,con=con)
    
    #Keep trip_id in the ways dataframe, to keep track and later on concat results
    ways_table['trip_id'] = [trip_id]*len(ways_table)
  
    return(ways_table)  
#END of nearestNeighbourWay()

#### getBoundingBox(): Get a bounding box around lat/lon for some distance
def getBoundingBox(radLat, radLon, distance):
    eRadius = 6371.01
    
    pRadius = eRadius*math.cos(radLat)
    
    minLat = radLat - distance/eRadius
    maxLat = radLat + distance/eRadius
    minLon = radLon - distance/pRadius
    maxLon = radLon + distance/pRadius
    
    boundingBox = list([minLon, minLat, maxLon, maxLat])
    
    return(boundingBox)
#END of getBoundingBox()

#### Calculate the great circle distance between two points on the earth
def haversine(lon1, lat1, lon2, lat2):
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r
#END of haversine
  
#### getMatchingPointData(): Get the matching point to the GPS point on the road
def getMatchingPointData(road, gpsLat, gpsLon):
    distances = []
    
    for node in road:
        latlon = node.split(" ")
        lat = float(latlon[1])
        lon = float(latlon[0])
        
        d = haversine(lon, lat, gpsLon, gpsLat)
        distances.append(d)
    #END of for loop
        
    minimumNodeIndex = distances.index(min(distances))
    matchedNode = road[minimumNodeIndex]
    latlon = matchedNode.split(" ")
    matchedLat = float(latlon[1])
    matchedLon = float(latlon[0])
  
    returnList = []
    returnList = [[matchedLat, matchedLon], distances[minimumNodeIndex]]
  
    return(returnList)
#END of getMatchingPointData()
        
#### getSwitchArray(): Array used later to classify tricky and non-tricky trippoints
def getSwitchArray(ways_df):
    wayDistances = ways_df['mylinedistance']
    meanValue = wayDistances.mean()
    returnArray = [0]*len(ways_df)
    
    if len(wayDistances) > 1:
        returnArray[0] = 1
        for x in range(len(wayDistances)):
            #Check to take care of last x value
            if x < len(wayDistances)-1:
                diff = wayDistances[x+1] - wayDistances[0]
                #Checks to classify tricky and non-tricky trippoints    
                if diff == 0:
                    break
                elif diff < meanValue:
                    returnArray[x+1] = 1
                else:
                    break
        #END of FOR loop
    else:
        returnArray[0] = 1
    #END of OUTER ELSE
    
    return(returnArray)
#END of getSwitchArray()
