## R-Solution to process trippoints, route trippoints and plot 
## the route with road segments on map
## Author: Amir Abu Jandal

# ==> getSampleDF(): Function to get sample data in a dataframe
getSampleDF = function(file1, file2){
  
  #Reading the trippoints sample data
  sample_df <-read.csv(file1, header=FALSE, stringsAsFactors=TRUE,sep="\t",quote="")
  
  #Reading column names of the trippoints data
  names(sample_df) = read.table(file2, header = FALSE)[[1]]
  names(sample_df) = gsub(",", "", names(sample_df))
  
  return(sample_df)
}# END of getSampleDF()

# ==> getDatabaseConnection(): R connectivity with PostgreSQL
getDatabaseConnection = function(dbName, hostName, portNumber, userName, password){
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname=dbName, host=hostName, port=portNumber, user=userName, password=password)
  
  return(con)
}# END of getDatabaseConnection()

# ==> nearestNeighbourWay(): Get nearest way from PostGIS using the getBoundingBox() function
nearestNeighbourWay = function(tripPoint, dbConnection){
  lat = tripPoint$latitude
  lon = tripPoint$longitude
  boundingBox <- getBoundingBox(lat, lon, 1.2)
  
  minLon = boundingBox[1]
  minLat = boundingBox[2] 
  maxLon = boundingBox[3]
  maxLat = boundingBox[4]

  spatialQuery = paste("with distance_query as(SELECT DISTINCT ON (way)
                       osm_id, name, ST_Distance(ST_Transform(ST_GeomFromText('POINT(",lon," ",lat,")',4326),3857), ST_Transform(way, 3857)) AS myLineDistance, ST_AsText(ST_Transform(way, 4326)) as way
                       FROM public.planet_osm_line 
                       WHERE public.planet_osm_line.way && ST_Transform(ST_MakeEnvelope(",minLon,",", minLat,",", maxLon,",", maxLat,", 4326), 3857)
                       AND ST_Distance(ST_Transform(ST_GeomFromText('POINT(",lon," ",lat,")',4326),3857), ST_Transform(way, 3857)) < 10.0
                       AND osm_id > 0)
                       SELECT *FROM distance_query 
                       ORDER BY myLineDistance LIMIT 3;")
  
  #way: Execute the query and get the results in a simple df with geometry as text $$$ PLEASE UPDATE THE ABOVE QUERY $$$
  way = dbGetQuery(conn = con, spatialQuery)
  
  return(way)  
}# END of nearestNeighbourWay()

# ==> getMatchingPointData(): Get the matching point to the GPS point on the road
getMatchingPointData = function(road, gpsLat, gpsLon){
  distances = vector()
  
  for(node in road){
    latlon = strsplit(node, " ")[[1]]
    lat = as.numeric(latlon[2])
    lon = as.numeric(latlon[1])
    
    #d = distm(c(lat, lon), c(gpsLat, gpsLon), fun = distHaversine)[1][1]
    d = distGeo(matrix(c(lat, lon), ncol = 2), matrix(c(gpsLat, gpsLon), ncol = 2))/1609.34
    distances = c(distances, d)
  }
  
  minimumNodeIndex = which(distances == min(distances))
  matchedNode = road[minimumNodeIndex]
  latlon = strsplit(matchedNode, " ")[[1]]
  matchedLat = as.numeric(latlon[2])
  matchedLon = as.numeric(latlon[1])
  
  matchedCoordinates = c(matchedLat, matchedLon)
  
  return(list(matchedCoordinates=matchedCoordinates, minimumDistance=distances[minimumNodeIndex]))
}#END of getMatchingPointData()

# ==> viewMatchedPointsOnLeaflet(): View gps points and matched points on leaflet
viewMatchedPointsOnLeaflet = function(m, tripDF){
  
  gpxData = cbind(tripDF$latitude, tripDF$longitude)
  matchedData = cbind(tripDF$matchedLat, tripDF$matchedLon)
  
  for(k in 1:nrow(gpxData)){
    if(k == 1){
      m <- addCircleMarkers(m, color = "#31a354", radius = 4, lng = gpxData[k, 2], lat = gpxData[k, 1])   
      m <- addCircleMarkers(m, color = "red", radius = 5, lng = matchedData[k, 2], lat = matchedData[k, 1])
      
    }else if( k == nrow(gpxData)){
      m <- addCircleMarkers(m, color = "#31a354", radius = 4, lng = gpxData[k, 2], lat = gpxData[k, 1]) 
      m <- addCircleMarkers(m, color = "red", radius = 5, lng = matchedData[k, 2], lat = matchedData[k, 1])    
    }else{
      m <- addCircleMarkers(m, color = "#31a354", radius = 4, lng = gpxData[k, 2], lat = gpxData[k, 1]) 
      m <- addCircleMarkers(m, color = "red", radius = 5, lng = matchedData[k, 2], lat = matchedData[k, 1])
    }
  }
  
  return(m)
}#END of viewMatchedPointsOnLeaflet()

### SEGEMENTS can be EXTRACTED from THIS function!!!
# ==> viewRoadSnapOnLeaflet(): View road snap for a trip on leaflet
viewRoadSnapOnLeaflet = function(m, roadSnapDF, tripMatchDF){
  
  if(nrow(roadSnapDF) > 0){
    roadSegment = ""
    lonSegmentVector = vector()         # CAN BE USED TO COLLECT THE SEGMENTS
    latSegmentVector = vector()         # CAN BE USED TO COLLECT THE SEGMENTS
    startLon = vector()
    startLat = vector()
    endLon = vector()
    endLat = vector()
    
    for(w in 1:nrow(roadSnapDF)){
      lonVector = vector()
      latVector = vector()
      
      roadSegmentList = substr(roadSnapDF[w, ]$st_astext, 12, nchar(roadSnapDF[w, ]$st_astext)-1)
      roadSegment = strsplit(roadSegmentList, ",")[[1]]
      
      for(u in 1:length(roadSegment)){
        latlon = strsplit(roadSegment[u], " ")[[1]]
        lonVector = c(lonVector, as.numeric(latlon[1]))
        latVector = c(latVector, as.numeric(latlon[2]))
      }
      
      m <- addCircleMarkers(m, color = "blue", radius = 2, lng = lonVector, lat = latVector)
      m <- addPolylines(m, color = "black", lng = lonVector, lat = latVector)
      '
      if(roadSnapDF[w, ]$seq == 0){
        endLon = lonVector[1]
        endLat = latVector[1]
        
        if(length(startLon)>0 && length(startLon)>0){
          m <- addPolylines(m, color = "black", lng = lonSegmentVector, lat = latSegmentVector)
        }
        startLon = vector()
        startLat = vector()
      }else{
        startLon = lonVector[u]
        startLat = latVector[u]
      }
      
      lonSegmentVector = c(startLon, endLon)
      latSegmentVector = c(startLat, endLat)'
    }# END of outer FOR LOOP
  }
  
  return(m)
}#END of viewRoadSnapOnLeaflet()

# ==> getTripMatchDF()
getTripMatchDF = function(dfSubset, con){
  n = list("osm_id", "names", "mylinedistance", "way")
  way <- data.frame(matrix(nrow = 0, ncol = length(n)))
  names(way) <- n
  class(way$osm_id) <- "numeric"
  class(way$names) <- "character"
  class(way$mylinedistance) <- "numeric"
  class(way$way) <- "character"
  
  n = list("matchedLat", "matchedLon")
  matchedNode <- data.frame(matrix(nrow = 0, ncol = length(n)))
  names(matchedNode) <- n
  class(matchedNode$matchedLat) <- "numeric"
  class(matchedNode$matchedLon) <- "numeric"
  
  tripDF = merge(dfSubset[0, ], matchedNode[0, ])
  tripDF = merge(tripDF[0, ], way[0, ])
  class(tripDF$osm_id) <- "numeric"
  class(tripDF$names) <- "character"
  class(tripDF$mylinedistance) <- "numeric"
  class(tripDF$way) <- "character"
  
  # For Snapping the road segments
  roadSnapDF <- data.frame()
  
  # Tracking previous way's Name
  previousWayId = ""
  
  for(i in 1:nrow(dfSubset)){
    tripPoint <- dfSubset[i, ]
    way <- nearestNeighbourWay(tripPoint, con)
    
    if(i == 1){
      way <- way[1, ]
      roadString <- substr(way$way, 12, nchar(way$way)-1)
      road <- strsplit(roadString, ',')[[1]]
      matchedCoordinates <- getMatchingPointData(road, as.numeric(tripPoint$latitude), as.numeric(tripPoint$longitude))[[1]]
      matchedNode[i, ] <- c(matchedCoordinates[1], matchedCoordinates[2])
      
      tripDF[i, ] <- merge(merge(tripPoint, matchedNode), way)
    }else if(nrow(way) > 0){
      gpsLat = as.numeric(tripPoint$latitude)
      gpsLon = as.numeric(tripPoint$longitude)
      
      # w is the estimation value for TOP 3 ways. Where as, (w = ä.d + ß.t)
      distanceCost = vector(length = nrow(way))
      pathCost = vector(length = nrow(way))
      nearestWay = vector(length = nrow(way))
      wayEstimationValue = vector(length = nrow(way))
      shortestPathList = list(length = nrow(way))
      
      alpha = 0.35
      beta = 0.45
      gama = 1 - (alpha+beta)
      
      # Initial check to avoid caluclation of further ways
      #switchArray <- getSwitchArray(way)
      
      for(q in 1:nrow(way)){
        #if(q > 1 && switchArray[q-1] == 0){
        #  for(r in q:nrow(way)){
        #    wayEstimationValue[r] = 10000
        #  }
        #  break
        #}
        
        roadString <- substr(way[q, ]$way, 12, nchar(way[q, ]$way)-1)
        road <- strsplit(roadString, ',')[[1]]
        
        # Get matching point for GPS on the nearest road or way
        matchingPointData <- getMatchingPointData(road, gpsLat, gpsLon)
        matchedCoordinates <- matchingPointData[[1]]
        
        sourceLon <- format(round(matchedNode[i-1, ][2], 5), nsmall = 5)
        sourceLat <- format(round(matchedNode[i-1, ][1], 5), nsmall = 5)
        targetLon <- format(round(matchedCoordinates[2], 5), nsmall = 5)
        targetLat <- format(round(matchedCoordinates[1], 5), nsmall = 5)
        
        sourceNodeQuery <- paste("SELECT source FROM hh_2po_4pgr 
                                ORDER BY geom_way <-> ST_SetSRID(ST_Makepoint(", sourceLon, ", ", sourceLat,"), 4326) LIMIT 1;")
        targetNodeQuery <- paste("SELECT target FROM hh_2po_4pgr 
                                ORDER BY geom_way <-> ST_SetSRID(ST_Makepoint(", targetLon, ", ", targetLat,"), 4326) LIMIT 1;")
          
        sourceNode <- dbGetQuery(conn = con, sourceNodeQuery)[1, 1]
        targetNode <- dbGetQuery(conn = con, targetNodeQuery)[1, 1]
        
        shortestPathQuery <- paste("SELECT seq, id1 AS node, id2 AS edge, di.cost, ST_AsText(ST_Transform(geom_way, 4326)) 
                                   FROM pgr_dijkstra('SELECT id, source, target, st_length(geom_way) as cost FROM hh_2po_4pgr',",
                                   sourceNode, ", ", targetNode,", false, false ) as di JOIN hh_2po_4pgr ON di.id2 = hh_2po_4pgr.id")
        
        shortestPathList[[q]] <- dbGetQuery(conn = con, shortestPathQuery)

        matchingPointDistance <- matchingPointData[[2]]
        shortestPathCost <- sum(shortestPathList[[q]]$cost)
        nearestWayDistance <- way[q, ]$mylinedistance
        
        # FINALLY, the decision value can be stated as w = (ä.x) + (ß.y) + (g.z)
        distanceCost[q] <- matchingPointDistance
        pathCost[q] <- shortestPathCost
        nearestWay[q] <- nearestWayDistance
        
        if(is.na(way$name[q])){
          way$name[q] <- ""
        }
        
        if(previousWayId == way$name[q]){
          wayEstimationValue[q] <- -1
          break
        }else{
          wayEstimationValue[q] <- (alpha * matchingPointDistance) + (beta * shortestPathCost) + (gama * nearestWayDistance)   
        }
      }# END of q FOR LOOP
        
      # Getting the index of the way with least cost
      bestIndex <- which.min(wayEstimationValue)
      print(paste("Selected way = ", bestIndex))
      print(paste("Way Estimation values = ", wayEstimationValue))
      
      # Extracting the road segmetn
      roadString <- substr(way[bestIndex, ]$way, 12, nchar(way[bestIndex, ]$way)-1)
      road <- strsplit(roadString, ',')[[1]]
      
      # Getting the matched points on the road segment
      matchedCoordinates <- getMatchingPointData(road, gpsLat, gpsLon)[[1]]
      matchedNode[i, ] <- c(matchedCoordinates[1], matchedCoordinates[2])
      
      # Snapping the road segments 
      roadSnapDF <- rbind(roadSnapDF, shortestPathList[[bestIndex]])
      
      tripDF[i, ] <- merge(merge(tripPoint, matchedNode[i, ]), way[bestIndex, ])
    }else{
      matchedNode[1, ] <- c(0, 0)
      
      n <- list("osm_id", "names", "mylinedistance", "way")
      way <- data.frame(matrix(nrow = 0, ncol = length(n)))
      names(way) <- n
      class(way$osm_id) <- "numeric"
      class(way$names) <- "character"
      class(way$mylinedistance) <- "numeric"
      class(way$way) <- "character"
      way[1, ] <- c(0, "", 0, "")
      
      tripDF[i, ] <- merge(merge(tripPoint, matchedNode), way)
    }
    
    previousWayId = tripDF[i, ]$names
    if(is.na(previousWayId)){
      previousWayId <- ""
    }
    print(paste(100*i/nrow(dfSubset), "% has been completed!"))
  }# END of OUTER FOR LOOP
  
  # Plot the matched points and the GPS points
  m <- viewMatchedPointsOnLeaflet(m, tripDF)
  
  # Plotting the road segment
  m <- viewRoadSnapOnLeaflet(m, roadSnapDF, tripDF)
  
  return(list(trip=tripDF, map=m))
}# END of getTripMatchDF()

# ==> getBoundingBox(): Get a bounding box around lat/lon for some distance
getBoundingBox = function(radLat, radlon, distance){
  
  eRadius <- 6371.01
  
  pRadius = eRadius*cos(radLat)
  
  minLat <- radLat - distance/eRadius
  maxLat <- radLat + distance/eRadius
  minLon <- radlon - distance/pRadius
  maxLon <- radlon + distance/pRadius
  
  boundingBox <- c(minLon, minLat, maxLon, maxLat)
  
  return(boundingBox)
}# END of getBoundingBox()

# ==> getSwitchArray(): Initial check to avoid calculation of further ways and to mostly select first way
getSwitchArray = function(way){
  wayDistances <- way$mylinedistance
  returnArray <- vector(length = length(wayDistances)-1)
  
  if(length(wayDistances) > 1){
    for(j in 1:length(wayDistances)-1){
      diff <- abs(wayDistances[1] - wayDistances[j+1])
      
      print(paste("DIFFERENCE = ", diff))
      
      if(diff < 10){
        returnArray[j] = 1
      }else{
        returnArray[j] = 0
        break
      }
    }    
  }else{
    returnArray[1] = 0
  }
  print(paste("returnArray = ", returnArray))
  return(returnArray)
}# END of getSwitchArray()

# ==> mainFunction(): MAIN function
mainFunction = function(){
  library(sp)
  library(geosphere)
  library(leaflet)
  library(RPostgreSQL)
  
  # Reading the GPS data
  #tripPoints_file <- "/trippoint_sample2.tsv"
  tripPoints_file <- "/testAlgo.csv"
  header_file <- "/trippoint_headerInfo.txt"
  #sample_df <- getSampleDF(tripPoints_file, header_file)
  sample_df <- read.csv(tripPoints_file, header=FALSE, stringsAsFactors=FALSE,sep=",")
  names(sample_df) <- read.table(header_file, header = FALSE)[[1]]
  names(sample_df) <- gsub(",", "", names(sample_df))
  
  # Opening a connection to the database
  con <- getDatabaseConnection("trieste", "localhost", 5432, "postgres", "123456")
  
  # Initializing the Leaflet map variable
  options(viewer = NULL)
  m <- leaflet()
  m <- setView(m, lng = 13.7, lat = 45.6, zoom = 11)
  m <- addTiles(m)
  
  # Final dataframe with trippoints and their matching coordinates and roads
  tripMatchDF <- data.frame()
  
  # Unique trips
  trip_id <- unique(as.character(sample_df[,c("trip=trip_id")]))
  
  for(k in 1:length(trip_id)){
    startTime <- Sys.time()
    
    id <- trip_id[k]
    dfSubset <- sample_df[sample_df$'trip=trip_id' == id, ]
    
    # Ordering the trip from start to the end
    dfSubset <- dfSubset[order(dfSubset$corrected_timestamp), ]
    dfSubset <- dfSubset[order(dfSubset$point_timestamp), ]
    
    # Getting the road segments and matching points for the provided GPS trippoints
    returnedList <- getTripMatchDF(dfSubset, con)
    tripMatchDF <- rbind(tripMatchDF, returnedList$trip)
    m <- returnedList$map
    
    endTime <- Sys.time()
    executionTime <- endTime - startTime
    executionTime
  }# END of FOR LOOP
  
  View(tripMatchDF)
  
  ### UNCLEAR TRIPPOINTS
  # Reading the unclear trippoints
  tripPoints_file = "/unclearTrippoints.csv"
  unclear_df <-read.csv(tripPoints_file, header=TRUE, stringsAsFactors=FALSE, sep=",")
  
  # Initializing the Leaflet map variable
  options(viewer = NULL)
  m <- leaflet()
  m <- setView(m, lng = 13.7, lat = 45.6, zoom = 11)
  m <- addTiles(m)
  
  # Final dataframe with trippoints and their matching coordinates and roads
  tripMatchDF <- data.frame()
  
  for(i in 1:nrow(unclear_df)){
    # Getting the roads for provided GPS trippoints
    tripMatchDF <- rbind(tripMatchDF, getTripMatchDF(unclear_df[i, ], con))
  }
  
  
  #Disconnect from the server and the database
  dbDisconnect(con)
  dbUnloadDriver(drv)
  
}# END of mainFunction()