import numpy as np
import osmnx as ox
import matplotlib.pyplot as plt
# from shapely.geometry import MultiPoint
from scipy.spatial import Delaunay
from shapely.geometry import Polygon

#%% -----------------------------------------------------------------

def dPointLine(x0,x1,x2):
    # imputs have to be x,y,z numpy arrays
    ## https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    x0 = np.append(x0, 0)
    x1 = np.append(x1, 0)
    x2 = np.append(x2, 0)
    d = np.linalg.norm(np.cross(x0-x1,x0-x2))/np.linalg.norm(x2-x1)
    return d

#%% -----------------------------------------------------------------

def extractUTMscenery(centerOfScene, tags, radiusOfScene):

    gdf = ox.geometries_from_point(centerOfScene, tags, dist=radiusOfScene) # using lon,lat coords
    gdf_proj = ox.project_gdf(gdf)  # using UTM coords
    bbox = ox.utils_geo.bbox_from_point(point=centerOfScene, dist=radiusOfScene, project_utm=True)
    return gdf_proj, bbox

#%% ----------------------------------------------

def extractUTMscenery2(BBOX, tags):
    
    north = BBOX[0]
    south = BBOX[1]
    east = BBOX[2]
    west = BBOX[3]
    
    # gdf = ox.geometries.geometries_from_bbox(BBOX[0], BBOX[1], BBOX[2], BBOX[3], tags)
    
    gdf = ox.geometries.geometries_from_bbox(north, south, east, west, tags)
    
    wWGS84, sWGS84, eWGS84, nWGS84 = gdf.total_bounds
    
    estimatedUTMZone = gdf.estimate_utm_crs(datum_name='WGS 84')
    gdf_proj = gdf.to_crs(estimatedUTMZone)
    # gdf_proj = ox.project_gdf(gdf)  # using UTM coords (alternative)
    
    wUTM, sUTM, eUTM, nUTM = gdf_proj.total_bounds
    
    print('W', west, wWGS84, wUTM)
    print('S', south, sWGS84, sUTM)
    print('E', east, eWGS84, eUTM)
    print('N', north, nWGS84, nUTM)
    
    # gdf_proj = gdf.to_crs(epsg=32630)
    print('  ')
    print(gdf.crs)
    print(gdf_proj.crs)
    print('   ')
    
    return gdf_proj, gdf    # remove one of projected gdfs

#%% ----------------------------------------------

def extractUTMscenery3(fileName):   # read faulty buildings file
    # https://stackoverflow.com/questions/62046624/how-to-save-shapely-geometry-into-a-file-and-load-it-later-into-a-variable
    import pickle
    # polygon from disc
    # with open('faultyBuilds2', 'rb') as resultsFile:
    with open(fileName, 'rb') as resultsFile:
        gdf = pickle.load(resultsFile)
    # gdf = ox.geometries.geometries_from_bbox(BBOX[0], BBOX[1], BBOX[2], BBOX[3], tags)
    gdf_proj = ox.project_gdf(gdf)  # using UTM coords
    return gdf_proj, gdf

#%% -----------------------------------------------------------------

def detectAndReversePoints(points, nPoints):
 
     zerosMat = np.zeros([nPoints,1])
     points3D = np.hstack((points, zerosMat))
     
     alongVects3D = np.empty([nPoints-1,3])
     
     for ii in range(nPoints-1):
         auxx1 = (points3D[ii+1,:]-points3D[ii,:])
         alongVects3D[ii,:] = auxx1/np.linalg.norm(auxx1)
     
     # alongVects2D = alongVects3D[:,:2]
     
     tot = 0
     for ii in range(nPoints-2):
         auxx = np.cross(alongVects3D[ii,:],alongVects3D[ii+1,:])
         print(auxx[2]/np.linalg.norm(auxx))
         tot+=auxx[2]/np.linalg.norm(auxx)   # denominator is magnitude of vector
         
         # mag = np.sqrt(x.dot(x))    # faster calculation of magnitude
     
     if tot < 0:
         # print('+++++++ Reverse order of points ++++++++++++')
         points = np.flip(points, 0)
         
     return points    
 
 #%% --------------------------------------------------------------------
 # from
 # https://gist.github.com/ChrisBarker-NOAA/6fc7de7557856c1351a27df2bb281df5

"""
code for checking if a polygon is cockwise or counter-clockwise

There are two versions:

is_clockwise_convex only works for convex polygons -- but is faster,
it only needs to check three points.

is_clockwise checks all points, but works for convex and cocave
  (note: that's the same as the area calculation)

from:
http://paulbourke.net/geometry/clockwise/

"""

def xrange(x):

    return iter(range(x))


def is_clockwise(poly):
    """
    returns True if the polygon is clockwise ordered, false if not

    expects a sequence of tuples, or something like it (Nx2 array for instance),
    of the points:

    [ (x1, y1), (x2, y2), (x3, y3), ...(xi, yi) ]

    See: http://paulbourke.net/geometry/clockwise/
    """

    total = poly[-1][0] * poly[0][1] - poly[0][0] * poly[-1][1]  # last point to first point
    for i in xrange(len(poly) - 1):
        total += poly[i][0] * poly[i + 1][1] - poly[i + 1][0] * poly[i][1]

    if total <= 0:
        return True
    else:
        return False


def is_clockwise_convex(poly):
    """
    returns True if the polygon is clockwise ordered, False if not

    expects a sequence of tuples, or something like it, of the points:

    [ (x1, y1), (x2, y2), (x3, y3), ...(xi, yi) ]

    This only works for concave polygons. See:

    http://paulbourke.net/geometry/clockwise/
    """

    x0 = poly[0][0]
    y0 = poly[0][1]
    x1 = poly[1][0]
    y1 = poly[1][1]
    x2 = poly[2][0]
    y2 = poly[2][1]

    cp = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1)
    if cp <= 0:
        return True
    else:
        return False
   
 
#%% -----------------------------------------------------------------
 
def fitRotRectangle(points, nPoints):
    
    MINRECT = Polygon(points).minimum_rotated_rectangle
    
    xR,yR = MINRECT.exterior.xy
    xR = np.array(xR)
    yR = np.array(yR)
        
    return MINRECT

#%% -----------------------------------------------------------------

def simplifyPolygonPoints(points, nPoints, minSide):
 
    toBeMerged = []
    for ii in range(nPoints-1):
        v1 = points[ii+1,:] - points[ii,:]
        modv1 = np.linalg.norm(v1)
        if modv1 < minSide:
            toBeMerged.append((ii, ii+1))
             
    casee = -1        
    if ((0,1) in toBeMerged) & (not((nPoints-2,nPoints-1) in toBeMerged)): 
        casee = 0  
    if (not((0,1) in toBeMerged)) & ((nPoints-2,nPoints-1) in toBeMerged):
        casee = 1
    if ((0,1) in toBeMerged) & ((nPoints-2,nPoints-1) in toBeMerged):        
        xAux = (points[0,0]+points[nPoints-2,0]+points[nPoints-1,0])/3
        yAux = (points[0,1]+points[nPoints-2,1]+points[nPoints-1,1])/3
        newPnt2 = np.array([xAux,yAux])
        casee = 2
        
    newPoints = np.empty((0, 2), 'float')
    
    ii = 0    
    while ii < nPoints:
        if (ii,ii+1) in toBeMerged:
            xAux = (points[ii,0]+points[ii+1,0])/2
            yAux = (points[ii,1]+points[ii+1,1])/2
            newPnt = np.array([xAux,yAux])
            newPoints = np.vstack((newPoints, newPnt))
            ii+=1
        else:
            newPoints = np.vstack((newPoints,points[ii,:]))
        ii+=1
                
    if casee == 0:
        newPoints[-1,:] = newPoints[0,:]      
    if casee == 1:
        newPoints[0,:] = newPoints[-1,:]
    if casee == 2:
        newPoints[0,:] = newPnt2
        newPoints[-1,:] = newPnt2 
        
    nNewPoints = len(newPoints)    
    
    return newPoints, nNewPoints

#%% -----------------------------------------------------------------

def simplifyPolygonSides(points, nPoints, minSide, minAnglDeg):

    removePnt = np.empty(0,'int')
    for ii in range(nPoints-2):
        
        x0 = points[ii,:]
        x1 = points[ii+1,:]
        x2 = points[ii+2,:]
        v1 = x1 - x0
        v2 = x2 - x0
        v3 = x2 - x1
        modv1 = np.linalg.norm(v1)
        modv2 = np.linalg.norm(v2)
        modv3 = np.linalg.norm(v3)
        v1hat = v1/modv1
        v2hat = v2/modv2
        v3hat = v3/modv3
        
        angv1v2 = 180/np.pi*np.arccos(np.dot(v1hat,v2hat))
        angv1v3 = 180/np.pi*np.arccos(np.dot(v1hat,v3hat))
        angv2v3 = 180/np.pi*np.arccos(np.dot(v2hat,v3hat))
        
        d = dPointLine(x0,x1,x2)

        if (angv1v2 < minAnglDeg):    
            removePnt = np.append(removePnt, ii+1)
    
    
    removePntList = removePnt.tolist()
    newPoints = np.empty((0, 2), 'float')
    for ii in range(nPoints):
        if ii not in removePntList:
            newPoints = np.vstack((newPoints, points[ii, :]))
    
    points = newPoints  
    
    rowscols = np.shape(points)
    nPoints = rowscols[0]    # no of points in 'points' (last is first repeated)  
            
    return points, nPoints

#%% -----------------------------------------------------------------

def normalsetc(points, nPoints, ref_z):
        
    zerosMat = np.zeros([nPoints,1])
    points3D = np.hstack((points, zerosMat))
    
    outNormals3D = np.empty([nPoints-1,3])
    alongVects3D = np.empty([nPoints-1,3])
    midPoints3D = np.empty([nPoints-1,3])
    
    for ii in range(nPoints-1):
        auxx1 = (points3D[ii+1,:]-points3D[ii,:])
        alongVects3D[ii,:] = auxx1/np.linalg.norm(auxx1)
        auxx = np.cross(alongVects3D[ii,:],ref_z)
        outNormals3D[ii,:] = auxx/np.linalg.norm(auxx)
        midPoints3D[ii,:] = (points3D[ii+1,:] + points3D[ii,:])/2
    
    outNormals2D = outNormals3D[:,:2]
    alongVects2D = alongVects3D[:,:2]
    midPoints2D = midPoints3D[:,:2]
    
    
    # verify by plotting 
    
    # fig, ax = plt.subplots()
    # ax.set_aspect('equal', 'box')
    
    # plt.plot(points[:,0], points[:,1], 'o')
    # plt.plot(points[:,0], points[:,1])
    
    # for ii in range(nPoints):
    #     plt.text(points[ii,0], points[ii,1], ' '+str(ii))     
    # for ii in range(nPoints-1):
    #     plt.text(midPoints2D[ii,0], midPoints2D[ii,1], ' '+str(ii)) 
            
    return outNormals3D, alongVects3D, midPoints3D

#%% -----------------------------------------------------------------

def triangulate(points, OVERALL):

    triLIST = []

    tri = Delaunay(points)  
    
    # create back TRIANGLE POLYGONS to work with SHAPELY
    
    noTri = len(tri.simplices)
    pointerr = [];    
    triList = []
    centroidList = []
    centroidCoords = np.empty([0,2])
    for ii in range(noTri):
        pntr = tri.simplices[ii]
        PLYN = Polygon([
            (points[pntr[0],0], points[pntr[0],1]),
            (points[pntr[1],0], points[pntr[1],1]),
            (points[pntr[2],0], points[pntr[2],1]),
            (points[pntr[0],0], points[pntr[0],1])
            ])
        triList.append(PLYN)
        centroidList.append(PLYN.centroid)
        centroidCoords = np.vstack((centroidCoords, np.array(list(PLYN.centroid.coords))))
        # print(ii)
        auux = PLYN.within(OVERALL)
        if auux == False:
            pointerr.append(ii)
        else:
           triLIST.append(PLYN) 
        # print( PLYN.centroid)
    
    # plotting centroids as well
    
    
    # plt.triplot(points[:,0], points[:,1], tri.simplices)
    # plt.plot(points[:,0], points[:,1], 'o')
    # plt.plot(centroidCoords[:,0],centroidCoords[:,1],'+')
    # for ii in range(noTri):
    #     plt.text(centroidCoords[ii,0], centroidCoords[ii,1], ' '+str(ii))
    # plt.show()
    
    
    #
    
    newTriList = []
    newSimplices = np.empty((0,3))
    currSimplices = tri.simplices
    
    for ii in range(noTri):
        if not(np.isin(ii, pointerr)):
            newTriList.append(triList[ii])
            newSimplices = np.vstack((newSimplices, currSimplices[ii,:]))
    
    # triList = newTriList
    # del(newTriList)
    
    newSimplices = newSimplices.astype(np.int32)
    tri.simplices = newSimplices
    
    # plot again
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal', 'box')
    plt.triplot(points[:,0], points[:,1], newSimplices)
    plt.plot(points[:,0], points[:,1], 'o')
    plt.show()
        
    return triLIST
           
#%% ------------------------------------------------
# from  https://stackoverflow.com/questions/21824157/how-to-extract-interior-polygon-coordinates-using-shapely

def extract_poly_coords(geom):
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
        interior_coords = []
        for interior in geom.interiors:
            interior_coords += interior.coords[:]
    elif geom.type == 'MultiPolygon':
        exterior_coords = []
        interior_coords = []
        for part in geom:
            epc = extract_poly_coords(part)  # Recursive call
            exterior_coords += epc['exterior_coords']
            interior_coords += epc['interior_coords']
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))
    return {'exterior_coords': exterior_coords,
            'interior_coords': interior_coords}

#%% ----------------------------------------------------
def createActualTriangles(points, triangles):

    triLIST = []
    
    nTri, xxx = np.shape(triangles)
    
    for ii in range(nTri) :
        P0 = triangles[ii,0]
        P1 = triangles[ii,1]
        P2 = triangles[ii,2]
        
        PLYN = Polygon([
                (points[P0,0], points[P0,1]),
                (points[P1,0], points[P1,1]),
                (points[P2,0], points[P2,1]),
                (points[P0,0], points[P0,1])
                ])
        
        triLIST.append(PLYN)
        
    return triLIST

#%% ----------------------------------------------------

def isItCW(pointsIN):
    
    # from https://stackoverflow.com/questions/61757304/how-to-determine-if-a-point-lies-inside-a-concave-hull-alpha-shape
    
    # import matplotlib.pyplot as plt
    # from descartes import PolygonPatch
    # import alphashape
    import random
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon
    
    nPointsIN, xxx = np.shape(pointsIN)

    Flagg = True
    
    pointsOUT = pointsIN[0:-1,:]
    
    minX = min( pointsOUT[:,0])
    maxX = max( pointsOUT[:,0])
    minY = min( pointsOUT[:,1])
    maxY = max( pointsOUT[:,1])

    # points = list(map(tuple, pointsIN))
    points = list(map(tuple, pointsOUT))
    
    # # points = [(0., 0.), (0., 1.), (1., 1.), (1., 0.),
    # #           (0.5, 0.25), (0.5, 0.75), (0.25, 0.5), (0.75, 0.5)]

    # # alpha_shape = alphashape.alphashape(points, 2.0) # Create the alpha shape
    # alpha_shape = alphashape.alphashape(points) # Create the alpha shape

    # # Plotting the alpha shape over the input data
    # fig, ax = plt.subplots()
    # ax.set_aspect('equal', 'box')
    # ax.scatter(*zip(*points), c='green')
    # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))

    # N = 1000 # number of random points
    # for i in range(N):
    #     xx = round(random.uniform(minX, maxX), 2)
    #     yy = round(random.uniform(minY, maxY), 2)
    #     pnt = Point(xx,yy) # analysis point
    #     print(alpha_shape.contains(pnt))
    #     if alpha_shape.contains(pnt) == True:
    #         plt.scatter(xx,yy,c='blue')
    #         break
    #     else:
    #         plt.scatter(xx,yy,c='red')
    
    
    # insidePoint = np.empty([1, 2])  
    # insidePoint[0,0]=pnt.x
    # insidePoint[0,1]=pnt.y
    # insidePoint = insidePoint[0,:]
    
    # alternative find point inside polygon ------
    # from https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python

    Pol = Polygon(points)
    stopWhile = False
    while stopWhile == False:
        xx = round(random.uniform(minX, maxX), 2)
        yy = round(random.uniform(minY, maxY), 2)
        pnt = Point(xx,yy) # analysis point
        stopWhile = Pol.contains(pnt)
 
 
 
    insidePoint = np.empty([1, 2])  
    insidePoint[0,0]=pnt.x
    insidePoint[0,1]=pnt.y
    insidePoint = insidePoint[0,:]
    
    
    # compute angles
    
    totSUMpos = 0
    totSUMneg = 0
    for ii in range(nPointsIN-1):
        # print(ii)
        # print(pointsIN[ii,:])
        # print(pointsIN[ii+1,:])
        v0 = pointsIN[ii,:] - insidePoint
        v1 = pointsIN[ii+1,:] - insidePoint
        mag0 = np.sqrt(v0.dot(v0))
        mag1 = np.sqrt(v1.dot(v1))
        v0n = v0/mag0
        v1n = v1/mag1
        
        v0n3D = np.hstack([v0n, 0])
        v1n3D = np.hstack([v1n, 0])
        
        
        xprod = np.cross(v0n3D,v1n3D)
        
        ang = np.degrees(np.arcsin(xprod[2]))
        # print(ang)
        if ang >= 0:
            totSUMpos =totSUMpos + ang
        else:
            totSUMneg =  totSUMneg + ang
                
    if abs(totSUMpos) > abs(totSUMneg):
        Flagg = False
    elif abs(totSUMpos) < abs(totSUMneg):
        Flagg = True
    else:
        print('something is very wrong')
    return Flagg

#%% ---------------------------------------------------------------------

def polygonFromArrayArrays(x,y):
    
    # based on https://www.matecdev.com/posts/shapely-polygon-from-points.html
    
    import shapely.geometry as sg
    
    points = []
    for ii in range(len(x)):
        currPnt = sg.Point(x[ii], y[ii])
        points.append(currPnt)
        
    out_poly = sg.Polygon([[p.x, p.y] for p in points])
    
    return out_poly














def polygonFrom2array_arrays(x,y):
    # x and y are two array.array same length inputs containing 
    # the output closed polygon vertices
    
    import shapely.geometry as sg
    
    points = []
    for ii in range(len(sx)):
        currPnt = sg.Point(x[ii], y[ii])
        points2.append(currPnt)
        
    poly2 = sg.Polygon([[p.x, p.y] for p in points2])
    
    
    
    
    
    return poly_out