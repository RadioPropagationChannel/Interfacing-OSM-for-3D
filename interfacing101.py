# from IPython import get_ipython;   
# get_ipython().magic('reset -sf')    

# import osmnx as ox
# from IPython.display import Image
import numpy as np
# from shapely.geometry import Point
from shapely.geometry import Polygon
# from shapely.geometry import mapping
from shapely.geometry import MultiPoint
import pandas as pd

from scipy.spatial import ConvexHull, convex_hull_plot_2d

import matplotlib.pyplot as plt
# from scipy.spatial import ConvexHull
# from scipy.ndimage import rotate
# from scipy.spatial import Delaunay

from osgeo import gdal

# import pandas as pd
# import geopandas as gpd
# import pickle

# import math
# import time
# import sys

# import concavity

import contextily as cx

import os
# %matplotlib inline

# import my help functions and constants and inputs

import various_functions as vf
import constants_and_inputs as cc
import triangulation_John_Burkardt as tr
import plotting3D as pl3
# import baseMaps as bM

# -----------------------------------------------------------------------
    
def main():
        
    plt.close('all')    # Close all open figures
    print()
    print('---------------------------------------------------')
    print('----------- INTERFACING WITH OSM DATA -------------')
    print('---------------------------------------------------')
    print()
    
    gdf_proj, gdf, CENTER = vf.extractUTMscenery2(cc.BBOX, cc.tags)
    gdf_proj = vf.filterNonPolygonsOut(gdf_proj)
    
    manualDeleteList = [174]
    gdf_proj = vf.removePolygonManually(gdf_proj, manualDeleteList)
     
    crsInUse = gdf_proj.crs   # to be used later on when loading baseMap
    print('No. of builings originally is',gdf_proj.shape[0])
    
    sigFile = str(cc.BBOX[0])+str(cc.BBOX[1])+str(cc.BBOX[2])+str(cc.BBOX[3])
    
    sigFile_Contained = sigFile+'_Contained_'+'.npy'
    if os.path.exists(sigFile_Contained) == False:
        Contained = vf.removeTotContainedPols(gdf_proj)     
        with open(sigFile_Contained, 'w') as f:
            np.save(sigFile_Contained, Contained, allow_pickle=True, fix_imports=True)       
            f.close()
    else:
        Contained = np.load(sigFile_Contained)
             
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box')    
    
    List1 = list(Contained[:,0])
    Counter = int(0)    
    for ii in List1:
        currPOLY = gdf_proj.iloc[ii].geometry
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        currCentr = np.array(currPOLY.centroid.coords.xy)
        AXIS_1.plot(currCentr[0],currCentr[1], 'o')
        AXIS_1.text(currCentr[0],currCentr[1], str(ii))
        
        jj = Contained[Counter,1]
        currPOLY = gdf_proj.iloc[jj].geometry
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        currCentr = np.array(currPOLY.centroid.coords.xy)
        AXIS_1.plot(currCentr[0],currCentr[1], 'o')
        AXIS_1.text(currCentr[0],currCentr[1], str(jj))
        Counter += 1
       
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box')   
    
    List1 = list(Contained[:,0])
    Counter = int(0)    
    for ii in List1:
        currPOLY = gdf_proj.iloc[ii].geometry
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        currCentr = np.array(currPOLY.centroid.coords.xy)
        AXIS_1.plot(currCentr[0],currCentr[1], 'o')
        AXIS_1.text(currCentr[0],currCentr[1], str(ii))
                
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box')   
    
    List1 = list(Contained[:,0])
    Counter = int(0)    
    for ii in List1:
         jj = Contained[Counter,1]
         currPOLY = gdf_proj.iloc[jj].geometry
         x,y = currPOLY.exterior.xy
         plt.plot(x,y)
         currCentr = np.array(currPOLY.centroid.coords.xy)
         AXIS_1.plot(currCentr[0],currCentr[1], 'o')
         AXIS_1.text(currCentr[0],currCentr[1], str(jj))
         Counter += 1    
    
    List2 = list(Contained[:,1])   # removal list
    
    gdf_proj = gdf_proj.drop(gdf_proj.index [List2])
    print('No. of builings after removing duplicates is',gdf_proj.shape[0])
    
    # note the original gdf_proj may have changed size
            
    # OLMatrix is overlay matrix where boxes touching ared detected
        
    sigFile_Groupings = sigFile+'_Groupings_'+'.npy'
    if os.path.exists(sigFile_Groupings) == False:
        OLMatrix = vf.overlayMatrix(gdf_proj)     
        with open(sigFile_Groupings, 'w') as f:
            np.save(sigFile_Groupings, OLMatrix, allow_pickle=True, fix_imports=True)       
            f.close()
    else:
        OLMatrix = np.load(sigFile_Groupings)
        
    nPolys = OLMatrix.shape[0]   

    indxs = np.argwhere(OLMatrix == 1)
    
    groupings = vf.groupBoxes(indxs)
    
        
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box')    
    for ii in range(nPolys):
        currPOLY = gdf_proj.iloc[ii].geometry
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        currCentr = np.array(currPOLY.centroid.coords.xy)
        AXIS_1.plot(currCentr[0],currCentr[1], 'o')
        AXIS_1.text(currCentr[0],currCentr[1], str(ii))
        
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box')    
    for ii in range(nPolys):
        currPOLY = gdf_proj.iloc[ii].geometry
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        currCentr = np.array(currPOLY.centroid.coords.xy)
        Location = vf.findInlistofArrays(ii,groupings)
        AXIS_1.plot(currCentr[0],currCentr[1], 'o')
        AXIS_1.text(currCentr[0],currCentr[1], str(Location))    
           
    print() 
    print('------------- reading no of building levels -----------------')
    print()
    # deal with building levels / floor -----------------------------  
    building_levels = np.array(gdf_proj['building:levels']) # no of floors
    print('Building levels:')
    print(building_levels)
    nLevels = len(building_levels)
    # levels_x = np.linspace(0, nLevels-1, nLevels, 'int')
    
    for ii in range(nLevels):
        if type(building_levels[ii])==float:
            building_levels[ii] = 'E' # empty 
    
    import random    
    new_building_levels = np.zeros(nLevels, int)
    for ii in range(nLevels):
        if building_levels[ii] == 'E':
            new_building_levels[ii] = random.randint(1, 6) # no of storeys 
        else:
            new_building_levels[ii] = int(building_levels[ii]) # no of storeys 
            
    building_heights =  np.zeros(nLevels, float)       
    # convert storeys to heights
    for ii in range(nLevels):
        building_heights[ii] = 3.5 + 3*new_building_levels[ii] + 1.5 
              
    #%% loop over all buildings
       
    LENGTHS = []
    ROTATIONS = []
    ALLBBOXES = []
    ALLTRIANGLES = []
    ALLSIMPLPOLYGONS = []
    ALLINTERNALS = []
    print() 
    print('------------- POLYGON LOOP -----------------')
    print()
    for ii in range(len(gdf_proj)):  # loop
           
        print(f'-------------- NEW POLYGON ---------------no. {ii}')
        
        aa = gdf_proj.geometry.iloc[ii]
        resultt = vf.extract_poly_coords(aa)
                
        points = np.array(resultt['exterior_coords'])        
        nPoints, xxx = np.shape(points)
        
        INTPOLYS = []
        nHoles = len(aa.interiors)
        ind2 = int
        curPolyArray = np.empty((0,2))
        if nHoles > 0:
            
            print(f'number of holes : {nHoles}')
            internalPoints = np.array(resultt['interior_coords'])
            for jj in range(nHoles):
                curPolyArray = np.empty((0,2))
                
                if jj == 0:
                    firstPoint = internalPoints[0,:]
                else:
                    firstPoint = internalPoints[ind2+1,:]
                           
                testVar = abs(internalPoints - firstPoint)
                testVar2 = testVar < 1e-5
                indexx = np.where(testVar2 == [True, True])
                indexx2 = indexx[0]
                if len(indexx2) > 4:
                    print('problem with internal hole(s) in building {ii}')
                    continue
                else:
                    ind1 = indexx2[0]
                    ind2 = indexx2[3]
                    curPolyArray = np.vstack((curPolyArray, internalPoints[ind1:ind2+1,:]))
                    INTPOLYS.append(curPolyArray)        
        
        ALLINTERNALS.append(INTPOLYS)
        
        #%% detect CW or CCW
        
        Flagg = vf.isItCW(points)
        if Flagg == True:               # triangulation works with CCW concave polygons
            print('Polygon no. '+str(ii)+' needs rotating')
            points = np.flip(points, 0) # reverse sense of rotation if necessary to CCW               
        
        #%% bounding box calculation
        
        points = vf.centerBox(points, CENTER)  # from here ALLBBOXES are centered at origin
        nPoints = points.shape[0]
        MINRECT, Rot, Len = vf.fitRotRectangle(points, nPoints)    
        ALLBBOXES.append(MINRECT) 
        ROTATIONS.append(Rot)
        LENGTHS.append(Len)
        
        #%% Simplify original
                
        points,nPoints = vf.simplifyPolygonPoints(points, nPoints, cc.minSide) 
        
        points, nPoints = vf.simplifyPolygonSides(points, nPoints, cc.minSide, cc.minAnglDeg)
        
        # points = vf.centerBox(points, CENTER)  # from here on OVERALL are centered at origin
        
        auxOverall = points.tolist()
        OVERALL = Polygon(auxOverall) # reconstruct overall simplified, CCW, polygon  
        
        ALLSIMPLPOLYGONS.append(OVERALL)
            
        # %% Apply alternative triangulation 

        x = points[0:-1,0]     # remove last point since it is equal to first
        y = points[0:-1,1]
       
        n = np.size(y)
                    
        print ( '' )
        print ( 'polygon_triangulation .......' )
        print ( '' )

        triangles = tr.polygon_triangulate (n, x, y)
               
        triLIST0 = vf.createActualTriangles(points, triangles)
        ALLTRIANGLES.append(triLIST0)
        
        ('End of loop all polygos go ahead and plot ......................')
    
            
    cx.add_basemap(AXIS_1, crs = gdf_proj.crs);        
    
    
    fig, ax = plt.subplots(1, figsize=(9, 9))
    ax.set_aspect('equal', 'box')    
    plt.title('Rectable fitted builings')
    plt.ylabel('Northing wrt center of scene, m')
    plt.xlabel('Easting wrt center of scene, m')
    for ii in range(len(ALLBBOXES)):  # fit.-rot. rectangle
        currPOLY = ALLBBOXES[ii]
        
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        
        currCentr = np.array(currPOLY.centroid.coords.xy)
        plt.plot(currCentr[0],currCentr[1], 'o')
        plt.text(currCentr[0],currCentr[1], str(ii))
    plt.show() 
    
    
    fig, ax = plt.subplots(1, figsize=(9, 9))
    ax.set_aspect('equal', 'box')    
    plt.title('Groupings')
    plt.ylabel('Northing wrt center of scene, m')
    plt.xlabel('Easting wrt center of scene, m')
    for ii in range(len(ALLBBOXES)):  # fit.-rot. rectangle
        currPOLY = ALLBBOXES[ii]
        
        x,y = currPOLY.exterior.xy
        plt.plot(x,y)
        
        currCentr = np.array(currPOLY.centroid.coords.xy)
        Location = vf.findInlistofArrays(ii,groupings)
        plt.plot(currCentr[0],currCentr[1], 'o')
        plt.text(currCentr[0],currCentr[1], str(Location))
    plt.show() 
    
    
    fig, ax = plt.subplots(1, figsize=(9, 9))
    ax.set_aspect('equal', 'box')   
    for ii in range(len(ALLTRIANGLES)):
        currPOLY = ALLTRIANGLES[ii]  
        for jj in range(len(currPOLY)):      
            x,y = currPOLY[jj].exterior.xy
            plt.plot(x,y)     
            ctrid = currPOLY[jj].centroid
            plt.plot(ctrid.xy[0], ctrid.xy[1],'+')
            # plt.text(ctrid.x, ctrid.y, str(jj))         ##### does not work
            # plt.text([ctrid.x], [ctrid.y], str(jj))         ##### does not work
        currBOX = ALLBBOXES[ii]
        currCentr = np.array(currBOX.centroid.coords.xy)
        plt.plot(currCentr[0],currCentr[1], 'o')
        plt.text(currCentr[0],currCentr[1], str(ii))
    plt.show() 

       
    # now introduce polygon(triangle) height
    # and create 3D geometry
    # and plot
    
    # from: https://stackoverflow.com/questions/67410270/how-to-draw-a-flat-3d-rectangle-in-matplotlib
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
        
    for ii in range(len(ALLTRIANGLES)):   #(4,6): #
        currPOLY = ALLTRIANGLES[ii]  
        Z = np.zeros([len(currPOLY)*3,3])
        Z2 = np.zeros([len(currPOLY)*3,3])
        verts = [];
        for jj in range(len(currPOLY)):      
            x,y = currPOLY[jj].exterior.xy # x and y contain trangles with first vertice repeated
            
            Z[jj*3,:] = [x[0], y[0], 0.0]
            Z[jj*3+1] = [x[1], y[1], 0.0]
            Z[jj*3+2] = [x[2], y[2], 0.0]
            
            Z2[jj*3,:] = [x[0], y[0], building_heights[ii]]
            Z2[jj*3+1] = [x[1], y[1], building_heights[ii]]
            Z2[jj*3+2] = [x[2], y[2], building_heights[ii]]
            
            verts.append([Z[jj*3,:], Z[jj*3+1,:], Z[jj*3+2,:]])
            verts.append([Z2[jj*3,:], Z2[jj*3+1,:], Z2[jj*3+2,:]])  
            verts.append([Z[jj*3,:], Z[jj*3+1,:], Z2[jj*3+1,:], Z2[jj*3,:]]) 
            verts.append([Z[jj*3+1,:], Z[jj*3+2,:], Z2[jj*3+2,:], Z2[jj*3+1,:]]) 
            verts.append([Z[jj*3+2,:], Z[jj*3,:], Z2[jj*3,:], Z2[jj*3+2,:]])
            
            # end of single polygon loop  --------------
        
        Z = np.vstack((Z,Z2))   # contains bottom and top vertices from currPOLY
        
        ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2], zdir='z', s=0, c=None, depthshade=True)  #plots vertices of currPOLY
        ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=0.1, edgecolors='k', alpha=0.10))
        
        
        currBOX = ALLBBOXES[ii]
        currCentr = np.array(currBOX.centroid.coords.xy)
        ax.scatter3D(currCentr[0],currCentr[1], building_heights[ii],'o', s=4)
        ax.text3D(float(currCentr[0]),float(currCentr[1]), building_heights[ii], str(ii)) 
          
    pl3.set_axes_equal(ax)
    plt.show()   
    
    # plot the simplified scenario in 3D    --------------------------
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
        
    for ii in range(len(ALLSIMPLPOLYGONS)):   #(4,6): #
        currPOLY = ALLSIMPLPOLYGONS[ii]  
        x,y = currPOLY.exterior.xy
        xx = np.array(x)
        yy = np.array(y)
        l = np.shape(xx)
        lenPol = l[0] - 1   # last point is repeated
        Z1 = np.zeros([lenPol,3])
        Z2 = np.zeros([lenPol,3])
        Z = np.zeros([lenPol*2,3])
        verts = [];                 # list of vertices for plotting sidewalls
        for jj in range(lenPol):    # this loop should build the side faces index jj is face #
                        
            Z1[jj,:] = [xx[jj], yy[jj], 0.0]            
            Z2[jj,:] = [xx[jj], yy[jj], building_heights[ii]]
            
            if jj == lenPol - 1:
                verts.append([Z1[jj,:], Z1[0,:], Z2[0,:], Z2[jj,:]]) 
            else:
                verts.append([Z1[jj,:], Z1[jj+1,:], Z2[jj+1,:], Z2[jj,:]]) 
            
            # end of single polygon loop  --------------
        
        Z = np.vstack((Z1,Z2))   # contains bottom and top vertices from currPOLY
        
        ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2], zdir='z', s=0, c=None, depthshade=True)  #plots vertices of currPOLY
        ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=0.1, edgecolors='k', alpha=0.10))
        
        currBOX = ALLBBOXES[ii]
        currCentr = np.array(currBOX.centroid.coords.xy)
        ax.scatter3D(currCentr[0],currCentr[1], building_heights[ii],'o', s=4)
        ax.text3D(float(currCentr[0]),float(currCentr[1]), building_heights[ii], str(ii)) 
          
    pl3.set_axes_equal(ax)
    plt.show()      
    
    ROTALLBBOXES = []    
    for ii in range(len(ALLBBOXES)):        
        currPolygon = ALLBBOXES[ii]
        if ROTATIONS[ii][0] < 0:
            Rotation = 180 +  ROTATIONS[ii][0]
        else:
            Rotation = ROTATIONS[ii][0] 
        rotAngleDeg = 90 - Rotation
        new_points, rotPolygon = vf.rotatePolygon(currPolygon, rotAngleDeg)
        ROTALLBBOXES.append(rotPolygon)
        
    # ----- NEW ALGOR: GROUPING ORIENTED ::::::::::::::::::::::::::::::::::::: 
    print() 
    print('------------- GROUPINGS LOOP -----------------')
    print()    
    
    
    ALLGROUPS = []
    
    for currGrouping in range(len(groupings)): 
        
        if currGrouping >= 0:
            print(f'-------------- NEW GROUPING ---------------no. {currGrouping}')
            # plt.close('all')
            RotsInGroup, LensInGroup, perimeterPointsPlot, groupCenter = \
                vf.groupingStep1(ALLBBOXES, ROTATIONS, LENGTHS, groupings, currGrouping)
        
            # do the processing below HERE
            
            LOCALROTATED, refRot = vf.groupingStep2rotation(groupCenter, perimeterPointsPlot, \
                ALLBBOXES, ROTATIONS, LensInGroup, RotsInGroup, groupings, currGrouping)
            
     
            # Group Reference direction into vector form
            x_vect = np.array([np.cos(refRot*np.pi/180), np.sin(refRot*np.pi/180)])
            y_vect = np.array([np.cos((refRot+90)*np.pi/180), \
                                   np.sin((refRot+90)*np.pi/180)])
            
            if len(LOCALROTATED) > 0:
            
                matrixDist = vf.distMatrix(LOCALROTATED)
                
                SQREDperiemeter = vf.groupingStep3trim(x_vect, y_vect, groupCenter, LOCALROTATED, matrixDist, \
                                      perimeterPointsPlot, groupings, currGrouping)
                    
                ListofRects, ListofRectsPolygs = \
                    vf.GroupingExtractRects(x_vect, y_vect, SQREDperiemeter, groupCenter, currGrouping)  
                    
                ALLGROUPS.append(ListofRectsPolygs)       
                
                # vf.GroupVerifications(ListofRectsPolygs, x_vect, y_vect, perimeterPointsPlot, groupings, currGrouping)  
            
    # END OF NEW ALGOR: GROUPING ORIENTED :::::::::::::::::::::::::::::::::::
    
        
    for ii in range(len(groupings)):   # add individual isolated builings
        if len(groupings[ii]) == 1:
            pointt = groupings[ii][0]
            oneElementList = [ALLBBOXES[pointt]]
            ALLGROUPS.append(oneElementList)

    SHIFTEDGROUPS = vf.plotAllGroupoings(ALLGROUPS)
    
    exportData = vf.extractAndExport(SHIFTEDGROUPS)
    print(exportData)
    
if __name__ == '__main__':
    main()
    