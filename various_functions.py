import numpy as np
import osmnx as ox
import matplotlib.pyplot as plt
# from shapely.geometry import MultiPoint
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely import affinity

import geopandas as gpd
import pandas as pd

from math import sin, cos, radians
# from graphics import *

import array as arr

from scipy.spatial import ConvexHull, convex_hull_plot_2d

import random 

#%% -----------------------------------------------------------------
#%% extract and export to file
def extractAndExport(SHIFTEDGROUPS):
    exportData = []
    for ii in range(len(SHIFTEDGROUPS)):
        currGrList = []
        noStoryes = random.randint(3, 6) # no of storeys 
        HHz = 3.5 + 3*noStoryes + 1.5
        for jj in range(len(SHIFTEDGROUPS[ii])):
            currPoly = SHIFTEDGROUPS[ii][jj]
            center = np.array([currPoly.centroid.xy[0][0], \
                               currPoly.centroid.xy[1][0], 0.0])
            P0 = np.array([currPoly.exterior.coords.xy[0][0], currPoly.exterior.coords.xy[1][0]])
            P1 = np.array([currPoly.exterior.coords.xy[0][1], currPoly.exterior.coords.xy[1][1]])
            P2 = np.array([currPoly.exterior.coords.xy[0][2], currPoly.exterior.coords.xy[1][2]])
            v0 = P1 - P0
            v1 = P2 - P1
            side0 = np.linalg.norm(v0)
            side1 = np.linalg.norm(v1)
            if side0 >= side1:
                LLx = side0
                WWy = side1
                refVect = v0/side0
            else:
                LLx = side1
                WWy = side0
                refVect = v1/side1
            Rot = np.arctan2(refVect[1], refVect[0])*180/np.pi    
            currGrList.append([center[0], center[1], center[2],LLx,WWy, HHz, Rot])    
        exportData.append(currGrList)
    
    file = open('export.txt','w')       # save to text file --------
    for ii in range(len(exportData)):
        for jj in range(len(exportData[ii])):
            # file.write(str(exportData[ii][jj])+'\n')
            for kk in range(6):
                file.write(str((round(exportData[ii][jj][kk],2)))+'  ')
            file.write(str((round(exportData[ii][jj][6],2)))+'\n')
    file.close()
       

        
    # and finally, plot in 3D
    
    # from: https://stackoverflow.com/questions/67410270/how-to-draw-a-flat-3d-rectangle-in-matplotlib
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for ii in range(len(SHIFTEDGROUPS)):
        for jj in range(len(SHIFTEDGROUPS[ii])):
            currPOLY = SHIFTEDGROUPS[ii][jj]  
            currH = exportData[ii][jj][5] 
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
                Z2[jj,:] = [xx[jj], yy[jj], currH]
                
                if jj == lenPol - 1:
                    verts.append([Z1[jj,:], Z1[0,:], Z2[0,:], Z2[jj,:]]) 
                else:
                    verts.append([Z1[jj,:], Z1[jj+1,:], Z2[jj+1,:], Z2[jj,:]]) 
                
                # end of single polygon loop  --------------
            
            Z = np.vstack((Z1,Z2))   # contains bottom and top vertices from currPOLY
            
            ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2], zdir='z', s=0, c=None, depthshade=True)  #plots vertices of currPOLY
            ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=0.1, edgecolors='k', alpha=0.10))
  













    
    

    return exportData
#%% remove manuelly from list
def removePolygonManually(gdf, List): 
    gdf.drop(gdf.index[List], inplace=True)    
    return gdf
#%% test distance exceeds min limit
def exceeDistTest(group, minDist):
    passedTest = True
    
    for ii in range(len(group)):
        curPoly = group[ii]
        for jj in range(ii+1,len(group)):
            nextPoly = group[jj]
            dd = curPoly.distance(nextPoly)
            if dd < minDist:
                break
        if dd < minDist:
            passedTest = False
            break  
    
    return passedTest
#%% Group overlapping
def overlapTest(group):
    Overlap = False
    
    # DisjointMatrix = np.full((len(group),len(group)), True)
    for ii in range(len(group)):
        curPoly = group[ii]
        for jj in range(ii+1,len(group)):
            nextPoly = group[jj]
            result = curPoly.disjoint(nextPoly)
            if result == False:
                break
        if result == False:
            Overlap = True
            break

    return Overlap
#%%
def plotAllGroupoings(ALLGROUPS):

    overallCentroid = np.array([0.0,0.0])
    Centroids = []
    convexGroupPerimeters = []
    convexGroupPerimetersShply = []
    for ii in range(len(ALLGROUPS)):
        pointsArray = np.zeros((len(ALLGROUPS[ii]*4),2))
        for jj in range(len(ALLGROUPS[ii])):
            currPoly = ALLGROUPS[ii][jj]
            points = currPoly.exterior.coords.xy
            x = points[0]
            y = points[1]
            for kk in range(4):
                pointsArray[4*jj+kk,]=np.array([x[kk],y[kk]])
        Perimeter = ConvexHull(pointsArray)
        pointter = Perimeter.vertices
        newPerimeter = pointsArray[pointter,:]
        currCentroid = sum(newPerimeter)/newPerimeter.shape[0]
        Centroids.append(currCentroid)
        newPerimeter = np.vstack((newPerimeter, newPerimeter[0,:]))  # repeat first point to create polygon
        convexGroupPerimeters.append(newPerimeter)   
        newPerimeterShply = Polygon(newPerimeter)
        convexGroupPerimetersShply.append(newPerimeterShply)
        
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box') 
    plt.title('Overall scene') 
    for ii in range(len(ALLGROUPS)):
        currCentroid = Centroids[ii]
        plt.plot(currCentroid[0],currCentroid[1],'o')
        plt.text(currCentroid[0],currCentroid[1],str(ii))
        for jj in range(len(ALLGROUPS[ii])):
            currPoly = ALLGROUPS[ii][jj]
            points = currPoly.exterior.coords.xy
            x = points[0]
            y = points[1]
            plt.plot(x,y)
        currHull = convexGroupPerimeters[ii]
        plt.plot(currHull[:, 0], currHull[:, 1], 'k-')
        plt.plot(overallCentroid[0],overallCentroid[1],'*')
    plt.grid(which='major', axis='both')  
    
        
    # verify superpositions and make measurements: dinstances, ...
    
    DisjointMatrix = np.full((len(ALLGROUPS),len(ALLGROUPS)), True)
    centroidDistances = np.zeros((len(ALLGROUPS),len(ALLGROUPS)))
    distToCenter = np.zeros((len(ALLGROUPS)))
    for ii in range(len(ALLGROUPS)):
        curPoly = convexGroupPerimetersShply[ii]
        currCentroid = Centroids[ii]
        distToCenter[ii] = np.linalg.norm(currCentroid)
        for jj in range(ii+1,len(ALLGROUPS)):
            nextPoly = convexGroupPerimetersShply[jj]
            DisjointMatrix[ii,jj] = curPoly.disjoint(nextPoly)
            DisjointMatrix[jj,ii] = curPoly.disjoint(nextPoly)
            nextCentroid = Centroids[jj]
            currDist = np.linalg.norm(currCentroid - nextCentroid)
            centroidDistances[ii,jj] = currDist
            
    print(np.any(np.invert(DisjointMatrix)))

    Overlap = overlapTest(convexGroupPerimetersShply)    
    
    # shift perimeters ------------------------------------------
    shftdPeri = convexGroupPerimetersShply       
    maxDist = np.max(centroidDistances)
    OrderOfCentroids = np.argsort(distToCenter)        
    refCentrPointr = OrderOfCentroids[0]
    
    xyGrTranslations = np.zeros((len(ALLGROUPS),2))
    
    shiftFactor = 10   # in meters  # shift parameter
    minDist = 10      # m specifed min distance between perimeters
    Exponent = 1.5
    
    allDisjoint = False
    while allDisjoint == False:
          
        RefCentroid = Centroids[refCentrPointr]
        translatedSet = []
        for ii in range(len(ALLGROUPS)):
            if ii != refCentrPointr:
                currCentroid = Centroids[ii]     
                transVect =  currCentroid - RefCentroid
                lenTransVect = np.linalg.norm(transVect)
                transVectn = transVect/lenTransVect
                shftX = transVectn[0]*shiftFactor*(lenTransVect/maxDist)**Exponent
                shftY = transVectn[1]*shiftFactor*(lenTransVect/maxDist)**Exponent
                auxx = \
                    affinity.translate(shftdPeri[ii], shftX, shftY)
                translatedSet.append(auxx)
                xyGrTranslations[ii,] = xyGrTranslations[ii,] + np.array([shftX, shftY]) 
                # if ii == 0:
                #     print(xyGrTranslations[ii,])
            else:
                translatedSet.append(shftdPeri[ii])
         
        
        # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
        # AXIS_1.set_aspect('equal', 'box') 
        # plt.title('Overall scene')
            
        # for ii in range(len(translatedSet)):
        #     currPoly = translatedSet[ii]
        #     points = currPoly.exterior.coords.xy
        #     x = points[0]
        #     y = points[1]
        #     plt.plot(x,y,'k')
        # plt.grid(which='major', axis='both')  

        shftdPeri = translatedSet  # supply shifted group to next iteration
        
        Overlap = overlapTest(translatedSet)
        
        if Overlap == False: 
            passedTest = exceeDistTest(translatedSet, minDist)
            if passedTest == True: 
                allDisjoint = True
                
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box') 
    plt.title('Overall scene')
        
    for ii in range(len(translatedSet)):
        currPoly = translatedSet[ii]
        points = currPoly.exterior.coords.xy
        x = points[0]
        y = points[1]
        plt.plot(x,y,'k')
    plt.grid(which='major', axis='both')  
    
    SHIFTEDGROUPS = []
    for ii in range(len(ALLGROUPS)):
        shftX, shftY = xyGrTranslations[ii,]
        currShiftedGroup = []
        for jj in range(len(ALLGROUPS[ii])):
            currPoly = ALLGROUPS[ii][jj]
            auxx = affinity.translate(currPoly, shftX, shftY)
            currShiftedGroup.append(auxx)
        SHIFTEDGROUPS.append(currShiftedGroup)
        
    
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box') 
    plt.title('Overall scene') 
    for ii in range(len(SHIFTEDGROUPS)):
       # currCentroid = Centroids[ii]
       # plt.plot(currCentroid[0],currCentroid[1],'o')
       # plt.text(currCentroid[0],currCentroid[1],str(ii))
       for jj in range(len(SHIFTEDGROUPS[ii])):
           currPoly = SHIFTEDGROUPS[ii][jj]
           points = currPoly.exterior.coords.xy
           x = points[0]
           y = points[1]
           plt.plot(x,y)
       # currHull = convexGroupPerimeters[ii]
       # plt.plot(currHull[:, 0], currHull[:, 1], 'k-')
       # plt.plot(overallCentroid[0],overallCentroid[1],'*')
    plt.grid(which='major', axis='both')   

    
    return SHIFTEDGROUPS
#%%
def GroupVerifications(ListofRectsPolygs, x_vect, y_vect, perimeterPointsPlot, groupings, currGrouping)  :
    
    
    fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    AXIS_1.set_aspect('equal', 'box') 
    plt.title('Current grouping: '+ str(currGrouping))
    for ii in range(len(ListofRectsPolygs)):
        currPoly = ListofRectsPolygs[ii]
        points = currPoly.exterior.coords.xy
        x = points[0]
        y = points[1]
        plt.plot(x,y)
    
    # also plot original perimeter
    plt.plot(perimeterPointsPlot[:,0], perimeterPointsPlot[:,1], linestyle='dashdot',color='green', linewidth=1.0)
    
    
    return
#%% Find order of points acording to x and y axes within current grouping
def findOrderPoints(POINTS, x_vect, y_vect, refPoint):
       
    nPoints = POINTS.shape[0]
    dotResultsX = np.zeros(nPoints)
    dotResultsY = np.zeros(nPoints)
    pointsArray = np.zeros((nPoints,2))
    for ii in range(nPoints):
        currPoint = POINTS[ii,:]
        # pointsArray[ii,] = currPoint
        currVect = currPoint - refPoint
        dotResultsX[ii] = np.dot(x_vect, currVect)
        dotResultsY[ii] = np.dot(y_vect, currVect)
    posResultsX = np.argsort(dotResultsX)
    posResultsY = np.argsort(dotResultsY)
     
    score = np.zeros(nPoints)
    for ii in range(nPoints):
        res = np.where(posResultsX == ii)
        score[ii] = score[ii] + res[0][0]
        res = np.where(posResultsY == ii)
        score[ii] = score[ii] + res[0][0]
        
    pos = np.where(score == min(score))
    pos = pos[0][0]
    
    # Square = False
    # if posResultsX[0] != pos or posResultsX[1] != pos:
    #     Square = True
        
    return pos
# #%%
# def setPointZero(pointsArray, x_vect, y_vect,groupCenter):
#     nPoints = pointsArray.shape[0]
#     # vectors = np.zeros((nPoints,2))
    
#     print(pointsArray)
#     pos = findOrderPoints(pointsArray, x_vect, y_vect, groupCenter)
#     auxVects = pointsArray - groupCenter
       
#     cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
#     cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
    
#     cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
#     cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
    
#     if (cond1 or cond2) and (cond3 or cond4):
#         nShifts = pos % nPoints     
#         pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#     else:    
#         newX_vect = y_vect
#         newY_vect = -x_vect
#         pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
#         print(pos)
#         y_vect = newY_vect
#         x_vect = newX_vect
#         nShifts = pos % nPoints     
#         pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#         pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
 
#         auxVects = pointsArray - groupCenter    
#         cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
#         cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
        
#         cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
#         cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
        
#         if (cond1 or cond2) and (cond3 or cond4):
#             nShifts = pos % nPoints     
#             pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#         else:    
#             newX_vect = y_vect
#             newY_vect = -x_vect
#             pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
#             print(pos)
#             y_vect = newY_vect
#             x_vect = newX_vect
#             nShifts = pos % nPoints     
#             pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#             pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
                
#             auxVects = pointsArray - groupCenter
#             cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
#             cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
            
#             cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
#             cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
            
#             if (cond1 or cond2) and (cond3 or cond4):
#                 nShifts = pos % nPoints     
#                 pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#             else:    
#                 newX_vect = y_vect
#                 newY_vect = -x_vect
#                 pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
#                 print(pos)
#                 y_vect = newY_vect
#                 x_vect = newX_vect
#                 nShifts = pos % nPoints     
#                 pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#                 pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        
#                 auxVects = pointsArray - groupCenter    
#                 cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
#                 cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
                
#                 cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
#                 cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
                
#                 if (cond1 or cond2) and (cond3 or cond4):
#                     nShifts = pos % nPoints     
#                     pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
#                 else:    
#                     newX_vect = y_vect
#                     newY_vect = -x_vect
#                     pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
#                     print(pos)
#                     y_vect = newY_vect
#                     x_vect = newX_vect
#                     nShifts = pos % nPoints     
#                     pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))  
#                     pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        
        
      
#     # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
#     # AXIS_1.set_aspect('equal', 'box') 
#     # plt.title('set point 0')
#     # x = pointsArray[:,0]
#     # x = np.append(x,x[0])
#     # y = pointsArray[:,1]
#     # y = np.append(y,y[0])
#     # plt.plot(x,y)
#     # for jj in range(nPoints):  # do not enumerate last point = last
#     #     xx = pointsArray[jj,0]
#     #     yy = pointsArray[jj,1]
#     #     plt.text(xx,yy, str(jj), fontsize=8, color='red') 
        
#     # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
#     # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
#     # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
#     # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
#     # plt.plot(xaxis_x,xaxis_y, 'blue')
#     # plt.plot(yaxis_x,yaxis_y, 'red')
#     # plt.grid(which='major', axis='both')
    
    
#     return pointsArray, x_vect, y_vect
#%%
def setPointZero(pointsArray, x_vect, y_vect,groupCenter):
    nPoints = pointsArray.shape[0]
    # vectors = np.zeros((nPoints,2))
    
    # print(pointsArray)
    pos = findOrderPoints(pointsArray, x_vect, y_vect, groupCenter)
    # auxVects = pointsArray - groupCenter
       
    # cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
    # cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
    
    # cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
    # cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
     
        
    Indx0 = (pos-1) % nPoints
    Indx = pos % nPoints
    Indx1 = (pos+1) % nPoints
    Indx2 = (pos+2) % nPoints
    v0 = pointsArray[Indx0] - pointsArray[Indx]
    v1 = pointsArray[Indx1] - pointsArray[Indx]
    v2 = pointsArray[Indx2] - pointsArray[Indx1]
    
    condAng = np.cross(v1,v2) > 0 and np.cross(v1,v0) > 0
    Attempt = False
    while Attempt == False:    
        if condAng:
            Attempt = True
        else:
            pos = pos + 1
            Indx0 = (pos-1) % nPoints
            Indx = pos % nPoints
            Indx1 = (pos+1) % nPoints
            Indx2 = (pos+2) % nPoints
            v0 = pointsArray[Indx0] - pointsArray[Indx]
            v1 = pointsArray[Indx1] - pointsArray[Indx]
            v2 = pointsArray[Indx2] - pointsArray[Indx1]
            condAng = np.cross(v1,v2) > 0 and np.cross(v1,v0) > 0 

            
        
    nShifts = pos % nPoints     
    pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))    
        
        
        # nShifts = pos % nPoints     
        # pointsArray = np.roll(pointsArray,(-1,-1))
    
        
        # newX_vect = y_vect
        # newY_vect = -x_vect
        # pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        # print(pos)
        # y_vect = newY_vect
        # x_vect = newX_vect
        # nShifts = pos % nPoints     
        # pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        # pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
 
        # auxVects = pointsArray - groupCenter    
        # cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
        # cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
        
        # cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
        # cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
        
        # if (cond1 or cond2) and (cond3 or cond4):
        #     nShifts = pos % nPoints     
        #     pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        # else:    
        #     newX_vect = y_vect
        #     newY_vect = -x_vect
        #     pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        #     print(pos)
        #     y_vect = newY_vect
        #     x_vect = newX_vect
        #     nShifts = pos % nPoints     
        #     pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        #     pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
                
        #     auxVects = pointsArray - groupCenter
        #     cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
        #     cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
            
        #     cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
        #     cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
            
        #     if (cond1 or cond2) and (cond3 or cond4):
        #         nShifts = pos % nPoints     
        #         pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        #     else:    
        #         newX_vect = y_vect
        #         newY_vect = -x_vect
        #         pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        #         print(pos)
        #         y_vect = newY_vect
        #         x_vect = newX_vect
        #         nShifts = pos % nPoints     
        #         pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        #         pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        
        #         auxVects = pointsArray - groupCenter    
        #         cond1 = np.argsort(np.dot(auxVects, x_vect))[0] == pos
        #         cond2 = np.argsort(np.dot(auxVects, x_vect))[1] == pos
                
        #         cond3 = np.argsort(np.dot(auxVects, y_vect))[0] == pos
        #         cond4 = np.argsort(np.dot(auxVects, y_vect))[1] == pos
                
        #         if (cond1 or cond2) and (cond3 or cond4):
        #             nShifts = pos % nPoints     
        #             pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))
        #         else:    
        #             newX_vect = y_vect
        #             newY_vect = -x_vect
        #             pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        #             print(pos)
        #             y_vect = newY_vect
        #             x_vect = newX_vect
        #             nShifts = pos % nPoints     
        #             pointsArray = np.roll(pointsArray,(-nShifts,-nShifts))  
        #             pos = findOrderPoints(pointsArray, newX_vect, newY_vect, groupCenter)
        
        
      
    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('set point 0')
    # x = pointsArray[:,0]
    # x = np.append(x,x[0])
    # y = pointsArray[:,1]
    # y = np.append(y,y[0])
    # plt.plot(x,y)
    # for jj in range(nPoints):  # do not enumerate last point = last
    #     xx = pointsArray[jj,0]
    #     yy = pointsArray[jj,1]
    #     plt.text(xx,yy, str(jj), fontsize=8, color='red') 
        
    # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
    # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
    # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
    # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
    # plt.plot(xaxis_x,xaxis_y, 'blue')
    # plt.plot(yaxis_x,yaxis_y, 'red')
    # plt.grid(which='major', axis='both')
    
    
    return pointsArray, x_vect, y_vect
#%% from SQREDperiemeter extract final Rectangles
def GroupingExtractRects(x_vect, y_vect, SQREDperiemeter, groupCenterOverall, currGrouping):
    
    SQREDperiemeter.exterior.coords.xy
    points = SQREDperiemeter.exterior.coords.xy
    x = points[0]
    y = points[1]
    # plt.plot(x,y)
    nPoints = len(x)
    
    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('Grouping no.' + str(currGrouping))
    
    # for jj in range(nPoints-1):  # do not enumerate last point = last
    #     xx = points[0][jj]
    #     yy = points[1][jj]
    #     plt.text(xx,yy, str(jj), fontsize=8, color='red') 
        
    # xaxis_x = np.array([groupCenterOverall[0][0], groupCenterOverall[0][0] + \
    #                     10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
    # xaxis_y = np.array([groupCenterOverall[0][1], groupCenterOverall[0][1] + \
    #                     10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
    # yaxis_x = np.array([groupCenterOverall[0][0], groupCenterOverall[0][0] + \
    #                     10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
    # yaxis_y = np.array([groupCenterOverall[0][1], groupCenterOverall[0][1] + \
    #                     10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
    # plt.plot(xaxis_x,xaxis_y, 'blue')
    # plt.plot(yaxis_x,yaxis_y, 'red')
    # plt.grid(which='major', axis='both')
    
    nArrayPoints = nPoints-1
    pointsArray = np.zeros((nArrayPoints,2))  # going to work without the first point repeated
    for ii in range(nArrayPoints):
        pointsArray[ii,:] = np.array([points[0][ii],points[1][ii]])
    
    groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
    # pointsArray, x_vect, y_vect = firstPt2 = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
    pointsArray, x_vect, y_vect = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
    firstPt = 0

    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('Grouping no.' + str(currGrouping))
    # x = pointsArray[:,0]
    # x = np.append(x,x[0])
    # y = pointsArray[:,1]
    # y = np.append(y,y[0])
    # plt.plot(x,y)
    # for jj in range(nArrayPoints):  # do not enumerate last point = last
    #     xx = pointsArray[jj,0]
    #     yy = pointsArray[jj,1]
    #     plt.text(xx,yy, str(jj), fontsize=8, color='red') 
    
    
    # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
    # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
    # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
    # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
    # plt.plot(xaxis_x,xaxis_y, 'blue')
    # plt.plot(yaxis_x,yaxis_y, 'red')
    # plt.grid(which='major', axis='both')
      
    rect0 = np.zeros((4,2))    # RECTANGLE to be sustracted   
    rect0[0,:] = pointsArray[firstPt-1,:]
    rect0[1,:] = pointsArray[firstPt,:]
    rect0[2,:] = pointsArray[firstPt+1,:]
    rect0[3,:] = pointsArray[firstPt+2,:]
          
    len0 = np.linalg.norm(rect0[1,:]-rect0[0,:]) 
    len1 = np.linalg.norm(rect0[2,:]-rect0[1,:]) 
    len2 = np.linalg.norm(rect0[3,:]-rect0[2,:]) 
    
    ListofRects = []
    TheEnd = False    
    while TheEnd == False  :
        
        if nArrayPoints-4 > 0:       # check and exit while loop ----->
            
            if abs(len0 - len2) < 1e-6:  # check if equal side lengths
        
                ListofRects.append(rect0)
                
                # plt.plot(rect0[:,0],rect0[:,1], linestyle='solid',color='black', linewidth=2.0) 
                       
                newPointsArray = np.zeros((nArrayPoints-4,2))
                counter = 0
                for ii in range(nArrayPoints):
                    # print(ii)
                    cond1x = pointsArray[ii,0] != rect0[0,0]
                    cond2x = pointsArray[ii,0] != rect0[1,0]
                    cond3x = pointsArray[ii,0] != rect0[2,0]
                    cond4x = pointsArray[ii,0] != rect0[3,0]
                    cond1y = pointsArray[ii,1] != rect0[0,1]
                    cond2y = pointsArray[ii,1] != rect0[1,1]
                    cond3y = pointsArray[ii,1] != rect0[2,1]
                    cond4y = pointsArray[ii,1] != rect0[3,1]
                    condX = cond1x and cond2x and cond3x and cond4x
                    condY = cond1y and cond2y and cond3y and cond4y
                    if  condX and condY:
                        newPointsArray[counter,:] = pointsArray[ii,:]
                        counter = counter + 1
                        
                # plt.plot(newPointsArray[:,0],newPointsArray[:,1], linestyle='solid',color='red', linewidth=2.0)           
                
                nNewPointsArray = newPointsArray.shape[0]
               
                newPointsArray = np.roll(newPointsArray,(1,1))    # roll around one position
                
                pointsArray = newPointsArray
                
                # groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                # # firstPt = findOrderPoints(pointsArray, x_vect, y_vect,groupCenter)
                # firstPt = 0
                
                groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                # pointsArray, x_vect, y_vect = firstPt2 = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                pointsArray, x_vect, y_vect = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                firstPt = 0
                
                rect0 = np.zeros((4,2))    # RECTANGLE to be sustracted
                rect0[0,:] = pointsArray[firstPt-1,:]
                rect0[1,:] = pointsArray[firstPt,:]
                rect0[2,:] = pointsArray[firstPt+1,:]
                rect0[3,:] = pointsArray[firstPt+2,:]
                      
                len0 = np.linalg.norm(rect0[1,:]-rect0[0,:]) 
                len1 = np.linalg.norm(rect0[2,:]-rect0[1,:]) 
                len2 = np.linalg.norm(rect0[3,:]-rect0[2,:]) 
                
                nArrayPoints = pointsArray.shape[0]
                
                # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
                # AXIS_1.set_aspect('equal', 'box') 
                # plt.title('Grouping no.' + str(currGrouping))
                # x = newPointsArray[:,0]
                # x = np.append(x,x[0])
                # y = newPointsArray[:,1]
                # y = np.append(y,y[0])
                # plt.plot(x,y)
                # for jj in range(nNewPointsArray):  # do not enumerate last point = last
                #     xx = newPointsArray[jj,0]
                #     yy = newPointsArray[jj,1]
                #     plt.text(xx,yy, str(jj), fontsize=8, color='red') 
                    
                # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
                # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
                
                # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
                # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
                
                # plt.plot(xaxis_x,xaxis_y, 'blue')
                # plt.plot(yaxis_x,yaxis_y, 'red')
                # plt.grid(which='major', axis='both')
                # a = 1    
            else:            # we don't have equal side lengths
                if len0 < len2:
                    
                    testVect =  (rect0[3,:] - rect0[2,:])/np.linalg.norm(rect0[3,:] - rect0[2,:])
                    
                    if abs(np.dot(testVect,x_vect)) > abs(np.dot(testVect,y_vect)):
                        rect0[3,:] = rect0[2,:] + len0*x_vect
                    else:
                        rect0[3,:] = rect0[2,:] + len0*y_vect
                        
                    ListofRects.append(rect0)    
                        
                    newPointsArray = np.zeros((nArrayPoints-2,2))
                    counter = 0
                    for ii in range(nArrayPoints):
                        cond1x = pointsArray[ii,0] != rect0[0,0] 
                        cond2x = pointsArray[ii,0] != rect0[1,0]
                        cond3x = pointsArray[ii,0] != rect0[2,0]
                        # cond4x = pointsArray[ii,0] != rect0[3,0]  # do not remove point0 
                        cond1y = pointsArray[ii,1] != rect0[0,1] 
                        cond2y = pointsArray[ii,1] != rect0[1,1]
                        cond3y = pointsArray[ii,1] != rect0[2,1]
                        # cond4y = pointsArray[ii,1] != rect0[3,1]   # do not remove point0
                        condX = cond1x and cond2x and cond3x    # <---------------- cond4x removed
                        condY = cond1y and cond2y and cond3y    # <---------------- cond4y removed    
                        
                        if counter == 0:
                            newPointsArray[counter,:] = rect0[3,:]
                            counter = counter + 1
                        if  condX and condY:
                            newPointsArray[counter,:] = pointsArray[ii,:]
                            counter = counter + 1
                            
                    # newPointsArray = np.vstack((rect0[3,:], newPointsArray))   # add point where we cut
                    
                    nNewPointsArray = newPointsArray.shape[0]
                    
                    newPointsArray = np.roll(newPointsArray,(1,1))    # roll around one position
                    
                    pointsArray = newPointsArray
                    
                    # groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                    # # firstPt = findOrderPoints(pointsArray, x_vect, y_vect,groupCenter)
                    # firstPt = 0
                    
                    groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                    # pointsArray, x_vect, y_vect = firstPt2 = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                    pointsArray, x_vect, y_vect = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                    firstPt = 0
                    
                    rect0 = np.zeros((4,2))    # RECTANGLE to be sustracted
                    rect0[0,:] = pointsArray[firstPt-1,:]
                    rect0[1,:] = pointsArray[firstPt,:]
                    rect0[2,:] = pointsArray[firstPt+1,:]
                    rect0[3,:] = pointsArray[firstPt+2,:]
                          
                    len0 = np.linalg.norm(rect0[1,:]-rect0[0,:]) 
                    len1 = np.linalg.norm(rect0[2,:]-rect0[1,:]) 
                    len2 = np.linalg.norm(rect0[3,:]-rect0[2,:]) 
                    
                    nArrayPoints = pointsArray.shape[0]
                    
                    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
                    # AXIS_1.set_aspect('equal', 'box') 
                    # plt.title('Grouping no.' + str(currGrouping))
                    # x = newPointsArray[:,0]
                    # x = np.append(x,x[0])
                    # y = newPointsArray[:,1]
                    # y = np.append(y,y[0])
                    # plt.plot(x,y)
                    # for jj in range(nNewPointsArray):  # do not enumerate last point = last
                    #     xx = newPointsArray[jj,0]
                    #     yy = newPointsArray[jj,1]
                    #     plt.text(xx,yy, str(jj), fontsize=8, color='red')    
                    
                    # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
                    # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
                    
                    # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
                    # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
                    
                    # plt.plot(xaxis_x,xaxis_y, 'blue')
                    # plt.plot(yaxis_x,yaxis_y, 'red')
                    # plt.grid(which='major', axis='both')
                    
                    # a = 1
                else:   # len0 > len2  -----> keep rect0's side 2 of len2 
                    
                    testVect =  (rect0[0,:] - rect0[1,:])/np.linalg.norm(rect0[0,:] - rect0[1,:])
                    
                    if abs(np.dot(testVect,x_vect)) > abs(np.dot(testVect,y_vect)):
                        rect0[0,:] = rect0[1,:] + len2*x_vect
                    else:
                        rect0[0,:] = rect0[1,:] + len2*y_vect
                        
                    ListofRects.append(rect0)    
                        
                    newPointsArray = np.zeros((nArrayPoints-3,2))
                    counter = 0
                    for ii in range(nArrayPoints):
                        # print('------------')
                        # print(ii)
                        # cond1x = pointsArray[ii,0] != rect0[0,0]  # do not remove point0 
                        cond2x = pointsArray[ii,0] != rect0[1,0]
                        cond3x = pointsArray[ii,0] != rect0[2,0]
                        cond4x = pointsArray[ii,0] != rect0[3,0]
                        # cond1y = pointsArray[ii,1] != rect0[0,1]   # do not remove point0 
                        cond2y = pointsArray[ii,1] != rect0[1,1]
                        cond3y = pointsArray[ii,1] != rect0[2,1]
                        cond4y = pointsArray[ii,1] != rect0[3,1]
                        condX = cond2x and cond3x and cond4x     # <---------------- cond1x removed
                        condY = cond2y and cond3y and cond4y     # <---------------- cond1y removed
                        if  condX and condY:
                            newPointsArray[counter,:] = pointsArray[ii,:]
                            counter = counter + 1
                            # print(counter)
                    newPointsArray = np.vstack((newPointsArray, rect0[0,:]))   # add point where we cut
                    # newPointsArray = np.vstack((rect0[0,:], newPointsArray))   # add point where we cut
                    
                    nNewPointsArray = newPointsArray.shape[0]
                    
                    newPointsArray = np.roll(newPointsArray,(1,1))    # roll around one position
                    
                    pointsArray = newPointsArray
                    
                    # groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                    # # firstPt = findOrderPoints(pointsArray, x_vect, y_vect,groupCenter)
                    # firstPt = 0
                    
                    groupCenter = np.sum(pointsArray, axis=0)/pointsArray.shape[0]    
                    # pointsArray, x_vect, y_vect = firstPt2 = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                    pointsArray, x_vect, y_vect = setPointZero(pointsArray, x_vect, y_vect,groupCenter)
                    firstPt = 0
                    
                    rect0 = np.zeros((4,2))    # RECTANGLE to be sustracted               
                    rect0[0,:] = pointsArray[firstPt-1,:]
                    rect0[1,:] = pointsArray[firstPt,:]
                    rect0[2,:] = pointsArray[firstPt+1,:]
                    rect0[3,:] = pointsArray[firstPt+2,:]
                          
                    len0 = np.linalg.norm(rect0[1,:]-rect0[0,:]) 
                    len1 = np.linalg.norm(rect0[2,:]-rect0[1,:]) 
                    len2 = np.linalg.norm(rect0[3,:]-rect0[2,:]) 
                    
                    nArrayPoints = pointsArray.shape[0] 
                                        
                    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
                    # AXIS_1.set_aspect('equal', 'box') 
                    # plt.title('Grouping no.' + str(currGrouping))
                    # x = pointsArray[:,0]
                    # x = np.append(x,x[0])
                    # y = pointsArray[:,1]
                    # y = np.append(y,y[0])
                    # plt.plot(x,y)
                    # for jj in range(nNewPointsArray):  # do not enumerate last point = last
                    #     xx = pointsArray[jj,0]
                    #     yy = pointsArray[jj,1]
                    #     plt.text(xx,yy, str(jj), fontsize=8, color='red')    
                    
                    # xaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
                    # xaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
                    
                    # yaxis_x = np.array([groupCenter[0], groupCenter[0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
                    # yaxis_y = np.array([groupCenter[1], groupCenter[1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
                    
                    # plt.plot(xaxis_x,xaxis_y, 'blue')
                    # plt.plot(yaxis_x,yaxis_y, 'red')
                    # plt.grid(which='major', axis='both')
                    
                    # a = 1
        else:
            ListofRects.append(rect0)
            TheEnd = True  
    
    
    # # plot now the resulting rectangles -------
    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('Grouping no.' + str(currGrouping))    
    # for ii in range(len(ListofRects)):
    #     currRect = ListofRects[ii]
    #     x = currRect[:,0]
    #     x = np.append(x,x[0])
    #     y = currRect[:,1]
    #     y = np.append(y,y[0])
    #     plt.plot(x,y)
    # plt.grid(which='major', axis='both')
    
    # xaxis_x = np.array([groupCenterOverall[0][0], groupCenterOverall[0][0] + \
    #                     10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
    # xaxis_y = np.array([groupCenterOverall[0][1], groupCenterOverall[0][1] + \
    #                     10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
    # yaxis_x = np.array([groupCenterOverall[0][0], groupCenterOverall[0][0] + \
    #                     10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
    # yaxis_y = np.array([groupCenterOverall[0][1], groupCenterOverall[0][1] + \
    #                     10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
    # plt.plot(xaxis_x,xaxis_y, 'blue')
    # plt.plot(yaxis_x,yaxis_y, 'red')
    
    # STREIPPED_BLOCK = Polygon(newPointPerimeter)
    
    ListofRectsPolygs = []
    for ii in range(len(ListofRects)):
        ListofRectsPolygs.append(Polygon(ListofRects[ii]))
    
    
    return ListofRects, ListofRectsPolygs
#%% create distance between points matrix
def distMatrix(POLYLIST):
    matrixDist = np.zeros([len(POLYLIST)*4,len(POLYLIST)*4])
    for ii in range(len(POLYLIST)):
        poly1 = POLYLIST[ii]
        for jj in range(len(POLYLIST)):
            if ii == jj:
                continue
            else:
                poly2 = POLYLIST[jj]
                for kk in range(4):
                    for ll in range(4):
                        x1 = poly1.exterior.coords.xy[0][kk]
                        y1 = poly1.exterior.coords.xy[1][kk]
                        x2 = poly2.exterior.coords.xy[0][kk]
                        y2 = poly2.exterior.coords.xy[1][kk]
                        dis = np.sqrt((x1-x2)**2 +(y1-y2)**2)
                        matrixDist[ii*4+kk,jj*4+ll] = dis
                        a = 1
    return matrixDist   
#%% find point in list of polygons
def findInPolyList(POLYLIST, point):
    for ii in range(len(POLYLIST)):
        currPoly = POLYLIST[ii]
        for jj in range(4):
            if currPoly.exterior.coords.xy[0][jj] == point[0] and \
                currPoly.exterior.coords.xy[1][jj] == point[1]:
                # print('there a match poly no. '+str(ii))
                # print('             point no. '+str(jj))
                whatPoly = ii
                whatPoint = jj
    return whatPoly, whatPoint
#%% New approach to step 3. Trim rotated, overlapping rectangles
def groupingStep3trim(x_vect, y_vect, groupCenter, LOCALROTATED, matrixDist, perimeterPointsPlot, groupings, currGrouping):
    
    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('Grouping no.' + str(currGrouping))
    
    polygonsINGroup = []
    LensInGroup = []
    RotsInGroup = []
    
    aauxPerimeterX = []
    aauxPerimeterY = []
    
    pointCounter = 0
    for ii in range(len(LOCALROTATED)):        
                         
        currPolygon = LOCALROTATED[ii]
        points = currPolygon.exterior.coords.xy
        x = points[0]
        y = points[1]
        # plt.plot(x,y)
        
        for jj in range(4):  # enumerate rect. corners
            xx = points[0][jj]
            yy = points[1][jj]
            # plt.text(xx,yy, str(jj), fontsize=8, color='red') #plt.text(xx,yy, str(pointCounter), fontsize=8, color='red')
            pointCounter += 1
        
        # plt.plot(currPolygon.centroid.xy[0][0], currPolygon.centroid.xy[1][0],'o')
        # plt.text(currPolygon.centroid.xy[0][0], currPolygon.centroid.xy[1][0],str(ii))
        
        for kk in range(4):
            aauxPerimeterX.append(x[kk])
            aauxPerimeterY.append(y[kk])
            
    # aauxPerimeterX.append(aauxPerimeterX[0])
    # aauxPerimeterY.append(aauxPerimeterY[0])
    
    LenPeri = len(aauxPerimeterX)
    aauxPerimeter2 = np.zeros((LenPeri,2))
    aauxPerimeter2[:,0] = aauxPerimeterX
    aauxPerimeter2[:,1] = aauxPerimeterY
    
    Perimeter = ConvexHull(aauxPerimeter2)
    pointter = Perimeter.vertices                # pointer to segment point
    newPerimeter = aauxPerimeter2[pointter,:]
    newPerimeter = np.vstack((newPerimeter, newPerimeter[0,:]))
    
    # plt.plot(newPerimeter[:,0], newPerimeter[:,1], linestyle='solid',color='black', linewidth=1.0)
    # also plot original perimeter
    # plt.plot(perimeterPointsPlot[:,0], perimeterPointsPlot[:,1], linestyle='dashdot',color='green', linewidth=1.0)
    
    # for jj in range(newPerimeter.shape[0]-1):
    #     plt.text(newPerimeter[jj,0],newPerimeter[jj,1], str(jj), fontsize=12, color='blue')

    MINRECTbbox, Rotbbox, Lenbbox = fitRotRectangle(newPerimeter, newPerimeter.shape[0])  # I need this for plotting main axes
    BoundingBoxPoints = MINRECTbbox.exterior.coords.xy
    
    # plt.plot(BoundingBoxPoints[0][:], BoundingBoxPoints[1][:],linestyle='dotted',color='yellow', linewidth=1.0)
    # plt.plot(groupCenter[0][0], groupCenter[0][1],'*', markersize=20)
    
    # xaxis_x = np.array([groupCenter[0][0], groupCenter[0][0] + 10*np.cos(np.arctan2(x_vect[1],x_vect[0]))])
    # xaxis_y = np.array([groupCenter[0][1], groupCenter[0][1] + 10*np.sin(np.arctan2(x_vect[1],x_vect[0]))])
    
    # yaxis_x = np.array([groupCenter[0][0], groupCenter[0][0] + 10*np.cos(np.arctan2(y_vect[1],y_vect[0]))])
    # yaxis_y = np.array([groupCenter[0][1], groupCenter[0][1] + 10*np.sin(np.arctan2(y_vect[1],y_vect[0]))])
    
    # plt.plot(xaxis_x,xaxis_y, 'blue')
    # plt.plot(yaxis_x,yaxis_y, 'red')
    
    # plt.grid(which='major', axis='both')
   
    periVect = np.zeros([newPerimeter.shape[0]-1, 2])
    for jj in range(newPerimeter.shape[0]-1):   
        periVect[jj,] = (newPerimeter[jj+1,] - newPerimeter[jj,])/ \
            np.linalg.norm(newPerimeter[jj+1,] - newPerimeter[jj,])
        
    perimeter2Pointer = []   #[pointer to perimeter segment, pointer to segment point]
    for jj in range(newPerimeter.shape[0]-1):      # this loop is along perimeter points jj
        test_x = abs(np.dot(x_vect, periVect[jj,]))
        if abs(test_x) < 1e-6:
            test_x = 0
        if abs(test_x -1) < 1e-6:
            test_x = 1
        test_y = abs(np.dot(y_vect, periVect[jj,]))
        if abs(test_y) < 1e-6:
            test_y = 0
        if abs(test_y -1) < 1e-6:
            test_y = 1
        
        if (abs(test_x - 1.0) < 1e-6) or (abs(test_x) < 1e-6):  # side not along perimeter segment
            perimeter2Pointer.append([jj, pointter[jj]])  #[pointer to perimeter segment, pointer to segment point]
                       
    # run loop a second time to deaL with non attached perimeter segments find new attached segments       
    newPerimeterIndex = np.arange(newPerimeter.shape[0])  
    newPointPerimeter = np.array(np.reshape(newPerimeter[0,],(1,2)))  # include first point of Perimeter
    perimeter2Pointer = np.array(perimeter2Pointer)              # from list to npArray
    
    
    perimeter2Pointer = np.vstack((perimeter2Pointer,perimeter2Pointer[0,])) ########
    perimeter2Pointer[-1,0] = newPerimeter.shape[0] - 1 #############################
     
    perimeter2PointerFirstColumn = perimeter2Pointer[:,0]
    
    pointToTest = newPerimeter[0,]    
    jj = 0
    while jj < newPerimeter.shape[0]-1:            # this loop is along perimeter points jj
    
        whatPoly, whatPoint = findInPolyList(LOCALROTATED, pointToTest)
        
        if jj + 1 < newPerimeter.shape[0]:
            pointToTest2 = newPerimeter[jj+1,]    
            whatPoly2, whatPoint2 = findInPolyList(LOCALROTATED, pointToTest2)
    
            if whatPoly == whatPoly2:
                newPointPerimeter = np.vstack((newPointPerimeter, pointToTest2))
                jj = jj + 1
                pointToTest = pointToTest2
            
                if jj + 1 < newPerimeter.shape[0]:
                    pointToTest3 = newPerimeter[jj+1,]    
                    whatPoly3, whatPoint3 = findInPolyList(LOCALROTATED, pointToTest3)
            
                    if whatPoly2 == whatPoly3:
                        newPointPerimeter = np.vstack((newPointPerimeter, pointToTest3))
                        jj = jj + 1
                        pointToTest = pointToTest3
                        
        if jj < newPerimeter.shape[0]-1: 
            
            Finish = False  # Condition to start while loop
            while Finish == False:
                
                whatPoly, whatPoint = findInPolyList(LOCALROTATED, pointToTest)
                # whatPoint is computed locally to that polygon
                # whatPoly is computed locally in the current grouping
            
                polyUnderTest0 = LOCALROTATED[whatPoly]
                x2 = polyUnderTest0.exterior.coords.xy[0][whatPoint+1] # initial point and next one on same rectangle
                x1 = polyUnderTest0.exterior.coords.xy[0][whatPoint] 
                y2 = polyUnderTest0.exterior.coords.xy[1][whatPoint+1] 
                y1 = polyUnderTest0.exterior.coords.xy[1][whatPoint] 
                P1 = np.array([x1, y1])
                P2 = np.array([x2, y2])
                vP1P2 = (P2 - P1)/np.linalg.norm(P2 - P1)
                
                # auxq = np.where(perimeter2PointerFirstColumn > jj)    # MALLLLLLLLLLLL
                # nextSQRperiIndex = perimeter2PointerFirstColumn[auxq][0]
                # nextSQRperiPoint = newPerimeter[nextSQRperiIndex,]
                nextSQRperiPoint = newPerimeter[jj+1,]
                
                maxPossibleDist = np.dot(vP1P2, nextSQRperiPoint - P1)
                maxPossibleDist = maxPossibleDist*1e6
                maxPossibleDist = round(maxPossibleDist)
                maxPossibleDist = maxPossibleDist*1e-6
                
                # angT1 = np.inf
                # angT2 = np.inf
                P = np.array([np.inf,  np.inf]) # dummy incersection point
                decision2 = [] #[-1, -1, dist, angT1, angT2]
                decision2Aux = []
                for tt in range(len(LOCALROTATED)):
                    if tt != whatPoly:
                        polyUnderTest = LOCALROTATED[tt]
                        for rr in range(4):
                            x2 = polyUnderTest.exterior.coords.xy[0][rr+1] 
                            x1 = polyUnderTest.exterior.coords.xy[0][rr] 
                            y2 = polyUnderTest.exterior.coords.xy[1][rr+1] 
                            y1 = polyUnderTest.exterior.coords.xy[1][rr] 
                            T1 = np.array([x1, y1])
                            T2 = np.array([x2, y2])
                            vT1T2 = (T2-T1)/np.linalg.norm(T2-T1)  
                            dotProduct = np.dot(vP1P2, vT1T2)
                            P = line_intersect(P1, P2, T1, T2)
                            dist = np.linalg.norm(P - P1)
                            
                            signedDist = np.dot(P - P1, vP1P2)
                            
                            P1P = P - P1
                            PT2 = T2 - P
                            PT1 = T1 - P
                            
                            xx = dist
                            yy1 = np.dot(vT1T2, PT1)
                            yy2 = np.dot(vT1T2, PT2)
                            
                            if not(np.isnan(signedDist)) and signedDist > 0.0  and signedDist < 1e9:
                                # print(dist)
                                dist = dist*1e6                     
                                dist = round(dist)
                                dist = dist*1e-6
                                if dist <= maxPossibleDist:
                                    # auxxDecision2 = [tt, rr, signedDist, yy1, yy2, P]
                                    auxxDecision2 = [tt, rr, signedDist, yy1, yy2]
                                    decision2.append(auxxDecision2) 
                                    decision2Aux.append(P)
  
                decision2 = np.array(decision2)  ###### rephase this statement
                decision2Aux = np.array(decision2Aux)
                # print(decision2)
                # print(decision2Aux)
                
                distCol = decision2[:,2]
                sortedDist = np.argsort(distCol)
                Found = False
                Failed = False
                nn = 0
                while Found == False:        
                    pointt = int(sortedDist[nn])
                    side1 = decision2[pointt,3]
                    signside1 = np.sign(side1)
                    side2 = decision2[pointt,4]
                    signside2 = np.sign(side2)
                    if signside1 == signside2:
                        nn = nn + 1
                        if nn == distCol.shape[0]:
                            Failed = True      # means we need to softened criteria
                            break
                    else:
                        Found = True # Output is nn points to pointer list sortedDist
                        # alternative output is Failed flag         
 
                # should the selection fail because the interesection point P is not between T1 and T2
                # settle for the rectangle yielding the smallest angle
                if Failed == True: 
                    minSide = 1e9
                    pointt = int(-1)
                    for kk in range(decision2.shape[0]):
                        side1 = decision2[kk,3]
                        side2 = decision2[kk,4]
                        currMinSide = min(side1, side2)
                        
                        
                        tt1 = int(decision2[kk,0])
                        polyUnderTest = LOCALROTATED[tt1]
                        rr1 = int(decision2[kk,1])
                        xx1 = polyUnderTest.exterior.coords.xy[0][rr1] 
                        yy1 = polyUnderTest.exterior.coords.xy[1][rr1] 
                        v2 = np.array([xx1,yy1]) - P1
                        v2 = v2/np.linalg.norm(v2)
                        
                        v1 = vP1P2
                        
                        testAngle = np.cross(v1,v2)
                        
                        if minSide > currMinSide and testAngle < 0:
                            pointt = kk
                            minSide = currMinSide
                
                # endPoint = decision2[pointt][5]        # this is intersection point 
                endPoint = decision2Aux[pointt,]        # this is intersection point 
                auxx = np.reshape(endPoint,(1,2))                       
                newPointPerimeter = np.vstack((newPointPerimeter, auxx))
                
                whatPoint2 = int(decision2[pointt][1])+1
                whatPoly2 = int(decision2[pointt][0])
                
                polyUnderTest0 = LOCALROTATED[whatPoly2]
                xx = polyUnderTest0.exterior.coords.xy[0][whatPoint2]
                yy = polyUnderTest0.exterior.coords.xy[1][whatPoint2]         
                endPoint = np.array([xx, yy])
                auxx = np.reshape(endPoint,(1,2)) 
                newPointPerimeter = np.vstack((newPointPerimeter, auxx))
                
                test3 = np.where(newPerimeter[:,0] == auxx[0,0])  # does new perimeterpoint belong to perimeter???
                test4 = np.where(newPerimeter[:,1] == auxx[0,1])
                              
                if np.any(test3) == True and np.any(test4) == True:
                    Finish = True
                    jj = jj + 1    
                              
                pointToTest = endPoint    # in any case, the oputput of the while loop is pointToTest      
                a = 1
                    
            
    # plt.plot(newPointPerimeter[:,0], newPointPerimeter[:,1], linestyle='solid',color='red', linewidth=3.0)    
    
    SQREDperiemeter = Polygon(newPointPerimeter)
    
    
    return SQREDperiemeter
#%% check if two polygons touch
# From: https://stackoverflow.com/questions/38481437/python-shapely-how-to-determine-if-two-polygons-cross-each-other-while-allowi
def myTouches(poly1, poly2):
    return poly1.intersects(poly2) and not poly1.crosses(poly2) and not poly1.contains(poly2)
#%%
def checkRectIntersection(Poly1,Poly2, vector):
    InstersFlag = False
    if Poly1.intersects(Poly2) == True:
        shiftLen = measureShiftLen(Poly1, Poly2, vector)
        print(shiftLen)
        InstersFlag = True
        a = 1
        
    return InstersFlag

#%%
# def checkRectIntersection(Poly1,Poly2):
#     InstersFlag = True
#     if Poly1.intersects(Poly2) == False:
#         InstersFlag = False
#     else:
#         Touches = myTouches(Poly1, Poly2)
#         if Touches == True:
#             InstersFlag = False
#     return InstersFlag
#%% shift Rectangle along vector
def shiftAlongVector(Poly, vect, shiftLen):
    
    shiftVector = shiftLen*vect
    
    auxx = np.zeros((5,2))
    
    points = Poly.exterior.coords.xy
    auxx[:,0] = points[0]
    auxx[:,1] = points[1]
    # auxx = np.reshape(points, (5,2))
    newPoints = auxx + shiftVector
    
    
    shiftedPoly = Polygon(newPoints)
    
    return shiftedPoly
#%% Measure shiftLength necessary
def measureShiftLen(Poly1, Poly2, vect):
    shiftLen = 0.0
    points_1 = np.zeros((4,2))
    auxx = np.reshape(np.array(Poly1.exterior.coords.xy[0][0:4]), (4,1))
    points_1[:,0] = auxx[:,0]
    auxx = np.reshape(np.array(Poly1.exterior.coords.xy[1][0:4]), (4,1))
    points_1[:,1] = auxx[:,0]
    
    currCentr_1 = np.array(Poly1.centroid.coords.xy)
    
    points_2 = np.zeros((4,2))
    auxx = np.reshape(np.array(Poly2.exterior.coords.xy[0][0:4]), (4,1))
    points_2[:,0] = auxx[:,0]
    auxx = np.reshape(np.array(Poly2.exterior.coords.xy[1][0:4]), (4,1))
    points_2[:,1] = auxx[:,0]
    
    currCentr_2 = np.array(Poly2.centroid.coords.xy)
    
    vectCentr1_Centr2 = currCentr_2 - currCentr_1
    vectCentr1_Centr2 = np.reshape(vectCentr1_Centr2,(1,2))
    distCentr1_Centr2 = np.dot(vect, vectCentr1_Centr2[0])
    
    
    sides_1 = np.zeros((4,2)) 
    sides_2 = np.zeros((4,2)) 
    for ii in range(4):
        if ii < 3:
            sides_1[ii,:] = points_1[ii+1,:] - points_1[ii,:]
            sides_2[ii,:] = points_2[ii+1,:] - points_2[ii,:]
        else:
            sides_1[ii,:] = points_1[0,:] - points_1[ii,:]
            sides_2[ii,:] = points_2[0,:] - points_2[ii,:]
    
    side_Poly1 = 0
    side_Poly2 = 0
    for ii in range(4):
        side_Poly1 = np.linalg.norm(np.dot(vect, sides_1[ii,:]))
        if side_Poly1 > 1e-6:
            break
      
    for ii in range(4):
        side_Poly2 = np.linalg.norm(np.dot(vect, sides_2[ii,:]))
        if side_Poly2 > 1e-6:
            break
         
    shiftLen = (side_Poly1 + side_Poly2)/2 - distCentr1_Centr2
    return shiftLen
#%% 3rd Grouping-oriented processing step. SHIFT ACCORDING TO REFERENCE DIRECTION x or y
def groupingStep3shift(vector, vectOrder, LOCALROTATED, perimeterPointsPlot, groupings, currGrouping):
    LOCALSHIFTEDvect = []
    
    currGrMembers = groupings[currGrouping]
    AUXX = np.array(LOCALROTATED)   # arranged as in currGrouping
    
    if len(LOCALROTATED)> 0:
        for ii in range(0,len(currGrMembers)-1):
            auxx = currGrMembers.index(vectOrder[ii])
            pushingRectagle = AUXX[auxx]
            for jj in range(ii+1,len(currGrMembers)):
                auxx2 = currGrMembers.index(vectOrder[jj])
                pushedRectagle = AUXX[auxx2]
                
                InstersFlag = checkRectIntersection(pushingRectagle,pushedRectagle, vector)
                
                if InstersFlag == True: 
                   
                    shiftLen = measureShiftLen(pushingRectagle, pushedRectagle, vector) # Measure the shift
                                        
                    shiftedPoly = shiftAlongVector(pushedRectagle, vector, shiftLen) # do the shift

                    AUXX[auxx2] = shiftedPoly
     
        LOCALSHIFTEDvect = AUXX.tolist()              
                    
        fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
        AXIS_1.set_aspect('equal', 'box') 
        plt.title('Grouping no.' + str(currGrouping))
    
        for ii in range(len(LOCALSHIFTEDvect)):
            currPolygon = LOCALSHIFTEDvect[ii]
            points = currPolygon.exterior.coords.xy
            x = points[0]
            y = points[1]
            plt.plot(x,y)
            plt.plot(currPolygon.centroid.coords.xy[0][0],currPolygon.centroid.coords.xy[1][0],'o')
            plt.text(currPolygon.centroid.coords.xy[0][0],currPolygon.centroid.coords.xy[1][0], str(currGrMembers[ii]))
        plt.plot(perimeterPointsPlot[:,0], perimeterPointsPlot[:,1], \
                 linestyle='solid',color='red', linewidth=2)             
    

                
    return LOCALSHIFTEDvect
#%% 2nd Grouping-oriented processing step. ROTATION
def groupingStep2rotation(groupCenter, perimeterPointsPlot, \
    ALLBBOXES, ROTATIONS, LensInGroup, RotsInGroup, groupings, currGrouping):
 
    if len(groupings[currGrouping]) > 1:
        pointt = LensInGroup.index(max(LensInGroup))  # max side length determines rotation of all polygons
        refRot = RotsInGroup[pointt]
        LOCALROTATED = []
        
        
        # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
        # AXIS_1.set_aspect('equal', 'box') 
        # plt.title('Grouping no.' + str(currGrouping))
        
        for ii in groupings[currGrouping]:
            # print(ii)
            currPolygon = ALLBBOXES[ii]
            rotAngleDeg = refRot - ROTATIONS[ii][0]
            
            new_points, rotPolygon = rotatePolygon(currPolygon, rotAngleDeg)
            
            LOCALROTATED.append(rotPolygon)      # LOCALROTATED.append(new_points)
            
        #     plt.plot(new_points[0],new_points[1])
        #     currCentr = np.array(rotPolygon.centroid.coords.xy)
        #     plt.plot(currCentr[0],currCentr[1], 'o')
        #     plt.text(currCentr[0],currCentr[1], str(ii))
        #     for jj in range(4):  # enumerate rect. corners
        #         xx = new_points[0][jj]
        #         yy = new_points[1][jj]
        #         plt.text(xx,yy, str(jj), fontsize=10, color='red')
                
            
        # plt.plot(perimeterPointsPlot[:,0], perimeterPointsPlot[:,1], \
        #          linestyle='dashed',color='green', linewidth=0.5)    
        # plt.plot(groupCenter[0][0], groupCenter[0][1],'*', markersize=20)    
    else:
        LOCALROTATED = []
        refRot = 0.0
 
    
    
    
    return LOCALROTATED, refRot
#%% 1st Grouping-oriented processing step. 
def groupingStep1(ALLBBOXES, ROTATIONS, LENGTHS, groupings, currGrouping):
    
    # fig, AXIS_1 = plt.subplots(1, figsize=(9, 9))
    # AXIS_1.set_aspect('equal', 'box') 
    # plt.title('Grouping no.' + str(currGrouping))
    
    
    polygonsINGroup = []
    LensInGroup = []
    RotsInGroup = []
    
    aauxCenter = np.zeros([1,2])
    aauxPerimeterX = []
    aauxPerimeterY = []
    for ii in range(len(ALLBBOXES)):        
             
       Location = findInlistofArrays(ii,groupings)
       if Location == currGrouping:
                  
           currPolygon = ALLBBOXES[ii]
           points = currPolygon.exterior.coords.xy
           x = points[0]
           y = points[1]
           # plt.plot(x,y)
           
           aauxCenter += [sum(x[0:4]), sum(y[0:4])]
           
           for kk in range(4):
               aauxPerimeterX.append(x[kk])
               aauxPerimeterY.append(y[kk])
           
           currCentr = np.array(currPolygon.centroid.coords.xy)
           # plt.plot(currCentr[0],currCentr[1], 'o')
           # plt.text(currCentr[0],currCentr[1], str(ii))
           
           xaxis_x = np.array([currCentr[0][0], currCentr[0][0] + LENGTHS[ii][0]*np.cos(ROTATIONS[ii][0]*np.pi/180)])
           xaxis_y = np.array([currCentr[1][0], currCentr[1][0] + LENGTHS[ii][0]*np.sin(ROTATIONS[ii][0]*np.pi/180)])
           
           yaxis_x = np.array([currCentr[0][0], currCentr[0][0] + LENGTHS[ii][1]*np.cos(ROTATIONS[ii][1]*np.pi/180)])
           yaxis_y = np.array([currCentr[1][0], currCentr[1][0] + LENGTHS[ii][1]*np.sin(ROTATIONS[ii][1]*np.pi/180)])
           
           # plt.plot(xaxis_x,xaxis_y, 'blue')
           # plt.plot(yaxis_x,yaxis_y, 'red')

           
           out_poly = polygonFromArrayArrays(x,y)
           polygonsINGroup.append(out_poly)
           
           RotsInGroup.append(ROTATIONS[ii][0])
           LensInGroup.append(LENGTHS[ii][0])
    
    
    groupCenter = aauxCenter/len(groupings[currGrouping])/4
    # plt.plot(groupCenter[0][0], groupCenter[0][1],'*', markersize=20)
    # plt.plot(groupCenter[0][0], groupCenter[0][1], str(jj))
    
    LenPeri = len(aauxPerimeterX)
    # aauxPerimeter = np.array([aauxPerimeterX, aauxPerimeterY])
    aauxPerimeter2 = np.zeros((LenPeri,2))
    aauxPerimeter2[:,0] = aauxPerimeterX
    aauxPerimeter2[:,1] = aauxPerimeterY
    # aauxPerimeter2 = np.reshape(aauxPerimeter, (LenPeri,2))
    # aauxPerimeter3 = MultiPoint(aauxPerimeter2)  # this is a shapely POLTYGON
    Perimeter = ConvexHull(aauxPerimeter2)
    pointt = Perimeter.vertices
    perimeterPoints = aauxPerimeter2[pointt,:]
    perimeterPointsPlot = np.vstack((perimeterPoints, perimeterPoints[0,:]))
    
    # plt.plot(perimeterPointsPlot[:,0], perimeterPointsPlot[:,1], linestyle='dashed',color='green', linewidth=0.5)
       
    return RotsInGroup, LensInGroup, perimeterPointsPlot, groupCenter

#%% Find order of rectangles acording to x and y axes
def findOrderVector(ALLBBOXES, groupings, currGrouping, vect):
    VectOrder = []
        
    # From: https://numpy.org/doc/stable/reference/generated/numpy.sort.html
    dtype = [('Pointer', np.int64), ('Length',np.float64)]
    Values = []
    
    for ii in groupings[currGrouping]:
        currPolygon = ALLBBOXES[ii]
        auxCurrCentr = np.array(currPolygon.centroid.coords.xy)
        currCentr = np.reshape(auxCurrCentr, (1,2))
        auxLen = np.dot(currCentr, vect)
        
        Values.append((ii, auxLen[0]))
       
    a = np.array(Values, dtype=dtype)
    a = np.sort(a, order='Length')
    
    for kk in range(a.shape[0]):
        VectOrder.append(a[kk][0])
        
    return VectOrder

#%% Handle 2 rectangles given its vertices
def handle2rectangles(Poly1, Poly2, x_vect, y_vect):
    points_1 = np.zeros((4,2))
    auxx = np.reshape(np.array(Poly1.exterior.coords.xy[0][0:4]), (4,1))
    points_1[:,0] = auxx[:,0]
    auxx = np.reshape(np.array(Poly1.exterior.coords.xy[1][0:4]), (4,1))
    points_1[:,1] = auxx[:,0]
    
    points_2 = np.zeros((4,2))
    auxx = np.reshape(np.array(Poly2.exterior.coords.xy[0][0:4]), (4,1))
    points_2[:,0] = auxx[:,0]
    auxx = np.reshape(np.array(Poly2.exterior.coords.xy[1][0:4]), (4,1))
    points_2[:,1] = auxx[:,0]
    
    sides_1 = np.zeros((4,2)) 
    sides_2 = np.zeros((4,2)) 
    for ii in range(4):
        if ii < 3:
            sides_1[ii,:] = points_1[ii+1,:] - points_1[ii,:]
            sides_2[ii,:] = points_2[ii+1,:] - points_2[ii,:]
        else:
            sides_1[ii,:] = points_1[0,:] - points_1[ii,:]
            sides_2[ii,:] = points_2[0,:] - points_2[ii,:]
    
            
    # verify CCW orientation of sides
    # np.cross(sides_1[0,]/np.linalg.norm(sides_1[0,]), sides_1[1,]/np.linalg.norm(sides_1[1,]))
    # np.cross(sides_1[1,]/np.linalg.norm(sides_1[1,]), sides_1[2,]/np.linalg.norm(sides_1[2,]))
    # np.cross(sides_1[2,]/np.linalg.norm(sides_1[2,]), sides_1[3,]/np.linalg.norm(sides_1[3,]))
    # np.cross(sides_1[3,]/np.linalg.norm(sides_1[3,]), sides_1[0,]/np.linalg.norm(sides_1[0,]))
    # np.cross(sides_2[0,]/np.linalg.norm(sides_2[0,]), sides_2[1,]/np.linalg.norm(sides_2[1,]))
    # np.cross(sides_2[1,]/np.linalg.norm(sides_2[1,]), sides_2[2,]/np.linalg.norm(sides_2[2,]))
    # np.cross(sides_2[2,]/np.linalg.norm(sides_2[2,]), sides_2[3,]/np.linalg.norm(sides_2[3,]))
    # np.cross(sides_2[3,]/np.linalg.norm(sides_2[3,]), sides_2[0,]/np.linalg.norm(sides_2[0,]))
    
    
    coords_1x =  np.zeros((4,1))      
    coords_2x =  np.zeros((4,1)) 
    coords_1y =  np.zeros((4,1))      
    coords_2y =  np.zeros((4,1))
    for ii in range(4):
        coords_1x[ii] = np.dot(x_vect, sides_1[ii,:])
        coords_1y[ii] = np.dot(y_vect, sides_1[ii,:])
        coords_2x[ii] = np.dot(x_vect, sides_2[ii,:])
        coords_1y[ii] = np.dot(y_vect, sides_2[ii,:])
    
    return
#%% Distance point line (defined by 2 points
# From: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
def distPointLine(P0, P1, P2):
    # Point is O0 and line defined by two points P1 and P2
    num = np.linalg.norm((P1[0] - P2[0])*(P1[1] - P0[1]) - (P1[0] - P0[0])*(P2[1] - P1[1]))
    den = np.sqrt((P2[0] - P1[0])**2 + (P2[1] - P1[1])**2)
    dist = num/den   
    return dist

    
    
    
#%% find rectangle coordinates along reference directions
def findCoords(polygon, axis):
    
    
    return

#%% interesection of 2 lines on a plane
# FROM: https://www.cuemath.com/geometry/intersection-of-two-lines/
def line_intersect(P1, P2, T1, T2):
    a1 = P2[1] - P1[1]
    b1 = P1[0] - P2[0]
    c1 = -P1[0]*(P2[1] - P1[1]) + P1[1]*(P2[0] - P1[0])
    
    a2 = T2[1] - T1[1]
    b2 = T1[0] - T2[0]
    c2 = -T1[0]*(T2[1] - T1[1]) + T1[1]*(T2[0] - T1[0])
    # P = np.zeros((1,2))
    # P[0,0] = (b1*c2 - b2*c1)/(a1*b2 - a2*b1)
    # P[0,1] = (c1*a2 - c2*a1)/(a1*b2 - a2*b1)
    
    # NUM1 = b1*c2 - b2*c1
    # NUM2 = c1*a2 - c2*a1
    # DEN = a1*b2 - a2*b1
    
    # if DEN < 1e-10:
    #     P = np.array([np.inf, np.inf])
    # else:
        
    v1 = (P2 - P1)/np.linalg.norm(P2 - P1)    
    v2 = (T2 - T1)/np.linalg.norm(T2 - T1)     
    
    if abs(abs(np.dot(v1,v2)) - 1.0) < 1e-6:
        P = np.array([np.inf, np.inf])
    else:    
        P = np.array([(b1*c2 - b2*c1)/(a1*b2 - a2*b1), (c1*a2 - c2*a1)/(a1*b2 - a2*b1)])
    return P

#%% scale polygon symple
def scalePolySimple(points, scaleFactor):
    ncenter = np.array([sum(points[0:4,0])/4,sum(points[0:4,1])/4])
    # ncenter = np.reshape(center, (1,2))
    newPoints = points - ncenter
    newPoints *= scaleFactor
    newPoints += ncenter
    return newPoints

#%% center polygon
def centerBox(BOX, CENTER):
    for ii in range(BOX.shape[0]):
        BOX[ii,:] -= CENTER
    return BOX
#%% find in list of arrays
def findInlistofArrays(Element, ListOfArrays):
    Location = -1
    for ii in range(len(ListOfArrays)):
        if (Element in ListOfArrays[ii]) == True:
            Location = ii
            break
    return Location
#%% group connected boxes
def groupBoxes(indxs):
    nIndxs = indxs.shape[0]
    # locate isolated buildings
    Isolated = [];
    for ii in range(nIndxs-1):
        if indxs[ii+1,0] - indxs[ii,0] > 1:
            Isolated.append(indxs[ii,0]+1)
    nIsolated = len(Isolated)       # number of isolated buildingsd 
                   
    # ---------- loop through all possible groupings -------
    
    groupings = np.empty((1,), dtype=object)
    groupings[0] = []
    groupings[0].append(indxs[0,0])
    groupings[0].append(indxs[0,1])
    indxs[0] = [-1, -1]
    jj = 0
    Empty = False
    while Empty == False:

        # ---------- inside grouping entry  jj --------------
        Change = 1
        while Change == 1:
            Change = 0
            for ii in range(1, nIndxs):
                # print(ii)
                list1 = groupings[jj,]
                list2 = indxs[ii,:]
                # print(list2)
                ress = list(set(list1).intersection(list2))         
                if len(ress)>0:
                    for kk in range(len(indxs[ii,:])):
                        toBeAppended = indxs[ii,kk]
                        curList = groupings[jj]
                        if (toBeAppended in curList) == False:
                            groupings[jj].append(toBeAppended)
                            Change = 1
                    indxs[ii] = [-1, -1]
                # a=1
        # end -------- inside grouping entry  jj --------------
        
        Empty = True
        for mm in range(nIndxs):
            if indxs[mm,0] != -1:
                Empty = False
                break        
        if Empty == False:    
            jj = jj + 1   
            groupings = np.append(groupings, np.empty((1,), dtype=object))
            groupings[jj] = []
            groupings[jj].append(indxs[mm,0])
            groupings[jj].append(indxs[mm,1])
            indxs[mm] = [-1, -1]
    # end ---------- loop through all possible groupings ------- 
    
    nGroupings = groupings.shape[0] + 1  # number of building groupings
    
    for ii in range(nIsolated):
        jj = jj + 1
        groupings = np.append(groupings, np.empty((1,), dtype=object))
        groupings[jj] = []
        groupings[jj].append(Isolated[ii])

    return groupings

#%% -------------------------------------------------------------------
# losely based on  https://stackoverflow.com/questions/45508202/how-to-rotate-a-polygon
def rotatePolygon(polygon, degrees):
    """ Rotate polygon the given angle about its center. """
    theta = radians(degrees)  # Convert angle to radians
    cosang, sinang = cos(theta), sin(theta)
    
    points = polygon.exterior.coords.xy

    n = len(points[0])
    cx = polygon.centroid.x
    cy = polygon.centroid.y

    new_pointsX = arr.array('d',[])
    new_pointsY = arr.array('d',[])
     
    auxRotPolygon = np.empty((n,2))
    
    for ii in range(n):     
        x = points[0][ii]
        y = points[1][ii]
        tx, ty = x-cx, y-cy
        # right hand rule around z-axis
        new_x = (tx*cosang - ty*sinang) + cx
        new_y = (tx*sinang + ty*cosang) + cy
        new_pointsX.append(new_x)
        new_pointsY.append(new_y)
        
        auxRotPolygon[ii,0] = new_x
        auxRotPolygon[ii,1] = new_y
        
    new_points = (new_pointsX,new_pointsY)     # this is a tuple with two arrays
    rotPolygon = Polygon(auxRotPolygon)
    

    return new_points, rotPolygon
#%% -------------------------------------------------------------------
def removeTotContainedPols(gdf):
    
    nPolys = gdf.shape[0] 
    auxArray = np.zeros((nPolys,2), dtype=int)
    
    Counter = int(0)
    for ii in range(nPolys):
        print(ii)
        for jj in range(nPolys):
            if ii != jj:
                AA = gdf.iloc[ii].geometry
                BB = gdf.iloc[jj].geometry
                if AA.contains(BB) == True:
                    auxArray[Counter,0] = ii
                    auxArray[Counter,1] = jj
                    Counter += 1 
    Contained = auxArray[0:Counter-1,:]  
    return Contained
#%% ---------------------------------------------------------------------
def overlayMatrix(gdf):
    # OverLayMatrix = 0
    nPolys = gdf.shape[0]
    OverLayMatrix = np.zeros((nPolys,nPolys), dtype=np.int0)
    # OverLayMatrix = np.chararray((nPolys,nPolys))
    # OverLayMatrix = np.ndarray((nPolys,nPolys), dtype=object)
    for ii in range(nPolys):
        for jj in range(nPolys):
            if (ii != jj):
                if gdf.iloc[ii].geometry.intersects(gdf.iloc[jj].geometry) == True:
                    OverLayMatrix[ii,jj] = 1
    return OverLayMatrix
#%% -----------------------------------------------------------------
def filterNonPolygonsOut(gdf):
    # From: https://gis.stackexchange.com/questions/425326/delete-rows-in-geopandas-geodataframe-object-based-on-geometry-shapely-geometry
    gdf = gdf[gdf.geom_type != 'Point']
    return gdf

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
    
    
    UTMcenterX = (eUTM + wUTM)/2
    UTMcenterY = (nUTM + sUTM)/2
    sceneCENTER = np.array((1,2), dtype=float)
    sceneCENTER[0] = UTMcenterX
    sceneCENTER[1] = UTMcenterY
    
    # gdf_proj = gdf.to_crs(epsg=32630)
    print('  ')
    print(gdf.crs)
    print(gdf_proj.crs)
    print('   ')
    
    return gdf_proj, gdf, sceneCENTER    # remove one of projected gdfs

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
    
    Lengths = np.empty((4,1))
    Vector =  np.empty((4,3))
    Vector[:,2] = 0.0
    normVector = Vector
    for ii in range(4):
        Vector[ii,0] = xR[ii+1] - xR[ii]
        Vector[ii,1] = yR[ii+1] - yR[ii]
        Lengths[ii] = np.linalg.norm(Vector[ii,:])
        normVector[ii,:] = Vector[ii,:]/Lengths[ii]
        
    # testt = np.cross(normVector[0,:], normVector[1,:])
    # if testt[2] < 0.0:
    #     xR = np.flip(xR,0)
    #     yR = np.flip(yR,0)
    #     a=1
    # for ii in range(4):
    #     Vector[ii,0] = xR[ii+1] - xR[ii]
    #     Vector[ii,1] = yR[ii+1] - yR[ii]
    #     Lengths[ii] = np.linalg.norm(Vector[ii,:])
    #     normVector[ii,:] = Vector[ii,:]/Lengths[ii]
        
    if Lengths[0] >= Lengths[1]:
        Rot0 = np.arctan2(normVector[0,1], normVector[0,0])*180/np.pi
        # Rot1 = np.arctan2(normVector[1,1], normVector[1,0])*180/np.pi
        Rot1 = Rot0 + 90
        Len0 = Lengths[0]
        Len1 = Lengths[1]
        # XX = np.array([np.cos(Rot0), np.sin(Rot0), 0])
        # YY = np.array([np.cos(Rot1), np.sin(Rot1), 0])
        # testt = np.cross(XX, YY)
        # a = 1
    else:
        Rot0 = np.arctan2(normVector[1,1], normVector[1,0])*180/np.pi
        # Rot1 = np.arctan2(normVector[0,1], normVector[0,0])*180/np.pi
        Rot1 = Rot0 + 90
        Len0 = Lengths[1]
        Len1 = Lengths[0]
        # XX = np.array([np.cos(Rot0), np.sin(Rot0), 0])
        # YY = np.array([np.cos(Rot1), np.sin(Rot1), 0])
        # testt = np.cross(XX, YY)
        # a = 1
    if Rot0 > 90 and Rot0 <=180:
        Rot0 = Rot0 - 180
        Rot1 = Rot0 + 90
    if Rot0 >= -180 and Rot0 <= -90:
        Rot0 = Rot0 + 180
        Rot1 = Rot0 + 90
        
    Rot = np.array((2,1), dtype=float)
    Rot[0] = Rot0
    Rot[1] = Rot1
    Len = np.array((2,1), dtype=float)
    Len[0] = Len0
    Len[1] = Len1
        
    return MINRECT, Rot, Len

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

#%% ---------------------------------------------------------------------

# def polygonFromArrayArrays(x,y):
    
#     # based on https://www.matecdev.com/posts/shapely-polygon-from-points.html
    
#     import shapely.geometry as sg
    
#     points = []
#     for ii in range(len(x)):
#         currPnt = sg.Point(x[ii], y[ii])
#         points.append(currPnt)
        
#     out_poly = sg.Polygon([[p.x, p.y] for p in points])
    
#     return out_poly














# def polygonFrom2array_arrays(x,y):
#     # x and y are two array.array same length inputs containing 
#     # the output closed polygon vertices
    
#     import shapely.geometry as sg
    
#     points = []
#     for ii in range(len(sx)):
#         currPnt = sg.Point(x[ii], y[ii])
#         points2.append(currPnt)
        
#     poly2 = sg.Polygon([[p.x, p.y] for p in points2])
    
    
    
    
    
#     return poly_out