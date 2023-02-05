import numpy as np

ref_x = np.array([1, 0, 0])
ref_y = np.array([0, 1, 0])
ref_z = np.array([0, 0, 1])   #  up pointing unitary vect

ref_x2D = np.array([1, 0])
ref_y2D = np.array([0, 1])

#

minSide = 1.10;   # m, max distance to merge points
minAnglDeg = 1.50   # deg, max angle to merge aligned sides


# configure the inline image display

fp = "./{img_folder}/madrid_bldgs.{extension}"
img_folder = "images"
extension = "png"
size = 240

# specify that we're retrieving building footprint geometries
tags = {"building": True,
        "building:levels": True
        } 


# # Location coords
# centerOfScene = (40.410439, -3.704929)   # deg, deg geographic coords
# radiusOfScene = 50   # m


# BBOX = (40.4127, 40.41249, -3.70431, -3.70430)
# BBOX = (40.41337, 40.41251, -3.70383, -3.70445)
# BBOX = (40.41435, 40.41251, -3.70404, -3.70445)

# BBOX = (40.41398, 40.41185, -3.69847, -3.70236)   # Center of Madrid
BBOX = (51.50871, 51.50553, -0.13112, -0.13739)   # Center of London
# BBOX = (North_deg, South_deg, East_deg, West_deg)


bboxHeight = 5  # m
triangHeight = 10 # m


