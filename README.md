# Cesium 3D Project

In order to solve the problem, I chose to think of the surface distance as an arc length and used numerical integration in order to find the distance. Height intervals were determined using grid marching and linear interpolation.

To run the program use the following:

`python3 solution.py [STARTPOINT.x] [STARTPOINT.z] [ENDPOINT.x] [ENDPOINT.z] [.data FILE]`
