from os import close
import sys
import math

# This function reads the .data files and returns a 2D 512 * 512 array 
# that contains the heightmap
def read_file (filename):
    my_file = open(filename, "rb")
    arr = [[] for _ in range(512)]
    byte = my_file.read(1)
    counter = 0
    while byte:
        arr[counter // 512].append(byte[0])
        byte = my_file.read(1)
        counter += 1
    return arr

# Vector class makes dealing with vectors and positions easier
class Vector (object):
    def __init__(self, vec):
        self.vec = vec
        self.dim = len(vec)

    # Returns the length of the vector
    def length (self):
        total_length = 0
        for comp in self.vec:
            total_length += math.pow(comp, 2)
        return math.sqrt(total_length)
    
    # Returns a normalized version of the vector
    def normalize (self):
        vec_curr = []
        vec_len = self.length()
        for comp in self.vec:
            vec_curr.append(comp/vec_len)
        return Vector(vec_curr)

    # Returns the difference between two vectors
    def __sub__ (self, other):
        diff = []
        for i in range(self.dim):
            self_comp = self.vec[i]
            other_comp = other.vec[i]
            diff.append(self.vec[i] - other.vec[i])
        return Vector(diff)

    # Returns the sum of two vectors
    def __add__ (self, other):
        summation = []
        for i in range(self.dim):
            self_comp = self.vec[i]
            other_comp = other.vec[i]
            summation.append(self.vec[i] + other.vec[i])
        return Vector(summation)
    
    # Scalar multiplication - only works when right-multiplying
    def __mul__ (self, other):
        prod = []
        for comp in self.vec:
            prod.append(other * comp)
        return Vector(prod)

    # Returns components of the vector
    def __getitem__(self, key):
        return self.vec[key]
    
    # Returns a string version of the vector
    def __str__ (self):
        return str(self.vec)

# Finds the distance in units of the (dir) vector until
# we intersect (axis)
def dist (origin, dir, axis):
    dist = math.floor(origin[axis]) - origin[axis] + 1
    return dist/dir[axis] 

# Checks which axis intersects first along the (dir) vector
# from the (origin) position
def gridmarch (origin, dir):
    x_dist = dist(origin, dir, 0)
    y_dist = dist(origin, dir, 1)
    if (x_dist < y_dist):
        return 0
    else:
        return 1

# Linear interpolation - used for heights
def lerp (a, b, u):
    return (1 - u) * a + u * b

# Solution function, takes in a start point and an end point in format (x, z)
# also takes in the heighmap (data), the spatial resolution, and the height value
def solution (startpoint, endpoint, data, spatial_res, height_val):

    # Setting the current origin and direction (trajectory) for our motion
    pt_a = Vector(startpoint)
    pt_b = Vector(endpoint)
    trajectory = pt_b - pt_a
    trajectory = trajectory.normalize()
    curr_origin = pt_a

    # Variable to be returned, holding the distance so far
    curr_dist = 0

    # Holds the previous closest direction (see line 112)
    previous_direction = -1

    # Continues until we reach the general vicinity of the endpoint
    while (pt_b - curr_origin).length() > 0.0001:

        # We find the axis that we'll intersect soonest based on our trajectory
        closest_direction = gridmarch(curr_origin, trajectory)

        # The amount by which we need to travel along our trajectory vector to intersect
        # an axis
        scalar = dist(curr_origin, trajectory, closest_direction)

        # Vector representing our movement
        delta = trajectory * scalar

        # The current target position
        target = curr_origin + delta

        # Finding the height at our origin position by linearly interpolating
        ya_pt_a = data[math.floor(curr_origin[1])][math.floor(curr_origin[0])]
        if (previous_direction == 0):
            ya_pt_b = data[math.floor(curr_origin[1])][math.ceil(curr_origin[0])]
        else:
            ya_pt_b = data[math.ceil(curr_origin[1])][math.floor(curr_origin[0])]
        ya_u = curr_origin[closest_direction] % 1
        height_a = lerp(ya_pt_a, ya_pt_b, ya_u)
            
        # Finding the height at our target position by linearly interpolating
        yb_pt_a = data[math.floor(target[1])][math.floor(target[0])]
        if (closest_direction == 0):
            yb_pt_b = data[math.floor(target[1])][math.ceil(target[0])]
        else:
            yb_pt_b = data[math.ceil(target[1])][math.floor(target[0])]
        yb_u = target[closest_direction] % 1
        height_b = lerp(yb_pt_a, yb_pt_b, yb_u)

        # Finding the distance traveled and adding that to our distance so far
        delta_y = height_b - height_a
        delta_3d = Vector([delta[0], delta_y ,delta[1]])
        delta_3d_metric = Vector([delta_3d[0] * spatial_res, delta_3d[1] * height_val, delta_3d[2] * spatial_res])
        curr_dist += delta_3d_metric.length()
        previous_direction = closest_direction
        curr_origin += trajectory * scalar
    return curr_dist

start = (int(sys.argv[1]), int(sys.argv[2]))
end = (int(sys.argv[3]), int(sys.argv[4]))
data = read_file(sys.argv[5])
spatial_res = 30
height_val = 11
soln = solution(start, end, data, spatial_res, height_val)

print("The surface distance between " + str(start) + " and " + str(end) + " is " + str(soln))