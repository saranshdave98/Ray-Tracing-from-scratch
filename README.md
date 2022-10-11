# Ray Tracing from scratch
 This is an implementation of a basic Ray-Tracing from algorithm scratch using pure C++. No graphics libraries whatsoever were used while building this project. It was created using only the basic C++ header files iostream, fstream and cmath.

# Optimization
 As the object shapes used in this implementation are basic objects like pyramids, spheres and cubes, the use of Barycentric Coordinate computations was replaced with basic cross product and dot product calculations to avoid heavy operations like matrix inversion. This results in very fast rendering only within a few seconds (2-6 seconds on my PC) without any parallelization whatsoever.
 
