# mandelbulb
Compute mandelbul (3d fractal inspired from Mandelbrot set) animations
usage: mandelbulb.exe <size> <output> [nb]
ex: mandelbulb.exe 512 mb1 64
Save a bmp image or sequence of bmp images you can use with external tool like VirtualDub (free) to create a movie.

The original code is js code from Roy van Rijn : https://github.com/royvanrijn/mandelbulb.js/blob/master/mandelbulb.html

tip: compile with C/C++ > Optimization > Fast /Ot

mandelbulb.ini file contains all parameters, static (xxx) or dynamic (xxxDelta) ex:
viewAngle 150
viewAngleDelta 2
lightAngle 170
# lightAngleDelta 3
MAX_ITER 100
eyeDistanceFromNearField 2.2
# eyeDistanceFromNearFieldDelta 1.0
DEPTH_OF_FIELD 2.5
Power 8.0
# PowerDelta 0.1
shade 1
shadowOn 1
