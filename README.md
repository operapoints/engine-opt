# Building
This script uses the pagmo optimization library's gaco optimizer to optimize a small 70N jet engine for Isp. To build run this in the terminal:
```
mkdir build
cd build
cmake ..
cmake --build .
```
If pagmo does not work on your system just replace the install. I used the spack package manager and it worked fine. To do this, install spack if you don't have it, navigate to the project dir (the one above the build dir) and then run:
```
spack install pagmo2
cp -R $(spack location -i pagmo2) ./pagmo2
```
Then, cd back into the build directory and build as normal.