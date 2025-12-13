# Building
This script uses the pagmo optimization library's gaco optimizer to optimize a small 70N jet engine for thrust and Isp. To build run this in the terminal:
```
mkdir build
cd build
cmake ..
cmake --build .
```
I installed pagmo2 using spack, which worked:
```
spack install pagmo2
```
Then, cd back into the build directory and build with:

```
spack load pagmo2
```