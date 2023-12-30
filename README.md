setup

to install OpenGL for plotting on Linux

```sudo apt-get install libglfw3-dev```

also install Eigen library, might be needed to change header in main.cpp

```#include <eigen3/Eigen/Eigenvalues>```

to the location of Eigen, might add this liberary to source files in the future

now what should work to build and run the program

```
mkdir -p build
cd build
cmake ..
make
./MES
```

