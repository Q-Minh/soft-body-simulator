# Soft Body Simulator for Fast Prototyping

## Dependencies

- [CMake](https://cmake.org)
- [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [glm](https://github.com/g-truc/glm)
- [JSON for Modern C++](https://github.com/nlohmann/json)
- [OpenGL 3.3+](https://www.opengl.org//)

To download and install the dependencies, open your preferred terminal in 
administrator mode and run:
```
$ git clone https://gitlab.com/libeigen/eigen
$ cd eigen
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --target install --config Release
$ cd ..
$ git clone https://github.com/g-truc/glm
$ cd glm
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --target install --config Release
$ cd ..
$ git clone https://github.com/nlohmann/json
$ cd json
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --target install --config Release
```

To install Debug binaries/artifacts, replace all "Release" arguments with "Debug".

## Building and installing

```
$ cd path/to/sbs
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --target sbs --config Release
$ cmake --build build --target install --config Release
```

Again, replace "Release" by "Debug" for debug builds and installs.
