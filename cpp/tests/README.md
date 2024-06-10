# C++ Testing

To set up testing on x86:
```
cmake -S . -B build
```

on Arm:
```
cmake -S . -B build -DTARGET_ARM=ON
```

Then to build and run tests:
```
cd build
cmake --build .
ctest
```

