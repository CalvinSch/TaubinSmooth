### Calvin Schaul

### CS657/718 Advanced Computer Graphics and Animation

Code written in the `taubin_smooth.cpp` file. Based on smoothing.cpp demo code

Using 100 iterations of Taubin Smoothing, the result is seen in `output.off`

The noisy, "rough" mesh is included as a `.off` file as well for further testing.

- Compile code with GCC
`gcc taubin_smooth.cpp`
- Run the exe
`a.exe` or `./a.out` depending on windows/linux

See the 1995 Research Paper by Gabriel Taubin https://graphics.stanford.edu/courses/cs468-01-fall/Papers/taubin-smoothing.pdf
See explicit steps for vector averaging and delta p calculations for Taubin Smoothing http://crl.ethz.ch/teaching/shape-modeling-18/lectures/07_RemeshingSmoothing.pdf
    - Slide 103
