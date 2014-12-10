Roseanne Feng
Belinda Zeng
Assignment 5


List of all files we've modified include:
- asst5.cpp
All other files submitted were exactly the same as the files given in the starter code. 

We've both done our work on Windows, and the compile instructions are standard. We haven't changed any of that. 
We believe we've met all problem set requirements, as per the specification. We've also compared to the solution binaries and believe that our code runs with the same result. 

In terms of code design, we've implemented hot keys, interpolation, and animation. Hot keys were implemented as described in Tasks 1 and 3. In terms of the interpolation, we simply interpolated the translation component of the RBTs with lerp and the rotation component with slerp, and animation as outlined in Task 3.

In terms of running the program, as per the spec, you can use the 'c' key to copy current key frame RBT data, the 'u' key to update, the '>' key to advance to the next frame, the '<' frame to go to the previous frame, the 'd' key to delete, the 'n' to create a new key frame, the 'i' to input key frames from the input file example.txt, the 'w' key to overwrite example.txt with the current key frames, the 'y' key to play/stop the animation, the '+' key to make the animation go faster, and the '-' key to make the animation go slower.