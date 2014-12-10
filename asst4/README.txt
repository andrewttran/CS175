Roseanne Feng
Belinda Zeng
Assignment 4


For the files we submitted, we really only modified picker.cpp, scenegraph.cpp and asst4.cpp. All other files were included in last weeks submission and/or were included in the starter code for this week's pset. Our platform we used are Windows laptops, and instructions for compiling and running the code are extremely standard.

We're pretty sure we met all the problem set requirements, including implementing the scene graph, the part picking, the transformation and the translation fix. 

In terms of the code design for the two major tasks:
 * We implemented picking in picker.cpp, following the formulas discussed in lecture, and used the picker visitor to assign ids/colors to objects whichallowed the user to select and activate them.
 * For the transformation hierarchy, we worked within the scenegraph.cpp to get the accumulated rbt. We then rewrote our previous transformations (in the motion function) to happen with respect to the parent frame rather than the world frame. 


Running the program is extremely standard as well. To reiterate, as per the spec, you can use the 'p' key to indicate that the next click means picking an object and the 'v' key to change the view point. For the mouse button usage, you can click and drag the arcball to rotate the the part of the robot. The right button is for x and y translations and the middle button is for translation in the z direction. 