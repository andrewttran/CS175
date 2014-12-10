Roseanne Feng
Belinda Zeng
Assignment 3


For the files we submitted, we really only modified rigtform.h and asst3.cpp. All other files were included in last week's submission and/or were included in the starter code for this week's pset. 

Our platform we used is Windows laptops, and instructions for compiling and running the code are extremely standard.

We're pretty sure we met all the problem set requirements, including implementing the new operations dealing with the RBTs, the arcball rotation, and the translation fix. 

In terms of the code design:
 * We did all the rigid body transforms in rigtform.h, following the formulas discussed in lecture and the textbook. 
 * For the arcball, we worked within the motion function to detect where the user has clicked. We then used GetScreenSpaceCoord() to find the center of the arcball in screen space coordinates, and called a helper function to calculate the "z" value for both points. We then calculated the vectors, normalized them, and used them to build the quaternion representing our desired rotation. Arcball radius and scale were handled as discussed in the specs, using globals g_arcballScreenRadius and g_arcballScale. 
 * For the translation fix, we used the g_arcballScale to create a scale matrix, and multiply that by the MVM before drawing the sphere in the DrawStuff() function.

Running the program is extremely standard as well. To reiterate, as per the spec, you can use the 'm' key to change frames, the 'o' key to change the active object, and the 'v' key to change the view point. For the mouse button usage, you can click and drag the arcball to rotate the cube. The right button is for x and y translations and the middle button is for translation in the z direction. 

One note is that during any user interaction besides translation, the getScreenSpaceCoord throws a warning, but we're not sure how to address that. Otherwise all operations perform as normal.