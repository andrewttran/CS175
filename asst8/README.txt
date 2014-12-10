Roseanne Feng
Belinda Zeng
Assignment 8


For the files we submitted, we really only modified asst8.cpp. All other files were included in last weeks submision and/or were included in the starter code for this week's pset. 

Our platform we used are Windows laptops, and instructions for compiling and running the code are extremely standard.

We're pretty sure we met all the problem set requirements, including working with the bunny mesh, implementing straight fur, animating the fur, and combing the hair

In terms of the code design:
 * For working with the cube mesh, we implemented a function that uploads a Mesh to a SimpleGeometryPN
 * To implement straight fur, we modified updateShellGeometry to calculate n,p,s, and t values which we then upload to the mesh. We left in a function called updateShellGeometryProto which is not called by the program itself at any time, but can be used to check correctness of Task 1 and 2.
 * To animate the fur, we calculated the force acting on the tip (which can be broken down into gravity and spring force), converted s to the appropriate frame, then updated the tip position.
 * To create curved fur we calculated the constant d value, which we then used to find the position of each point. We also then adjusted the value of the normal.

In terms of running the program, as per the spec, you can use the left/right arrow to change fur length, the up/down arrow to change hairiness, the 'f' key to toggle between smooth and flat shading, the '0' key to increase the number of subdivision steps, the '9' key to decrease the number of subdivision steps, the '7' key to halve the speed of deformation, the '8' key to double the speed of deformation, the 'c' key to copy current key frame RBT data, the 'u' key to update, the '>' key to advance to the next frame, the '<' frame to go to the previous frame, the 'd' key to delete, the 'n' to create a new key frame, the 'i' to input key frames from the input file, the 'w' key to output key frames to the output file, the 'y' key to play/stop the animation, the '+' key to make the animation go faster, and the '-' key to make the animation go slower. You can now also move the light as you would move any other object in the scene.

