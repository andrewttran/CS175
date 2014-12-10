Roseanne Feng
Belinda Zeng
Assignment 6


For the files we submitted, we really only modified asst6.cpp and normal-gl3.fshader. All other files were included in last weeks submision and/or were included in the starter code for this week's pset. 

Our platform we used are Windows laptops, and instructions for compiling and running the code are extremely standard.

We're pretty sure we met all the problem set requirements, including implementing the new operations dealing creating a smoother animation, changing the material infrastructure to deal with uniforms, make the light movable, and add texture bump mapping for the ground. 

In terms of the code design:
 * For our Catmull-Rom interpolation, we modified slerp and lerp. We evaluated Bezier curves at intermediate times to display intermediate frames. 
 * In terms of modifying our material infrastructure, we simply followed the instructions from asst6_5-snippets.cpp to deal with the uniforms class.
 * For the bump map texture, we modified the fragment shader to read the appropriate pixel from the texture, transform it eye space, and use it for shading calculation. 


In terms of running the program, as per the spec, you can use the 'c' key to copy current key frame RBT data, the 'u' key to update, the '>' key to advance to the next frame, the '<' frame to go to the previous frame, the 'd' key to delete, the 'n' to create a new key frame, the 'i' to input key frames from the input file, the 'w' key to output key frames to the output file, the 'y' key to play/stop the animation, the '+' key to make the animation go faster, and the '-' key to make the animation go slower. You can now also move the light as you would move any other object in the scene.