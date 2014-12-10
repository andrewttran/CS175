Roseanne Feng
Belinda Zeng
Assignment 7


For the files we submitted, we modified asst7.cpp and geometry.h All other files were included in last weeks submision and/or were included in the starter code for this week's pset. 

Our platform we used are Windows laptops, and instructions for compiling and running the code are extremely standard.

We're pretty sure we met all the problem set requirements, including working with the cube mesh, implementing smoothshading, anmating the vertices of the cube, and implementing subdivision.

In terms of the code design:
 * For working with the cube mesh, we implemented a function createMeshToSimple that creates a SimpleGeometryPN from a Mesh
 * To implement smooth shading, we implemented a function called average_normals which computed an "average normal" for each vertex by averaging the normals of all the faces incident to the vertex in the mesh, normalizing to length 1
 * For animating the vertices, we simply got the position of the vertices, multiplied each position by a sin function that takes in time t as a variable, then called GlutPostRedisplay()
 * To implement subdivision, we looped over the faces, calculating faceVertex values, then setting the new face vertices, looped over the edges, calculating edgeVertex values, then setting the new edge vertices, looped over the vertices, caluclating vertexVertex values, then setting the vertex vertices, and finally calling subdivision. 

Right now we've capped subdivision at 6, meaning we cap the subdivision steps at 7 including the original mesh (the cube).

In terms of running the program, as per the spec, you can use the 'f' key to toggle between smooth and flat shading, the '0' key to increase the number of subdivision steps, the '9' key to decrease the number of subdivision steps, the '7' key to halve the speed of deformation, the '8' key to double the speed of deformation, the 'c' key to copy current key frame RBT data, the 'u' key to update, the '>' key to advance to the next frame, the '<' frame to go to the previous frame, the 'd' key to delete, the 'n' to create a new key frame, the 'i' to input key frames from the input file, the 'w' key to output key frames to the output file, the 'y' key to play/stop the animation, the '+' key to make the animation go faster, and the '-' key to make the animation go slower. You can now also move the light as you would move any other object in the scene.

