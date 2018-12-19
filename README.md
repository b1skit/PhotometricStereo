# Photometric Stereo Recovery of Shape & Relative Distance

This project uses photometric stereo to recover the surface normal and overall shape of an unknown object from sets of 3 reference images.
It is written from scratch in MATLAB. An estimation of the light position for each reference image is used to generate an inverse lookup table that indexes gradient-space (p, q) coordinates as a function of image intensities. The program uses the lookup table to efficiently calculate an object’s surface orientation in terms of per-pixel gradients and normal directions. The program then performs integration along a path to recover the relative distance of each pixel in the image, completing an overall description of the object’s shape.

### Instructions:

Download the repository, and execute the **ShapeRecovery.m** file within MATLAB.
