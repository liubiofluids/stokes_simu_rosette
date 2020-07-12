# stokes_simu_rosette

Hydrodynamic simulation of a self-propelled sphere (by a helical flagellum) near a solid-liquid interface, using the resistance matrices.

To generate the 3D velocity and angular velocity in a body-fixed frame of reference (with flagellum aligned in x), use: [vV, vOmega]=propAnisoCell(-pi, mrstc, a, pangle); where "mrstc" is a resistence matrix of the sphere, "a" is the radius of the sphere, "pangle" is the deviation of the flagellar axis from the radial direction. 

To generate a movie for the simulation, use the script: mov06_trace_rosette.
