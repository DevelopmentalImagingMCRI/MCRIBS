MIRTK Deformable Package
========================

The Deformable module of the Medical Image Registration ToolKit (MIRTK) is a library
for the Euler integration of deformable meshes such as cortical surfaces. The MIRTK
[deform-mesh](https://mirtk.github.io/commands/deform-mesh) command
can be used to deform an initial mesh such as the convex hull of an input segmentation
or a bounding sphere based on internal and external point set/surface forces.
The integration is stopped when a suitable stopping criterion is fullfilled such as
a fixed number of iterations, target objective function value, or surface smoothness
(e.g., for cortical surface inflation). The internal forces can further be utilized by
the MIRTK Registration module to constrain the transformation, for example, to constrain
the cortical surface to remain smooth after transformation.

See [online documentation](https://mirtk.github.io/modules/deformable)
for more.


License
-------

The MIRTK Deformable module is distributed under the terms of the
[Apache License Version 2](http://www.apache.org/licenses/LICENSE-2.0).
