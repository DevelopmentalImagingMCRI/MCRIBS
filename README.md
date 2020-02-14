Trying to fix the vertex ordering issue in MIRTK.

When new surfaces are created by merging/cutting in Deformable the vertices are reordered. This branch will try to put vertex IDs into the point data arrays within the VTP files in deformable so that we can map the pial surface order to the white surface order.
