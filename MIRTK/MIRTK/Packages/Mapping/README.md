MIRTK Mapping Package
=====================

The Mapping module of the Medical Image Registration ToolKit (MIRTK) is a library
for the mapping and re-parameterization of surfaces and volumes. It provides the
following commands for computing such maps:

- [calculate-boundary-map](https://mirtk.github.io/commands/calculate-boundary-map):
  Compute boundary map for fixed-boundary surface maps.
- [calculate-surface-map](https://mirtk.github.io/commands/calculate-surface-map):
  Compute map value for each (interior) point of the surface.
- [calculate-volume-map](https://mirtk.github.io/commands/calculate-volume-map):
  Compute map value for each (interior) point of the volume.

with corresponding tools to evaluate the distortions resulting from surface
and volume reparameterizations computed with these tools:

- [evaluate-surface-map](https://mirtk.github.io/commands/evaluate-surface-map):
  Evaluate distortions and other quality measures given a surface parameterization.
- [evaluate-volume-map](https://mirtk.github.io/commands/evaluate-volume-map)
  Evaluate distortions and other quality measures given a volume parameterization.

See [online documentation](https://mirtk.github.io/modules/mapping)
for more.


License
-------

The MIRTK Mapping module is distributed under the terms of the
[Apache License Version 2](http://www.apache.org/licenses/LICENSE-2.0).
