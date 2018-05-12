# Bounding Volume Hierarchy
- medbvh: median based bounding box hierarchy  
- binbvh: binning based bounding box hierarchy
- medqbvh: median based quad bounding box hierarychy [1]  
- binqbvh: binning based quad bounding box hierarchy [1]  
- gridqbvh: grid based quad bounding box hierarchy.

# Results

## Sibenik
[Sibenik by Marko Dabrovic](http://hdri.cgtechniques.com/~sibenik2/download/), consists of 75284 triangles.

|method|depth|build time|render time|rays per second|
|:---|---:|---:|---:|---:|
|medbvh|14|0.082720|1.767554|173799|
|binbvh|18|0.157188|0.507048|605860|
|medqbvh|8|0.041134|0.398104|771657|
|binqbvh|11|0.159771|0.310184|990381|
|gridqbvh|12|0.099577|0.670912|457884|

## Conference room
Conference room model by Greg Ward, consists of 331179 triangles.

|method|depth|build time|render time|rays per second|
|:---|---:|---:|---:|---:|
|medbvh|16|2.594485|2.510323|122374|
|binbvh|21|2.767390|2.259608|135952|
|medqbvh|9|0.147858|1.105135|277975|
|binqbvh|15|1.040235|0.754667|407067|
|gridqbvh|12|0.148544|0.496506|618723|

# References
1. Holger Dammertz, Johannes Hanika, Alexander Keller, “Shallow bounding volume hierarchies for fast SIMD ray tracing of incoherent rays”, EGSR 2008.
2. Danilewski, Piotr & Popov, Stefan & Slusallek, Philipp, “Binned SAH Kd-Tree Construction on a GPU”, 2018.
3. Karras, Tero., “Maximizing Parallelism in the Construction of BVHs, Octrees, and k-d Trees”, High Performance Graphics 2012.
4. Jacopo Pantaleoni, David Luebke, “HLBVH: Hierarchical LBVH Construction for Real-Time Ray Tracing”, High Performance Graphics 2010.
5. Kirill GaranzhaEmail, Simon Premože, Alexander Bely, Vladimir Galaktionov, “Grid-based SAH BVH construction on a GPU”, 2011.

# License
This is free and unencumbered software released into the public domain.
