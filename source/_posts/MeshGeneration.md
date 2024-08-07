---
title: Mesh Generation
date: 2024-05-18 20:02:29
tags: [Mesh]
---

## Overview

On-demand 3D terrain generation in games is traditionally accomplished through several methods. One is voxelization which quantifies terrain as regularly sized voxels(cubes) and the other is polygonization which quantifies terrain as combinations of 2D polygons. 
<!-- more -->

The prior has its advantages in being straightforward to interact and generate but innevitably results in the choppy/cubic terrain. The term voxel is itself a 3D extention of the 2D pixel and similar to its counterpart, there are strategies such as sticky-cubes to counteract this limitation. However, further abstraction often leads to increased memory usage and reduced intuitivity. 

Polygon terrrain, meanwhile often requires more consideration in generation and modification. It is often lauded for being less intuitive then voxelization especially when modifying terrain. However, the terrain generated is smooth unlike its counterpart and, given the proper implementation, it has the potential for smoother and more natural interaction. 

Voxelization is a process explored thoroughly in the gaming industry, with many games, most notably of which <i>Minecraft</i>, employing the system in generation of completely voxelized terrain. Comparatively, though procedural polygon generation has gained much attention, it still lacks as much research as voxel generation especially in the area of liquid representation. Through this post, I highlight the considerations and solutions in generating a smooth procedural terrain.

## Background

The Marching Cubes Algorithm is a well-known algorithm for converting an implicit surface into a renderable and manageable mesh. An implicit surface refers to a surface defined by any generating function wherein the surface exists wherever the function evaluates to zero. 

For instance the implicit surface of a sphere would have a generating function of the form:
f(x,y,z) = (x + x<sub>1</sub>)<sup>2</sup> + (y + y<sub>1</sub>)<sup>2</sup> + (z + z<sub>1</sub>)<sup>2</sup>

In fact, evaluating each object as a generating functions in rendering is indeed a valid method of rendering primitives. However, often complex models necessitates the creation of complicated surfaces that are difficult to render through functions or collections of primitives. As such, the existence of meshes stands as a well researched and standardized way of rendering and handling implicit surfaces.

The generating function from which Marching-Cubes may extract a mesh from is represented as a scalar field. If each scalar-value contains information concering its specific point in space, such a representation allows for regionalized evaluation of the surface with no overlap as long as the nearest samples may be identified. In this case a grid representation is most common, where taking regular samples along grid corners allows for regions to be evaluated as cubes. Creating mesh components individually from each isolated cube, these cubes may be combined to produce cohesive meshes. Further detail regarding Marching-Cubes may be found [here](https://paulbourke.net/geometry/polygonise/).


## Data Representation

### Scalar Fields

To obtain a scalar field given a simple function, it is trivial to evaluate the function at regular intervals along the grid. Returning to the sphere example, given the function: 
f(x,y,z) = (x - 4)<sup>2</sup> + (y - 4)<sup>2</sup> + (z - 4)<sup>2</sup> - 4

We are able to extract such values by simply evaluating the functions at each integer coordinate 

|_| 0| 1| 2| 3| 4| 5| 6| 7| 8|
|-|--|--|--|--|--|--|--|--|--|
|0|28|21|16|13|12|13|16|21|28|
|1|21|14| 9| 6| 5| 6| 9|14|21|
|2|16| 9| 4| 1| 0| 1| 4| 9|16|
|3|13| 6| 1|-2|-3|-2| 1| 6|13|
|4|12| 5| 0|-3|-4|-3| 0| 5|12|
|5|13| 6| 1|-2|-3|-2| 1| 6|13|
|6|16| 9| 4| 1| 0| 1| 4| 9|16|
|7|21|14| 9| 6| 5| 6| 9|14|21|
|8|28|21|16|13|12|13|16|21|28|

By identifying zero as the surface threshold(commonly referred to as iso-value), it may be determined that given a continuous function, the surface must exist between two adjacent points if the samples straddle this value. Moreover, each sample may be said to represent the signed integral distance to the surface which may be recentered for isovalues other than zero. By the same property, the inverse may be applied to resolve a signed distance field capable of recreating the surface through Marching-Cubes. 

### Density 

In regards to terrain generation, the scalar field is most widely abstracted to represent density. Though this may be said of almost any application of the algorithm, it is especially appropriate in terrain generation as noise functions are ordinarily normalized in a way that is difficult to visualize as distance. Additionally, abstracting samples to represent density allows for a natural quantification of the terrain beneficial in potential modification of field, as through the equation:

M = D x V

Mass may be reinterpreted to symbolize the quantity of the terrain modified while regular samples allow for a consistent unit of volume to facilitate calculations. Thus with appropriate design, a conservative system essential for closed feedback modification may be developed while remaining intuitive.

### Material 

Employing this density map, it is possible through Marching-Cubes to generate meshes of uniform apperance. To add variation, a simple method is to add superficial details when rendering that are based on factors independent of data accessible to the user such as static noise maps, but this method heavily limits the user's control of the terrain's apperance. A better way is to define a secondary dynamic scalar field accessible to the user that may be considered when determining the apperance of the terrain, or a material field. 

As storing data regarding such details in apperance in the scalar field itself is obviously redundant given a finite number of materials, such variations may be indexed and stored in a look-up table, allowing for the scalar field to store only a reference index. Since Marching-Cubes dictates that a vertex be created on an edge between two samples that straddle the iso-value, it is guaranteed that for every vertex only one of its generating grid-aligned corners is underground. Hence by capturing the material of the underground sample in each vertex, this index may be referenced in rendering to alter the apperance of the terrain. 


### Viscosity 

With these two fields alone, it is enough to approximate any simple terrain by plainly editing the density and material fields alone. However, this is not enough to describe the layered meshes necessary for semi-transparent materials like liquids. 

A naiive solution is to define two density maps for both liquids and the base terrain. Applying Marching-Cubes to both maps would create two independent meshes that could achieve the desired effect. While this solution is simple, careful management is necessary to prevent both maps from overlapping, and even **more** meticulous planning is required to ensure that both meshes connect when close enough. 

A different solution is to define the <i>viscosity</i> of the terrain. Given a density map normalized in an arbitrary positive range <i>min < val < max</i>, an iso-value of exactly half the range guarantees that given a value, either val or max-val is greater than the iso-value. By extension, such an iso-value ensures that if liquid density resolves from the complement of the base density, only one may be underground for a given sample. 

Unfortunately, marking the complement as liquid contradicts the original purpose of the complement which demarcates the above-ground(air). As a consequence to directly replacing density with viscosity, all air is immediately replaced with water. To resolve this, viscosity is kept as a seperate field clamped in the range 0-1 symbolizing the **<i>percentage of density that is solid</i>**. Thus the solid and liquid density may be calculated as follows:
- Solid Density = Viscosity * Base Density
- Liquid Density = (1 - Viscosity) * Base Density
- Air => Solid Density & Liquid Density < (Max Base Density / 2)

Apart from liquid and solid mesh being seperate, this representation has several advantages: 
- The base density can be extracted easily which is expedient for systems that do not require information about its composition
- It allows the combined density to be easily confined in a pre-defined range by clamping Base Density 
- Solid and Liquid Mass Conservation still works

Finally, the problem of connecting both meshes may be solved by executing both March algorithms in parallel, and selectively choosing the solid geometry to keep when two vertices coexist on an edge.

## Considerations

When implementing this form of Marching-Cubes there are several considerations to keep in mind. As since most of these problems have been well established there exist many solutions online that provide insight on addressing them.
- Normals
    - For a smooth surface, it is necessary to resolve the normals based on the adjacent grid information; the resulting normal is simply the change in density along the 3 grid vectors. Information can be found [here](https://stackoverflow.com/questions/27726848/calculating-normals-indices-and-uv-coordinates-for-marching-cubes-mesh)
- Duplicate Vertices
    - When generating mesh geometry on a cube-by-cube basis, adjacent cubes will share edges, and vertices generated on said edges will be duplicated. One may resolve this through a hash-table but the vertex order becomes unpredictable and more importantly, such a solution is lacking for parallel generation. A different solution is to only keep vertices generated next to one corner, and maintain a custom lookup hash to keep track of the index. A later pass can then be arranged to ascertain the correct indexes--Information can be found [here](https://gamedev.net/forums/topic/614060-remove-duplicate-vertices/4878921/)
- Level of Detail
    - For applications where large terrains need to be generated the production & maintenance of all the data immediately may be impossible. This is often resolved using several several levels of details. When different levels of detail meet, there is innevitably a gap formed wherever one grid point fails to find a counterpart in a lower level. This issue may be resolved through several strategies detailed [here](https://transvoxel.org/Lengyel-VoxelTerrain.pdf)

## Code Source
<i>Optimized HLSL Parallel Implementation</i>

{% codeblock lang:C#%}

#pragma kernel March

#include "Assets/Resources/MapData/CCoordHash.hlsl"
#include "Assets/Resources/Utility/GetIndex.hlsl"
#include "Assets/Resources/TerrainGeneration/BaseGeneration/MarchingTables.hlsl"
#include "Assets/Resources/TerrainGeneration/BaseGeneration/MapNormalSampler.hlsl"
#define SAMPLE_TERRAIN 0
#define SAMPLE_WATER 1

static const float3 Epsilon = float3(1E-6, 1E-6, 1E-6);

struct vertex{
    float3 tri;
    float3 norm;
    int2 material; //x = base, y = water
};

struct CubeCorner{
    float solid;
    float liquid;
    int material;
};

struct MarchCube{
    CubeCorner corners[8];
};

struct pInfo{
    uint axis[3];
};


const static int numThreads = 8;



/*               
*  _____________  8  ------
* |  _________  | 7 6      |<-NumOfPoints, Density Map Size
* | |_|_|_|_|_| | 6 5 5 -  |
* | |_|_|_|_|_| | 5 4 4  | |
* | |_|_|_|_|_| | 4 3 3  |<+--ChunkSize / meshSkipInc
* | |_|_|_|_|_| | 3 2 2  | |  "Real Vertices"
* | |_|_|_|_|_| | 2 1 1 _| |
* |_____________| 1_^_^____|
*Purpose:           ^ | |
*   Border Normals⅃ | |
*          Last Edge⅃ |     
*      Number of Cubes⅃
*
*    7---6----6
*   7|       5|
*  / 11      /10
* 4--+-4---5  |
* |  |     9  |
* 8  3---2-+--2
* | 3      | 1
* 0----0---1/
* 
* z
* ^     .--------.
* |    /|       /|
* |   / |      / |    y
* |  .--+-----.  |   /\
* |  |  |     |  |   /
* |  |  .-----+--.  /
* |  | /      | /  /
* | xyz_______./  /
* +---------> x  /
*/

//Map Info
#ifndef MAP_SAMPLER
#define MAP_SAMPLER
StructuredBuffer<uint> _MemoryBuffer;
StructuredBuffer<uint2> _AddressDict;
int3 CCoord;
float meshSkipInc;
int numCubesPerAxis;

const static int POINT_STRIDE_4BYTE = 1;
#endif

float IsoLevel;

//Output
RWStructuredBuffer<uint> counter;
uint3 counterInd; //x -> vertex, y -> baseTri, z -> waterTri

RWStructuredBuffer<vertex> vertexes;

RWStructuredBuffer<pInfo> triangles;
RWStructuredBuffer<uint> triangleDict;
uint bSTART_dict;
uint bSTART_verts;
uint bSTART_baseT;
uint bSTART_waterT;



MarchCube ReadCube(int3 coord){
    MarchCube cube;
    [unroll]for(int i = 0; i < 8; i++){
        int3 sCoord = clamp(coord + positionFromIndex[i], 0, numCubesPerAxis);
        uint mapData = ReadMapData(sCoord, CCoord);
        float density = (mapData & 0xFF) / 255.0f;
        float viscosity = ((mapData >> 8) & 0xFF) / 255.0f;

        cube.corners[i].solid = density * viscosity;
        cube.corners[i].liquid = density * (1-viscosity);
        cube.corners[i].material = ((mapData >> 16) & 0x7FFF);
    }
    return cube;
}


float interpolate(float p1Val, float p2Val) {
    //If p1Val != p2Val, Epsilon is lost through float conversion, otherwise Epsilon prevents Nan
    return ((IsoLevel - p1Val) / (p2Val - p1Val + Epsilon.x)); 
}


[numthreads(numThreads,numThreads,numThreads)]
void March (uint3 id : SV_DispatchThreadID)
{
    if (id.x >= numPointsPerAxis || id.y >= numPointsPerAxis || id.z >= numPointsPerAxis)
        return;

    int3 oCoord = int3(id.xyz);
    int oIndex = indexFromCoord(id.xyz);
    MarchCube cube = ReadCube(oCoord);

    uint baseIndex = 0; uint waterIndex = 0;
    [unroll] for(uint u = 0; u < 8; u++){
        if ((cube.corners[u].solid) < IsoLevel) baseIndex |= (1 << u);
        if ((cube.corners[u].liquid) < IsoLevel) waterIndex |= (1 << u);
    }
    
    uint i;
    //Generate Water first so terrain can override its geometry
    for (i = 0; triangulation[waterIndex][i] != -1; i +=3) {
        [unroll] for(uint v = 0; v < 3; v++){
            //Get the point indexes
            int2 pInd = cornerIndexFromEdge[triangulation[waterIndex][i + v]].xy;

            if(pInd.x != 0) continue; //Duplicate Vertex

            //Figure out point positions
            int3 p1 = oCoord + positionFromIndex[pInd.x];
            int3 p2 = oCoord + positionFromIndex[pInd.y];

            if(any((uint)p2 >= numPointsPerAxis)) continue; //Out of Bounds

            vertex newVert;
            //Get point densities
            float p1Val = cube.corners[pInd.x].liquid;
            float p2Val = cube.corners[pInd.y].liquid;
            float interpFactor = interpolate(p1Val, p2Val);

            newVert.norm = normalize(GetVertexNormal(p1, p2, interpFactor, SAMPLE_WATER));
            newVert.tri = (p1 + interpFactor * (p2 - p1)) * meshSkipInc;
            newVert.material = uint2(cube.corners[p1Val > p2Val ? pInd.y: pInd.x].material, 
                                     cube.corners[p1Val > p2Val ? pInd.x: pInd.y].material);

            //Append Vertex
            int appendInd = 0;
            InterlockedAdd(counter[counterInd.x], 1, appendInd);
            vertexes[bSTART_verts + appendInd] = newVert;

            //Set To Dictionary
            triangleDict[3 * oIndex + cornerToDict[pInd.y] + bSTART_dict] = appendInd;
        }
    }

    for (i = 0; triangulation[baseIndex][i] != -1; i +=3) {
        [unroll] for(uint v = 0; v < 3; v++){
            //Get the point indexes
            int2 pInd = cornerIndexFromEdge[triangulation[baseIndex][i + v]].xy;

            if(pInd.x != 0) continue; //Duplicate Vertex

            //Figure out point positions
            int3 p1 = oCoord + positionFromIndex[pInd.x];
            int3 p2 = oCoord + positionFromIndex[pInd.y];

            if(any((uint)p2 >= numPointsPerAxis)) continue; //Out of Bounds

            vertex newVert;
            //Get point densities
            float p1Val = cube.corners[pInd.x].solid;
            float p2Val = cube.corners[pInd.y].solid;
            float interpFactor = interpolate(p1Val, p2Val);

            newVert.norm = normalize(GetVertexNormal(p1, p2, interpFactor, SAMPLE_TERRAIN));
            newVert.tri = (p1 + interpFactor * (p2 - p1)) * meshSkipInc;
            //If it borders liquid, the liquid will always be in y channel
            newVert.material = uint2(cube.corners[p1Val > p2Val ? pInd.x: pInd.y].material, 
                                     cube.corners[p1Val > p2Val ? pInd.y: pInd.x].material);
            
            //Append Vertex
            int appendInd = 0;
            InterlockedAdd(counter[counterInd.x], 1, appendInd);
            vertexes[bSTART_verts + appendInd] = newVert;

            //Set To Dictionary
            triangleDict[3 * oIndex + cornerToDict[pInd.y] + bSTART_dict] = appendInd;
        }
    }

    //Only add triangles if a valid cube
    if(id.x >= numCubesPerAxis || id.y >= numCubesPerAxis || id.z >= numCubesPerAxis) 
        return;

    for (i = 0; triangulation[baseIndex][i] != -1; i +=3) {
        pInfo baseTri;
        [unroll] for(uint v = 0; v < 3; v++){ //Get the dict index locations of the vertices
            int3 cornerInfo = cornerIndexFromEdge[triangulation[baseIndex][i + v]];
            baseTri.axis[v] = indexFromCoord(oCoord + positionFromIndex[cornerInfo.x]) * 3 + cornerInfo.z;
        }

        //Append Base Triangle
        int appendInd = 0;
        InterlockedAdd(counter[counterInd.y], 1, appendInd);
        triangles[bSTART_baseT + appendInd] = baseTri;
    }

    for (i = 0; triangulation[waterIndex][i] != -1; i +=3) {
        pInfo waterTri;
        [unroll] for(uint v = 0; v < 3; v++){ //Get the dict index locations of the vertices
            int3 cornerInfo = cornerIndexFromEdge[triangulation[waterIndex][i + v]];
            waterTri.axis[v] = indexFromCoord(oCoord + positionFromIndex[cornerInfo.x]) * 3 + cornerInfo.z;
        }

        //Append Base Triangle
        int appendInd = 0;
        InterlockedAdd(counter[counterInd.z], 1, appendInd);
        triangles[bSTART_waterT + appendInd] = waterTri;
    }
}

{% endcodeblock %}
