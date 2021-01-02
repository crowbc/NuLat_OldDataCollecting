# NuLatDec2020

## **To Run:**

1. Download and unpack files

2. Create build directory in 

    `$ mkdir build4`

3. Enter build directory and run cmake (never needs to be done again)

    `$ cd build4`
    
    `$ cmake ..`
    
4. Run make (needs to be run if any files are changed)

    `$ make`
    
5. Run the simulation (starts the simulation)

    `$ ./Nulat`
    
    
    
## **Results:**

* NuLatEvent (General run data including all initial conditions, total count of PMT hits, and total count of Voxel hits)

* ParticlesInVoxels (Every particle to hit a voxel: particle ID, energy deposit, precise position, time)

* PhotonsInPMTs (Every photon to hit a PMT: PMT location, Photon Energy (spectrum from scintillation), time)

* PMTHits (Counts total hits in each PMT)



## **Setting Initial Conditions, etc.**
- **Input Particle (in NuLatPrimaryParticleGenerator.cc)**
  - Choose between Simulated IBD event, Single Test Particle, or Photon Source in GeneratePrimaries Function
    - GenerateIBDevent (generates 2MeV Positron + 400 keV Neutron; position, time, momentum direction adjustable within function)
    - GeneratePhotonEvent (generates 10,000 4 eV Photons with random momentum direction)
    - GenerateTestParticle (generates test particle; particle type is function input; other variables edited within function)
  
- **Output format (in Analysis.hh)**
  - Uncomment prefered output source, make sure other options are commented
    - ROOT
    - .CSV (default)
    - XML

- **Parameters for if events get recorded or forgotten in both SensitiveDetector.cc files**
  - Can filter out useless events such as photon absorbsions in Voxels

- **Detector Geometry (NuLatDetectorConstuction.cc and Materials.hh)**
  - Material Properties (scintillation, density, molecular makeup etc.) in Materials.hh
  - Physical Geometries (Shapes of pieces, what each piece is made of, multipart construction) in NuLatDetectorConstruction.cc

