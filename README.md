# eSQG
Code for solving the effective Surface Quasi Geostrophy method to reconstruct vertical velocities in the open ocean


N.B. This method uses a spectral analysis, so there are a few things to keep in mind for your domain structure:
- The domain should be in the open ocean, no costs or islands 
- The x and y dimensions should be similar in dimension otherwise structures will have the tendency to align along one preferential dimension
- You should have an initial guess of the vertical structure of the water column with a model, to initialize the reconstruction. The two coefficients to optimize are the buoyancy frequency *N2* and the multiplicative coefficient *c* which
- The eSQG methodology is valid for structures between ~20 km and 400 km, from below the mixed layer to ~1000 m. This means that the SSH field, before being used for the reconstruction, needs to be pre-processed.
