[![DOI](https://zenodo.org/badge/951415482.svg)](https://doi.org/10.5281/zenodo.15088479)


This code has been created by Elisa Carli, Lia Siegelman and Patrice Klein.
If you use this code, please cite [the Zenodo repository for the code](https://doi.org/10.5281/zenodo.15088479) and Carli et al. 2024. [10.1029/2024JC021216](10.1029/2024JC021216)


# eSQG
How to apply the eSQG to altimetry surface fields and retrieve vertical velocities 

The effective version of the Surface Quasi Geostrophy theory is first instroduced by Lapeyre and Klein [(2006)](https://journals.ametsoc.org/view/journals/phoc/36/2/jpo2840.1.xml) and is based on the reduction of Eady's model for baroclinic instabilities from one boundary (here the ocean surface), used to describe the oceanic mesoscale eddy fields in the surface layers (Klein & Lapeyre, [2009](https://www.annualreviews.org/content/journals/10.1146/annurev.marine.010908.163704)).
This approach is valid for the surface intensified eddies. The reconstruction of other eddies, namely subsurface intensified eddies, requires the combination of the effective sQG with the interior sQG (iSQG) as in Wang et al. [(2013)](https://journals.ametsoc.org/view/journals/phoc/43/8/jpo-d-12-0204.1.xml), where both have a signature in SSH. 

The major assumption for effective sQG is that potential vorticity ($PV$) is uniform and constant in the ocean interior. Using $PV$  inversion, one can retrieve a geostrophic stream function ($ψ$) related to $PV$ using an elliptic operator (Klein et al.,  [2009](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2009GL038359)). Following Klein et al. [(2009)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2009GL038359) and assuming a doubly periodic domain, the geostrophic stream function is retrieved at all depths ($z$) from SSH ($η$) using an exponential decay as

$$\hat{\psi} (\mathbf{k}, z) = \frac{g}{f}\hat{η}(\mathbf{k})exp\left(\frac{N_0}{f}kz\right)$$

where  $\widehat{\left(.\right)}$   corresponds to the Fourier transform, $\mathbf{k} = (kx, ky)$ is the wavenumber vector and $k = \mathbf{k}$ is its norm, $g$ is the gravity constant, $f$ is the Coriolis frequency, $N_0$ is a Brunt‐Väissälä frequency that accounts for the contribution of the interior $PV$. Vertical motions are retrieved from the buoyancy equations as in Klein et al. [(2009)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2009GL038359)

$$\hat{w}(\mathbf{k}, z) = −\frac{c^2}{N^2_0}\left[−J\left(\widehat{ψ_s, b_s}\right)exp\left(\frac{N_0}{f}kz\right) + J\left(\widehat{ψ, b}\right)\right]$$

where the subscript $s$ refers to the surface values, $b$ is the buoyancy, $c$ is a unitless constant that helps optimize the reconstructed fields' amplitude and $J(A, B) ≡ (∂_xA∂_yB − ∂_yA∂_xB)$ . This methodology requires the optimization of the constants $N_0$ and $c$.

1. $N_0$ optimization: we want to minimize the mean squared error between the model and eSQG vorticity from the bottom of the mixed layer to 1000 m. 
2. $c$ optimization: we aim to align the highest values of the modeled and reconstructed w. So we create a mask to select the strongest $1\sigma$ w in the model. We then apply the same mask to select the strongest reconstructed eSQG w. The $c$ is computed so that the difference between the two selected fields is minimized in average from the bottom of the mixed layer and 1000 m . 
This process is done for every time step of the available time series and one overall average value per season is chosen.

The eSQG methodology is valid for structures between ~20 km and 400 km, from below the mixed layer to ~1000 m. This means that the SSH field, before being used for the reconstruction, needs to be pre-processed and spatially filtered.

This method uses a spectral analysis, so there are a few things to keep in mind for the structure of the domain:
- The domain should be in the open ocean, no costs or islands 
- The x and y dimensions should be similar otherwise structures will have the tendency to align along one preferential direction 
