# Infrared Guide Star Catalog (IRGSC)

## Introduction

This is a repository dedicated to generating the IRGSC for the A0 observations of the Thirty Meter Telescope. The code computes the expected NIR magnitudes of the sources in the PS1 optical data in g, r, i, z, and y filters.

### Installation
```
pip install tmt-IRGSC
```

# Usage
This python package is made to generate a catalog of the NIR magnitudes of the stellar sources in the sky. The class irgsc requires four inputes: R.A. (degrees, float), Decl. (degrees, float) , search radius (arcminutes, float), validate (bool). For a given field, the code performs a cone search around the coordinates with the given search radius and obtains the optical photometry and information flags of the sources from PANSTARRS DR2 (Chambers et.al. 2016) obtains the Gaia DR3 (Gaia collaboration et.al. 2022) for these sources. The code then seperates the stars and galaxies in the data, models the optical colors of the stellar sources using the Kurucz (Kurucz 1992a,b, 1993), CAstelli-Kurucz (Castelli & Kurucz 2003) and Phoenix/NextGen (Hauschildt 1999a,b) and computes the NIR magnitudes for them.

```
from tmt-irgsc import irgsc
ic = irgsc(ra = None, dec = None, search_radius = None, validate=False)
```

If validate is True then the code also obtains the observed NIR UKIDSS data for the given field (default is False). This data is then used to compare the observed and computed NIR magnitudes for the source.


# Acknowledgements
Please add the following acknowledgment if you use our package in your work.

"This work has made use of Infrared Guide Star Catalog (IRGSC) developed as a part of the Thirty Meter Telescope (TMT) project."

If you have any questions or suggestions for improvements to this repo, please contact the owners of the repository.


## References
1. Hauschildt, P. H., Allard, F., & Baron, E. 1999a, apj, 512, 377, doi: 10.1086/306745
2. Hauschildt, P. H., Allard, F., Ferguson, J., Baron, E., & Alexander, D. R. 1999b, apj, 525, 871,
      doi: 10.1086/307954
3. Chambers, K. C., Magnier, E. A., Metcalfe, N., et al. 2016, arXiv e-prints, arXiv:1612.05560.
      https://arxiv.org/abs/1612.05560
4. Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. 2022, arXiv e-prints, arXiv:2208.00211.
      https://arxiv.org/abs/2208.00211
5. Kurucz, R. L. 1992a, rmxaa, 23, 45
6. Kurucz, R. L. 1992b, in The Stellar Populations of Galaxies, ed. B. Barbuy & A. Renzini, Vol. 149, 225
7. Kurucz, R. L. 1993, in Astronomical Society of the Pacific Conference Series, Vol. 44, IAU Col-
      loq. 138: Peculiar versus Normal Phenomena in A-type and Related Stars, ed. M. M. Dworetsky,
       F. Castelli, & R. Faraggiana, 87
