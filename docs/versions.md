# Version history

To find out what version you're using, run `dreamcoat.hello()`.

## <span style='color:#653700; font-weight: 1000'>Brown</span> (0.4)

Incorporates dependencies [great_circle_calculator](https://github.com/seangrogan/great_circle_calculator) and [vptree](https://github.com/RickardSjogren/vptree) internally to allow adding dreamcoat to conda-forge.

## <span style='color:#15b01a; font-weight: 1000'>Green</span> (0.3)

> Developed for RV *Pelagia* cruises 64PE533 and 64PE534 in April-May 2024.

Simplifies the package to remove messy dependencies and focus on core cruise data processing capabilities.  Removes CMEMS, MODIS and emailing tools.

## <span style='color:#ffff14; font-weight: 500; text-shadow: -1px 1px 0 #ceb301aa, 1px 1px 0 #ceb301aa, 1px -1px 0 #ceb301aa, -1px -1px 0 #ceb301aa;'>Yellow</span> (0.2)

> Developed for RV *Pelagia* cruise 64PE517 (26<sup>th</sup> May to 14<sup>th</sup> June 2023) in the North Sea, a round trip from Texel, the Netherlands to the Norwegian Trench as part of the NoSE research project.

Adds mapping tools for working with different coordinate formats and calculating distance to shore.

## <span style='color:#e50000; font-weight: 1000'>Red</span> (0.1)

> Developed for RV *Pelagia* cruise 64PE513 (10<sup>th</sup> February to 5<sup>th</sup> March 2023) in the Atlantic Ocean, a transect from Cape Town, South Africa to Mindelo, Cape Verde as part of the BEYΩND research project.

Creates a set of tools for downloading, processing and plotting subsets of global surface ocean reanalysis / forecast datasets from CMEMS and NASA Ocean Color and sending the data and/or figures by email.
