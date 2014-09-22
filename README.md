# README
This repository is an implementation of the optics algorithm, combined with the automated cluster detection outlined in:
http://dx.doi.org/10.1007/3-540-36175-8_8

The code for the hierarchical clustering comes from AmyXZhang, who authored the original project.

## Inputs
### Required
1. A csv file of the data to be clustered
2. A destination directory for data to be written

### Optional
3. Minimum cluster size (defaults to 0.5% of the data set or 3 whichever is larger)
4. Minimum reachability distance for maxima (significantMin)
5. Maxima ratio (Maximum allowed ratio of the maxima's reachability distance / average reachability distance of bordering clusters)
6. Distance method, see http://hcluster.damianeads.com/cluster.html


## Dependencies
This python script requires numpy, matplotlib, and hcluster.
All are available on PyPi.
