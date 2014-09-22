# README
This repository is an implementation of the optics algorithm, combined with the automated cluster detection outlined in:
http://dx.doi.org/10.1007/3-540-36175-8_8

The code for the hierarchical clustering comes from AmyXZhang, who authored the original project.
Modifications include the original standard deviation cluster level criteria (see article), user defined parameters, and output of the cluster identity to .csv.

## Inputs
### Required
1. A csv file of the data to be clustered: no column names and arbitrary identifiers for the rows.
2. A destination directory for data to be written

### Optional
3. Minimum cluster size [0.5% of the data set or 3, whichever is larger]
4. Minimum reachability distance for maxima [0.003]
5. Maxima ratio (Maximum allowed ratio of the maxima's reachability distance / average reachability distance of bordering clusters) [0.75]
6. Distance method, see http://hcluster.damianeads.com/cluster.html

## Outputs
1. optics-clustering.csv -- a csv list of the data points and their respective cluster
2. Reachability_plot.png -- an alternative to a hierarchical clustering tree

## Dependencies
This python script requires numpy, matplotlib, and hcluster.
All are available on PyPi.
