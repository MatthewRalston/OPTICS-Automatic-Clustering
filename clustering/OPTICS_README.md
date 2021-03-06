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
3. <a href="http://en.wikipedia.org/wiki/OPTICS_algorithm">Minimum cluster size</a> [0.5% of the data set or 3, whichever is larger]
4. <a href="http://en.wikipedia.org/wiki/OPTICS_algorithm">Minimum reachability distance</a> for maxima [0.003]
5. <a href="http://en.wikipedia.org/wiki/OPTICS_algorithm#Extracting_the_clusters">Maxima ratio</a> (Mainimum accepted ratio of the average reachability distance of bordering clusters / maxima's reachability distance ) [0.75]
5. [Clustering distance method](http://docs.scipy.org/doc/scipy/reference/spatial.distance.html), used to calculate the distance matrix for clustering.
6. [Distance method](http://docs.scipy.org/doc/scipy/reference/spatial.distance.html), used to calculate all distances and dissimilarities
7. The <a href="http://en.wikipedia.org/wiki/OPTICS_algorithm#Extracting_the_clusters">percentage</a> of the middle of the cluster to sample to calculate the average reachability distance [0.75]

## Returns
Prints an array of summary statistics, calculated from output below.
* [Avg. Silhouette coefficient](http://en.wikipedia.org/wiki/Silhouette_(clustering))
* [Davies-Bouldin index](http://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index)
* [Dunn index](http://en.wikipedia.org/wiki/Dunn_index)
* Avg. intracluster dissimilarity
* Avg. minimum intercluster distance
* Avg. maximum intercluster distance

## Outputs
* 1\. optics-clustering.csv -- a csv list of the data points and their respective cluster
* 2\. clustering-metrics.csv -- an ordered csv list of metrics, per cluster 'i'
   * 1.1\. Silhouette coefficient = *avg { (b<sub>i</sub> - a<sub>i</sub>) / max(a<sub>i</sub>,b<sub>i</sub>) }*,

		where a is the average intracluster dissimilary, and b is the minimum inter-cluster dissimilarity.
   * 1.2\. Davies-Bouldin 'coefficient' = *(a<sub>i</sub> + a<sub>j</sub>)/d(i,j)*,

		where a<sub>x</sub> is the average intracluster dissimilarity of x, j refers to the farthest cluster from j, and d(i,j) is the distance between their centroids.
   * 1.3\.Dunn 'coefficient' = *min( d(i,j) ) / a<sub>k</sub>*,

		where d(i,j) is intercluster distance and a<sub>k</sub> is the maximum intracluster dissimilarity.
   * 1.4\. Intracluster dissimilarity
   * 1.5\. Minimum intercluster distance
   * 1.6\. Maximum intercluster distance
* 3\. Reachability_plot.png -- an alternative to a hierarchical clustering tree


## Dependencies
This python script requires numpy, ,scipy, sci-kit learn, matplotlib, and hcluster.
All are available on PyPi.
