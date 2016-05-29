#!/bin/env python

# Author: Matt Ralston
# Date: 5/27/16
# Description:
# This is an implementations of the optics algorithm

####################
# PACKAGES
####################
import numpy as np
from scipy.spatial import distance as dist
import operator
from itertools import *
#import sklearn as sk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv
import magic
import os
import argparse
import sys
import ConfigParser
import logging
import logging.config
import collections
import json
#import pdb

####################
# CONSTANTS
####################
min_neighborhood_size = 2
min_maxima_ratio = 0.001


####################
# CLASS DEFS
####################

class TreeNode(object):
    def __init__(self, points, start, end, parentNode):
        self.points = points
        self.start = start
        self.end = end
        self.parentNode = parentNode
        self.children = []
        self.splitpoint = -1

    def __str__(self):
        return "start: %d, end %d, split: %d" % (self.start, self.end, self.splitpoint)

    def assignSplitPoint(self,splitpoint):
        self.splitpoint = splitpoint

    def addChild(self, child):
        self.children.append(child)


####################
# FUNCTIONS
####################

def main():
    X, rownames = read_data(args.input, args.delim)
    print rownames
    RD, CD, order = optics(X, args.min_pts, args.dist_method)#Reachability distances, core distances
    RPlot = [RD[item] for item in order]
    RPoints = [X[item] for item in order]
    
    rootNode = cluster_extraction(RPlot, RPoints, args.min_reach, args.maxima_ratio, args.sample, args.min_pts)
    leaves = getLeaves(rootNode, [])
    graphTree(rootNode, RPlot)
    root_logger.info("Leaves:\n{0}".format(map(lambda x: str(x), leaves)))
    clusters = clusters_as_list(leaves, RPoints)
    tabular_clustering(clusters, X, rownames, args.outdir)
    data = map(lambda x: {"data" : x}, clusters)
    data = internal_indices(data, args.dist_method)
    tabular_metrics(data, args.outdir)
    print_json(data, args.outdir)

              #################################################
              # INTERNAL INDICES
              #################################################
def internal_indices(data, dist_method):
    metric = getattr(dist, dist_method)
    root_logger.info("Calculating centroids...")
    centroids = map(lambda x: np.mean(x["data"], axis=0), data)
    root_logger.info("Calculating inter and intra-cluster dissimilarities...")
    intercluster = [] # m x m  Numpy array of inter-cluster dissimilarities
    intracluster = [] # Numpy array of length m, where each element is an array of intra-cluster dissimilarities whose length is equal to the points in the cluster
    for i in range(len(data)): # For each cluster i
        # Calculate inter-cluster dissimilarity
        intercluster.append(np.zeros(len(data)))
        for j in range(len(data)):
            intercluster[i][j] = metric(centroids[i], centroids[j]) # Calculate the distance between centroids i and j
        root_logger.info("Centroid {0}: {1}".format(i, centroids[i]))
        # Calculate intra-cluster dissimilarity
        dissims = []
        for x in data[i]["data"]:
            dissims.append(metric(x, centroids[i]))
        intracluster.append(dissims)
    avg_intra_dissims = map(lambda x: np.mean(x, axis=0), intracluster) # Numpy array of length m of average intracluster dissimilarities
    root_logger.info("Calculating clustering metrics...")
    for i in range(len(data)):
        # C A L C U L A T E   M E T R I C S
        # intracluster dissimilarity is the average of intracluster dissimilarities
        intra = avg_intra_dissims[i]# np.mean(intracluster[i], axis=0)
        inters = intercluster[i][np.nonzero(intercluster[i])]
        closest = np.where(intercluster[i] == np.min(inters))[0][0]
        min_inter_dissim = intercluster[i][closest]
        furthest = np.where(intercluster[i] == np.max(intercluster[i]))[0][0]
        max_inter_dissim = intercluster[i][furthest]
        silhouette = (min_inter_dissim - intra) / np.max([min_inter_dissim, intra])
        davies_bouldin = (intra + avg_intra_dissims[furthest]) / intercluster[i][furthest]
        dunn = intercluster[i][closest] / np.max(avg_intra_dissims, axis=0)

        data[i]["metrics"] = {"silhouette": silhouette,
                        "davies_bouldin": davies_bouldin,
                        "dunn": dunn,
                        "intracluster_dissimilarity": intra,
                        "min_intercluster_dissimilarity": min_inter_dissim,
                        "max_intercluster_dissimilarity": max_inter_dissim
        }
    return data
        

    
              #################################################
              # CORE OPTICS ALGORITHM
              #################################################
def cluster_extraction(RPlot, RPoints, min_reach, maxima_ratio, sample, minpts):
    neighborhood_size = int(min_maxima_ratio * len(RPoints))
    if neighborhood_size < min_neighborhood_size: neighborhood_size = min_neighborhood_size
    localMaximaPoints = findLocalMaxima(RPlot, RPoints, neighborhood_size)
    root_logger.debug("Initial maxima: {0}".format(localMaximaPoints))
    rootNode = TreeNode(RPoints, 0, len(RPoints), None)
    clusterTree(rootNode, None, localMaximaPoints, RPlot, RPoints, min_reach, maxima_ratio, sample, minpts)
    return rootNode
    
def clusterTree(node, parentNode, localMaximaPoints, RPlot, RPoints, min_reach, maxima_ratio, sample, minpts):
    # 'node' is a TreeNode object. It is either the root with no parent or a child of another node
    # 'parentNode' is the parend node of 'node' or 'None' if node is the root node
    if len(localMaximaPoints) == 0:
        root_logger.debug("No more maxima points")
        return

    # take larget local maximum as possible separation
    s = localMaximaPoints[0]
    node.assignSplitPoint(s)
    root_logger.debug("Assigning split point/new node: {0}".format(node))
    localMaximaPoints = localMaximaPoints[1:]

    # create two new nodes and add to list of nodes
    Node1 = TreeNode(RPoints[node.start:s], node.start, s, node)
    Node2 = TreeNode(RPoints[s+1:node.end], s+1, node.end, node)
    root_logger.debug("Two child nodes:\nLeft: {0}\nRight: {1}".format(Node1, Node2))
    LocalMax1 = [i for i in localMaximaPoints if i < s]
    LocalMax2 = [i for i in localMaximaPoints if i > s]
    Nodelist = []
    Nodelist.append((Node1, LocalMax1))
    Nodelist.append((Node2, LocalMax2))
    if RPlot[s] < min_reach:
        root_logger.debug("Split point {0} is not significant: {1} < {2}".format(node, RPlot[s], min_reach))
        node.assignSplitPoint(-1)
        # if splitpoint is not significant, ignore this split and continue
        clusterTree(node, parentNode, localMaximaPoints, RPlot, RPoints, min_reach, maxima_ratio, sample, minpts)
        return
    # only check a certain ratio of points in the child nodes formed to the left and right of the maxima
    # Calculate the reachability distance for a region
    checkValue1 = int(np.round(sample*len(Node1.points)))
    checkValue2 = int(np.round(sample*len(Node2.points)))
    if checkValue2 == 0: checkValue2 = 1
    avgReachValue1 = float(np.average(RPlot[(Node1.end - checkValue1):Node1.end]))
    avgReachValue2 = float(np.average(RPlot[Node2.start:(Node2.start + checkValue2)]))

    # Use all points
    #avgReachValue1 = float(np.average(RPlot[Node1.start - Node1.end]))
    #avgReachValue2 = float(np.average(RPlot[Node2.start - Node2.end]))

    ##########################
    # REMOVE NODES THAT DONT PASS LOCAL MAXIMA TEST
    ##########################
    # low ratio = more stringency for local maximum (more difference between cluster)
    leftFail = float(avgReachValue1) / float(RPlot[s]) >= maxima_ratio
    if leftFail:
        root_logger.debug("Points on the left of the split point do not satisfy the maxima ratio")
        Nodelist.remove((Node1, LocalMax1))
    rightFail = float(avgReachValue2) / float(RPlot[s]) >= maxima_ratio
    if rightFail:
        root_logger.debug("Points on the right of the split point do not satisfy the maxima ratio")
        Nodelist.remove((Node2, LocalMax2))
    if leftFail and rightFail:
        node.assignSplitPoint(-1)
        clusterTree(node, parentNode, localMaximaPoints, RPlot, RPoints, min_reach, maxima_ratio, sample, minpts)
        return
    ##########################
    # REMOVE SMALL NODES
    ##########################
    # Retroactively remote clusters that are too small
    if len(Node1.points) < minpts:
        try:
            Nodelist.remove((Node1, LocalMax1))
        except Exception as e:
            root_logger.error("Nodes:\n{0}".format(Nodelist))
            root_logger.error("Failed to remove node! {0}".format(Node1))
            root_logger.error(e)
            sys.exc_clear()
    if len(Node2.points) < minpts:
        try:
            Nodelist.remove((Node2, LocalMax2))
        except Exception as e:
            root_logger.error("Nodes:\n{0}".format(map(lambda x: str(x), Nodelist)))
            root_logger.error("Failed to remove node! {0}".format(Node2))
            root_logger.error(e)
            sys.exc_clear()
    if len(Nodelist) == 0:
        node.assignSplitPoint(-1)
        return
    ##########################
    # MERGE NODES
    ##########################
    bypassNode = False
    if parentNode is not None:
        parentSD = np.std(RPlot[parentNode.start:parentNode.end])
        parentChildDist = float(np.average(RPlot[node.start:node.end]) - np.average(RPlot[parentNode.start:parentNode.end]))
        if parentChildDist < parentSD:
            parentNode.children.remove(node)
            bypassNode = True
    for n in Nodelist:
        if bypassNode:
            parentNode.addChild(n[0])
            clusterTree(n[0], parentNode, n[1], RPlot, RPoints, min_reach, maxima_ratio, sample, minpts)
        else:
            node.addChild(n[0])
            clusterTree(n[0], node, n[1], RPlot, RPoints, min_reach, maxima_ratio, sample, minpts)

# Find all maxima and return the index of the points sorted in terms of their reachability distance
def findLocalMaxima(RPlot, RPoints, neighborhood_size):
    localMaximaPoints = {}
    # 1st and last points cannot be maxima
    for i in range(1, len(RPoints)-1):
        if RPlot[i] > RPlot[i-1] and RPlot[i] >= RPlot[i+1] and isLocalMaxima(i, RPlot, RPoints, neighborhood_size):
            localMaximaPoints[i] = RPlot[i]
    return sorted(localMaximaPoints, key=localMaximaPoints.__getitem__ , reverse=True)

def isLocalMaxima(index, RPlot, RPoints, neighborhood_size):
    for i in range(1, neighborhood_size + 1): #Determine if the point has the maximum reachability distance from its neighbors with respect to the neighborhood size in both directions
        if index + i < len(RPlot) and (RPlot[index] < RPlot[index+1]):
            return False
        if index - i >= 0 and (RPlot[index] < RPlot[index-i]):
            return False
    return True
    
# X is a numpy array
def optics(X, min_pts, dist_method):
    m,n = X.shape
    root_logger.info("Calculating distance matrix...")
    D = dist.squareform(dist.pdist(X, dist_method)) # Calculate distance matrix
    CD = np.zeros(m) # Initialize the core distance
    RD = np.ones(m)*1E10 # Initialize the reachability distances

    root_logger.info("Initializing core distance as the distance to each points {0} nearest neighbor...".format(grammatical_number(min_pts)))
    for i in xrange(m): # Initialize core distance for point i as the min_pts'th nearest neighbor
        tempInd = D[i].argsort()
        tempD = D[i][tempInd]
        CD[i] = tempD[min_pts]

    order = []
    seeds = np.arange(m, dtype=np.int)
    root_logger.info("Arranging points by their reachability distance, the max of core distance or distance to the nearest neighbor...")
    ind = 0
    while len(seeds) != 1:
        ob = seeds[ind]
        seedInd = np.where(seeds != ob)
        seeds = seeds[seedInd]

        order.append(ob)
        tempX = np.ones(len(seeds))*CD[ob] # Create a vector of len(seeds) of the core distance of the object
        tempD = D[ob][seeds] # Create a vector of distances for the ob'th point/row with the columns specified by seed

        temp = np.column_stack((tempX, tempD))
        mm = np.max(temp, axis=1) # Create a vector of maximum distances among the core distance and the distance between points
        ii = np.where(RD[seeds] > mm)[0] # Find the indices of the reachability distance array where the current value (initially huge) are larer than the values of the maximum vector 'mm'
        RD[seeds[ii]] = mm[ii] # Replace these values with the maxima
        ind = np.argmin(RD[seeds]) # Find the closest point and then consider its distances to its neighbors...
    order.append(seeds[0])
    RD[0] = 0 # The zeroth point has 0 distance to itself
    return (RD, CD, order)

              #################################################
              # ACCESSORY FUNCTIONS
              #################################################
def printTree(node, num):
    if node is not None:
        print "Level %d" % num
        print str(node)
        for n in node.children:
            printTree(n, num+1)

def writeTree(fileW, locationMap, RPoints, node, num):
    if node is not None:
        fileW.write("Level " + str(num) + "\n")
        fileW.write(str(node) + "\n")
        for x in range(node.start,node.end):
            item = RPoints[x]
            lon = item[0]
            lat = item[1]
            placeName = locationMap[(lon,lat)]
            s = str(x) + ',' + placeName + ', ' + str(lat) + ', ' + str(lon) + '\n'
            fileW.write(s)
        fileW.write("\n")
        for n in node.children:
            writeTree(fileW, locationMap, RPoints, n, num+1)

def grammatical_number(x):
    if x is 1:#st
        return "1st"
    elif x is 2:#nd
        return "2nd"
    elif x is 3:#rd
        return "3rd"
    else:
        return "{0}th".format(x)
            
def getArray(node,num, arr):
    if node is not None:
        if len(arr) <= num:
            arr.append([])
        try:
            arr[num].append(node)
        except:
            arr[num] = []
            arr[num].append(node)
        for n in node.children:
            getArray(n,num+1,arr)
        return arr
    else:
        return arr

def getLeaves(node, arr):
    if node is not None:
        if node.splitpoint == -1:
            arr.append(node)
        for n in node.children:
            getLeaves(n,arr)
    return arr

def graphTree(root, RPlot):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    a1 = [i for i in range(len(RPlot))]
    ax.vlines(a1, 0, RPlot)
    num = 2
    graphNode(root, num, ax)
    plt.savefig('RPlot.png', dpi=None, facecolor='w', edgecolor='w',
      orientation='portrait', papertype=None, format=None,
     transparent=False, bbox_inches=None, pad_inches=0.1)
    #plt.show()
            
def graphNode(node, num, ax):
    ax.hlines(num,node.start,node.end,color="red")
    for item in node.children:
        graphNode(item, num - .4, ax)

def clusters_as_list(leaves, RPoints):
    return map(lambda x: RPoints[x.start:x.end], leaves)

def tabular_clustering(clusters, rawdata, rownames, outdir):
    outfile = os.path.join(outdir, "clustering.txt")
    with open(outfile, 'w') as ofile:
        ofile.write("id\tcluster\n")
        for c in range(len(clusters)):
            for point in clusters[c]:
                # Ensure that there is only one index for this point
                indices = [i for i,x in enumerate(rawdata) if np.array_equal(x, point)] 
                if len(indices) is not 1:
                    root_logger.error("The following point occurs multiple times in this dataset: {0}".format(point))
                ofile.write(str(rownames[indices[0]]) + "\t" + str(c) + "\n")

def tabular_metrics(metrics, outdir):
    outfile = os.path.join(outdir, "metrics.txt")
    keys = metrics[0]["metrics"].keys()
    with open(outfile, 'w') as ofile:
        ofile.write("cluster\t" + "\t".join(keys) + "\n")
        for c in range(len(metrics)):
            values = map(lambda x: str(metrics[c]["metrics"][x]), keys)            
            ofile.write(str(c) + "\t" + "\t".join(values) + "\n")

def print_json(data, outdir):
    outfile = os.path.join(outdir, "optics.json")
    with open(outfile, 'w') as ofile:
        ofile.write(json.dumps(dictionary,indent=2))

def read_data(infile, delimiter):
    data = []
    with open(infile,'r') as csvfile:
        for row in csv.reader(csvfile, delimiter=delimiter):
            if len(row) > 0: data.append(row)
    header = False if will_it_float(data[0][1]) else True # If there is a 'float'able number in the 0,1 position of the matrix then there is no header
    rownames = False if will_it_float(data[1][0]) else True # If there is a 'float'able number in the 1,0 position of the matrix then there are no rownames

    if header: data = data[1:] 
    if rownames:
        rownames = map(lambda x: x[0], data) # Assign rownames to the 
        data = map(lambda x: x[1:], data) # Remove the rownames
    else:
        rownames = range(0,len(data))
    data = np.array(map(lambda x: map(lambda y: float(y), x), data))
    return (data, rownames)
    
def will_it_float(maybe_num):
    try:
        float(maybe_num)
        return True
    except ValueError:
        return False


def get_root_logger(loglevel):
   # Requires 'import logging' and 'import logging.config'
    def log_level(loglevel):
        case = {"DEBUG": logging.DEBUG,
                "INFO": logging.INFO,
                "WARNING": logging.WARNING,
                "ERROR": logging.ERROR}
        return case[loglevel.upper()]
    logging.basicConfig(level=log_level(loglevel),format="%(levelname)s: %(asctime)s %(funcName)s L%(lineno)s| %(message)s",datefmt="%Y/%m/%d %I:%M:%S %p")
    root_logger = logging.getLogger()
    log_format = root_logger.handlers[0].format
    return root_logger

def validate_args():
    invalid = False
    if not (os.path.isfile(args.input) and magic.from_file(args.input) == "ASCII text"):
        root_logger.error("Argument --input should be an existing, ASCII format matrix.")
        invalid = True
    if not os.path.isdir(args.outdir):
        root_logger.error("Argument --outdir should be an existing directory.")
        invalid = True
    try:
        args.min_reach = float(args.min_reach)
    except ValueError as e:
        root_logger.fatal(e)
        root_logger.error("Argument --min-reach could not be converted to floating point number...")
        invalid = True
    try:
        args.maxima_ratio = float(args.maxima_ratio)
        if args.maxima_ratio <= 0 or args.maxima_ratio >= 1:
            root_logger.error("Argument --maxima-ratio must be between 0 and 1...")
            invalid = True
    except ValueError as e:
        root_logger.fatal(e)
        root_logger.error("Argument --maxima-ratio cannot be converted to floating point number...")
        invalid = True
    try:
        args.min_pts = int(args.min_pts)
        if args.min_pts <= 2:
            root_logger.error("Argument --min-pts must be greater than 2")
            invalid = True
    except ValueError as e:
        root_logger.fatal(e)
        root_logger.error("Argument --min-pts could not be converted to an integer...")
        invalid = True
    try:
        args.sample = float(args.sample)
        if args.sample <= 0 or args.sample >= 1:
            root_logger.error("Argument --sample must be between 0 and 1...")
            invalid = True
    except ValueError as e:
        root_logger.fatal(e)
        root_logger.error("Argument --sample cannot be converted to floating point number...")
        invalid = True
    if invalid: sys.exit(2)
    root_logger.debug(args)



####################
# OPTIONS AND MAIN
####################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help="Data in delimited format",required=True)
    parser.add_argument('--outdir', help="Directory of output", default=os.getcwd())
    parser.add_argument('--min-reach', help="Minimum reachability distance. Beware the curse of dimensionality...",required=True)
    parser.add_argument('--dist-method', help="Distance method from NumPy",choices=["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "dice", "euclidean", "hamming", "jaccard", "julsinski", "matching", "rogerstanimoto", "russellrao", "sokalmichener", "sokalsneath", "sqeuclidean", "yule",],default="euclidean")
    parser.add_argument('--min-pts', help="Minimum number of points per cluster",default=3)
    parser.add_argument('--maxima-ratio', help="Maxima ratio: inter-cluster reachability distance / average intra-cluster reachability distance",default=0.5)
    parser.add_argument('--sample', help="Fraction of putative cluster to sample to calculate the average intra-cluster reachability distance",default=0.75)
    parser.add_argument('--delim', help="Delimiter of data matrix. Default is tab '\t'", default="\t")
    parser.add_argument('--log-level', help="Prints warnings to console by default",default="WARNING",choices=["DEBUG","INFO","WARNING","ERROR"])
    args = parser.parse_args()
    # Set up the root logger
    root_logger=get_root_logger(args.log_level)
    validate_args()
    # Main routine
    main()
