"""
@author: Amy X Zhang
amy.xian.zhang@gmail.com
http://amyxzhang.wordpress.com


Modified by Matthew T. Ralston
mrals89@gmail.com
http://matthewralston.github.com

This automatic optics cluster extraction was modified to produce
a list of cluster identities.


"""
import csv, sys, getopt
import numpy as np
import matplotlib.pyplot as plt
from itertools import *
from operator import itemgetter




def main(argv):
    ifile = ''
    odir = ''
    mpts = False
    minReach = False
    maxRatio = False
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:",["input_file=","output_dir=","min_pts=","min_reach=","maxima_ratio="])
    except getopt.GetoptError:
        print 'optics-automated.py -i <input file> -o <output directory> -n <minimum pts> -r <minimum reachability distance> -m <maxima ratio>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'optics-automated.py -i <input file> -o <output directory> -n <minimum pts> -r <minimum reachability distance> -m <maxima ratio>'
            sys.exit()
        elif opt in ("-i","--input_file"):
            ifile = arg
        elif opt in ("-o","--output_dir"):
            odir = arg
        elif opt in ("-n","--min_pts"):
            mpts = checkint(arg)
        elif opt in ("-r","--min_reach"):
            minReach = checkfloat(arg)
        elif opt in ("-m","--maxima_ratio"):
            maxRatio = checkfloat(arg)
        else:
            print 'optics-automated.py -i <input file> -o <output directory> -n <minimum pts> -r <minimum reachability distance> -m <maxima ratio>'
            sys.exit(2)
    autoclust(ifile,odir,mpts,minReach,maxRatio)

def checkint(arg):
    try:
        return int(arg)
    except ValueError:
        print "-n should be an integer"
        sys.exit(2)

def checkfloat(arg):
    try:
        return float(arg)
    except ValueError:
        print "-r and -m should be a floating point number"
        sys.exit(2)

def autoclust(infile,outdir,minpts,minReach,maxRatio):
    X=[]
    # Load some data, a NxM dataset of N objects, M variables
    with open(infile,'r') as csvfile:
        data = csv.reader(csvfile,delimiter=',')
        for row in data:
            X+=row
            
    m,n = X.shape


    #run the OPTICS algorithm on the points, using a smoothing value (0 = no smoothing)
    RD, CD, order = OP.optics(X,9)

    RPlot = []
    RPoints = []
        
    for item in order:
        RPlot.append(RD[item]) #Reachability Plot
        RPoints.append([X[item]]) #points in their order determined by OPTICS

    #hierarchically cluster the data
    if  maxRatio and minReach:
        rootNode = automaticCluster(RPlot, RPoints, minpts, minReach=minReach, maximaRatio=maxRatio)
    elif maxRatio:
        rootNode = automaticCluster(RPlot, RPoints, minpts, maximaRatio=maxRatio)
    elif minReach:
        rootNode = automaticCluster(RPlot, RPoints, minpts, minReach=minReach)
    else:
        rootNode = automaticCluster(RPlot, RPoints, minpts)

    #print Tree (DFS)
    printTree(rootNode, 0)

    #graph reachability plot and tree
    graphTree(rootNode, RPlot)

    #array of the TreeNode objects, position in the array is the TreeNode's level in the tree
    array = getArray(rootNode, 0, [0])

    #get only the leaves of the tree
    leaves = getLeaves(rootNode, [])


    clusters=cycle(''.join(str(c) for c in range(1,len(leaves)+1)))
    final_list = []
    for clust, c in zip(leaves,clusters):
        for x in range(clust.start,clust.end):
            final_list.append(RPoints[x]+[c])


# _________ A U T O C L U S T E R

def isLocalMaxima(index,RPlot,RPoints,nghsize):
    # 0 = point at index is not local maxima
    # 1 = point at index is local maxima
    
    for i in range(1,nghsize+1):
        #process objects to the right of index 
        if index + i < len(RPlot):
            if (RPlot[index] < RPlot[index+i]):
                return 0
            
        #process objects to the left of index 
        if index - i >= 0:
            if (RPlot[index] < RPlot[index-i]):
                return 0
    
    return 1

def findLocalMaxima(RPlot, RPoints, nghsize):
    
    localMaximaPoints = {}
    
    #1st and last points on Reachability Plot are not taken as local maxima points
    for i in range(1,len(RPoints)-1):
        #if the point is a local maxima on the reachability plot with 
        #regard to nghsize, insert it into priority queue and maxima list
        if RPlot[i] > RPlot[i-1] and RPlot[i] >= RPlot[i+1] and isLocalMaxima(i,RPlot,RPoints,nghsize) == 1:
            localMaximaPoints[i] = RPlot[i]
    
    return sorted(localMaximaPoints, key=localMaximaPoints.__getitem__ , reverse=True)
    


def clusterTree(node, parentNode, localMaximaPoints, RPlot, RPoints, min_cluster_size, minReach, maximaRatio):
    #node is a node or the root of the tree in the first call
    #parentNode is parent node of N or None if node is root of the tree
    #localMaximaPoints is list of local maxima points sorted in descending order of reachability
    if len(localMaximaPoints) == 0:
        return #parentNode is a leaf
    
    #take largest local maximum as possible separation between clusters
    s = localMaximaPoints[0]
    node.assignSplitPoint(s)
    localMaximaPoints = localMaximaPoints[1:]

    #create two new nodes and add to list of nodes
    Node1 = TreeNode(RPoints[node.start:s],node.start,s, node)
    Node2 = TreeNode(RPoints[s+1:node.end],s+1, node.end, node)
    LocalMax1 = []
    LocalMax2 = []

    for i in localMaximaPoints:
        if i < s:
            LocalMax1.append(i)
        if i > s:
            LocalMax2.append(i)
    
    Nodelist = []
    Nodelist.append((Node1,LocalMax1))
    Nodelist.append((Node2,LocalMax2))
    
    #set a lower threshold on how small a significant maxima can be


    if RPlot[s] < minReach:
        node.assignSplitPoint(-1)
        #if splitpoint is not significant, ignore this split and continue
        clusterTree(node,parentNode, localMaximaPoints, RPlot, RPoints, min_cluster_size, minReach, maximaRatio)
        return
        
        
    #only check a certain ratio of points in the child nodes formed to the left and right of the maxima
    checkRatio = .8
    checkValue1 = int(np.round(checkRatio*len(Node1.points)))
    checkValue2 = int(np.round(checkRatio*len(Node2.points)))
    if checkValue2 == 0:
        checkValue2 = 1
    avgReachValue1 = float(np.average(RPlot[(Node1.end - checkValue1):Node1.end]))
    avgReachValue2 = float(np.average(RPlot[Node2.start:(Node2.start + checkValue2)]))


    '''
    To adjust the fineness of the clustering, adjust the following ratios.
    The higher the ratio, the more generous the algorithm is to preserving
    local minimums, and the more cuts the resulting tree will have.
    '''

    #the maximum ratio we allow of average height of clusters on the right and left to the local maxima in question
	
    if float(avgReachValue1 / float(RPlot[s])) >= maximaRatio and float(avgReachValue2 / float(RPlot[s])) >= maximaRatio:
        Nodelist.remove((Node1,LocalMax1))
        Nodelist.remove((Node2,LocalMax2))
        node.assignSplitPoint(-1)
        clusterTree(node,parentNode,localMaximaPoints,Rplot,Rpoints,min_cluster_size,minReach,maximaRatio)
        return
    elif float(avgReachValue1 / float(RPlot[s])) >= maximaRatio:
        Nodelist.remove((Node1,LocalMax1))
    elif float(avgReachValue2 / float(RPlot[s])) >= maximaRatio:
        Nodelist.remove((Node2,LocalMax2))

 
    #remove clusters that are too small
    if len(Node1.points) < min_cluster_size:
        #cluster 1 is too small"
        try:
            Nodelist.remove((Node1, LocalMax1))
        except Exception:
            sys.exc_clear()
    if len(Node2.points) < min_cluster_size:
        #cluster 2 is too small
        try:
            Nodelist.remove((Node2, LocalMax2))
        except Exception:
            sys.exc_clear()
    if len(Nodelist) == 0:
        #parentNode will be a leaf
        node.assignSplitPoint(-1)
        return
    
    '''
    Check if nodes can be moved up one level - the new cluster created
    is too "similar" to its parent, given the similarity threshold.
    Similarity is determined by the distribution of the reachibility of the cluster. If the cluster's average RD is less than one standard deviation from the parent node's average RD, then the new node becomes a child of the parent node.
    '''
    bypassNode = 0
    if parentNode != None:
        parentSD = np.std(RPlot[parentNode.start:parentNode.end])
        dist = np.average(RPlot[node.start:node.end]) - np.average(RPlot[parentNode.start:parentNode.end])
        if float(dist) < parentSD:
            parentNode.children.remove(node)
            bypassNode = 1
        
    for nl in Nodelist:
        if bypassNode == 1:
            parentNode.addChild(nl[0])
            clusterTree(nl[0], parentNode, nl[1], RPlot, RPoints, min_cluster_size,minReach,maximaRatio)
        else:
            node.addChild(nl[0])
            clusterTree(nl[0], node, nl[1], RPlot, RPoints, min_cluster_size,minReach,maximaRatio)
        

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

    plt.xlabel('Order of points')
    plt.ylabel('Reachability-distance')
    
    num = 2
    graphNode(root, num, ax)

    plt.savefig('Reachability_plot.png', dpi=None, facecolor='w', edgecolor='w',
      orientation='portrait', papertype=None, format=None,
     transparent=False, bbox_inches=None, pad_inches=0.1)
    plt.show()

            
def graphNode(node, num, ax):
    ax.hlines(num,node.start,node.end,color="red")
    for item in node.children:
        graphNode(item, num - .4, ax)

def automaticCluster(RPlot, RPoints, minpts, minReach=0.003, maximaRatio=0.75):

    min_cluster_size_ratio = .005
    min_neighborhood_size = 2
    min_maxima_ratio = 0.001
    
    if minpts == 0:
        min_cluster_size = int(min_cluster_size_ratio * len(RPoints))
        if min_cluster_size < 3:
            min_cluster_size = 3
    else:
        min_cluster_size == minpts
    
    
    nghsize = int(min_maxima_ratio*len(RPoints))

    if nghsize < min_neighborhood_size:
        nghsize = min_neighborhood_size
    
    localMaximaPoints = findLocalMaxima(RPlot, RPoints, nghsize)
    
    rootNode = TreeNode(RPoints, 0, len(RPoints), None)
    clusterTree(rootNode, None, localMaximaPoints, RPlot, RPoints, min_cluster_size,minReach,maximaRatio)
    return rootNode



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


if __name__ == "__main__":
    main(sys.argv[1:])
