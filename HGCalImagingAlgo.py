##############################################################################
# Implementation of HGCalImagingAlgo functionality (stand-alone)
# based on the CMSSW implementation mainly in RecoLocalCalo/HGCalRecAlgos
##############################################################################
# needed for ROOT funcs/types
import ROOT
import math
# needed for KDTree indexing & searches
import numpy as np
from scipy import spatial
# needed to extend the maximum recursion limit, for large data sets
import sys
sys.setrecursionlimit(100000)

## basic setup for testing
# 2D clustering
delta_c = 2.0
kappa = 10.
ecut = 0.060
# multi-clustering
multiclusterRadius = 0.015
realSpaceCone = False
minClusters = 3
# det. layers to consider
det_layers = 40
# others
verbosityLevel = 0 # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed

# definition of Hexel element
class Hexel:
    def __init__(self, rHit = None):
        self.eta = 0
        self.phi = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.isHalfCell = False
        self.weight = 0
        self.fraction = 1
        self.detid = None
        self.rho = 0
        self.delta = 0
        self.nearestHigher = -1
        self.isBorder = False
        self.isHalo = False
        self.clusterIndex = -1
        if rHit is not None:
            self.eta = rHit.eta
            self.phi = rHit.phi
            self.x = rHit.x
            self.y = rHit.y
            self.z = rHit.z
            self.weight = rHit.energy
            self.detid = rHit.detid
    def __gt__(self, other_rho):
        return self.rho > other_rho

# definition of basic cluster (based on a set of sub-clusters or set of hexels)
class BasicCluster:
    def __init__(self, energy = None, position = None, thisCluster = None, algoId = None, caloId = None):
        self.eta = 0
        self.phi = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.energy = 0
        self.thisCluster = None
        self.algoId = None
        self.caloId = None
        if energy is not None:
            self.energy = energy
        if position is not None:
            self.eta = position.eta()
            self.phi = position.phi()
            self.x = position.x()
            self.y = position.y()
            self.z = position.z()
        if algoId is not None:
            self.algoId = algoId
        if caloId is not None:
            self.caloId = caloId
        if thisCluster is not None:
            self.thisCluster = thisCluster

# distance squared (in eta/phi) between the two objects (hexels, clusters)
def distanceDR2(Hex1, Hex2):
    return (pow(Hex2.eta - Hex1.eta,2) + pow(Hex2.phi - Hex1.phi,2))

# distance squared (in x/y) between the two objects (hexels, clusters)
def distanceReal2(clust1, clust2):
    return (pow(clust2.x - clust1.x,2) + pow(clust2.y - clust1.y,2))

# position of the cluster, based on hexels positions weighted by the energy
def calculatePosition(cluster):
    total_weight = 0.
    x = 0.
    y = 0.
    z = 0.
    for iNode in cluster:
        if(not iNode.isHalo):
            total_weight += iNode.weight
            x += iNode.x*iNode.weight
            y += iNode.y*iNode.weight
            z += iNode.z*iNode.weight
    return ROOT.Math.XYZPoint( x/total_weight, y/total_weight, z/total_weight ) # return as ROOT.Math.XYZPoint
#    return x/total_weight, y/total_weight, z/total_weight

# calculate max local density in a 2D plane of hexels
def calculateLocalDensity(nd, lp):
    maxdensity = 0
    for iNode in nd:
        # search in a circle of radius delta_c (not identical to search in the box delta_c)
        found = lp.query_ball_point([iNode.x,iNode.y],delta_c)
        for j in found:
            if(distanceReal2(iNode,nd[j]) < delta_c*delta_c):
                iNode.rho += nd[j].weight
                if(iNode.rho > maxdensity):
                    maxdensity = iNode.rho
    return maxdensity

# calculate distance to the nearest hit with higher density (still does not use KDTree)
def calculateDistanceToHigher(nd, lp):
    #sort vector of Hexels by decreasing local density
    rs = sorted(range(len(nd)), key=lambda k: nd[k].rho, reverse=True)

    # intial values, and check if there are any hits
    maxdensity = 0.0
    nearestHigher = -1
    if(len(nd)>0):
        maxdensity = nd[rs[0]].rho
    else:
        return maxdensity # there are no hits

    #   start by setting delta for the highest density hit to the most distant hit - this is a convention
    dist2 = 2500.0
    for jNode in nd:
        tmp = distanceReal2(nd[rs[0]], jNode)
        if(tmp > dist2):
            dist2 = tmp
    nd[rs[0]].delta = pow(dist2,0.5)
    nd[rs[0]].nearestHigher = nearestHigher
            
    # now we save the largest distance as a starting point
    max_dist2 = dist2
    # calculate all remaining distances to the nearest higher density
    for oi in range(1,len(nd)): # start from second-highest density
        dist2 = max_dist2
        # we only need to check up to oi since hits are ordered by decreasing density
        # and all points coming BEFORE oi are guaranteed to have higher rho and the ones AFTER to have lower rho
        for oj in range(0,oi):
            tmp = distanceReal2(nd[rs[oi]], nd[rs[oj]])
            if(tmp <= dist2): #this "<=" instead of "<" addresses the (rare) case when there are only two hits
                dist2 = tmp
                nearestHigher = rs[oj]
        nd[rs[oi]].delta = pow(dist2,0.5)
        nd[rs[oi]].nearestHigher = nearestHigher #this uses the original unsorted hitlist

    return maxdensity

# find cluster centers that satisfy delta & maxdensity/kappa criteria, and assign coresponding hexels
def findAndAssignClusters(nd, points_0, points_1, lp, maxdensity):
    clusterIndex = 0
    #sort Hexels by decreasing local density and by decreasing distance to higher
    rs = sorted(range(len(nd)), key=lambda k: nd[k].rho, reverse=True) # indices sorted by decreasing rho
    ds = sorted(range(len(nd)), key=lambda k: nd[k].delta, reverse=True) # sort in decreasing distance to higher
    
#    for i in range(0,len(nd)):
#        print "index rs[i]: ", rs[i], ", rho: ", nd[rs[i]].rho, " delta: ", nd[rs[i]].delta, ", nearestHigher: ", nd[rs[i]].nearestHigher, ", eta: ", nd[rs[i]].eta, ", phi: ", nd[rs[i]].phi

    for i in range(0,len(nd)):
        if(nd[ds[i]].delta < delta_c):
            break # no more cluster centers to be looked at
        if(nd[ds[i]].rho < maxdensity/kappa):
            continue # skip this as a potential cluster center because it fails the density cut
        nd[ds[i]].clusterIndex = clusterIndex
        if (verbosityLevel>=2):
            print "Adding new cluster with index ", clusterIndex
            print "Cluster center is hit ", ds[i], " with density rho: ", nd[ds[i]].rho, "and delta: ", nd[ds[i]].delta, "\n"
        clusterIndex += 1
        
    # at this point clusterIndex is equal to the number of cluster centers - if it is zero we are done
    if(clusterIndex==0):
        return clusterIndex
    current_clusters = [[] for i in range(0,clusterIndex)]

    # assign to clusters, using the nearestHigher set from previous step (always set except for top density hit that is skipped)...
    for oi in range(1,len(nd)):
        ci = nd[rs[oi]].clusterIndex
        if(ci == -1):
            nd[rs[oi]].clusterIndex =  nd[nd[rs[oi]].nearestHigher].clusterIndex

    # assign points closer than dc to other clusters to border region and find critical border density
    rho_b = [0. for i in range(0,clusterIndex)]
    lp = spatial.KDTree(zip(points_0, points_1), leafsize=1000) # new KDTree
    # now loop on all hits again :( and check: if there are hits from another cluster within d_c -> flag as border hit
    for iNode in nd:
        ci = iNode.clusterIndex
        flag_isolated = True
        if(ci != -1):
            found = lp.query_ball_point([iNode.x,iNode.y],delta_c)
            for j in range(1,len(found)):
                # check if the hit is not within d_c of another cluster
                if(nd[j].clusterIndex!=-1):
                    dist2 = distanceReal2(nd[j],iNode)
                    if(dist2 < delta_c*delta_c and nd[j].clusterIndex!=ci):
                        # in which case we assign it to the border
                        iNode.isBorder = True
                        break
                    # because we are using two different containers, we have to make sure that we don't unflag the
                    # hit when it finds *itself* closer than delta_c
                    if(dist2 < delta_c*delta_c and dist2 != 0. and nd[j].clusterIndex==ci):
                    # this is not an isolated hit
                        flag_isolated = False
            if(flag_isolated):
                iNode.isBorder = True # the hit is more than delta_c from any of its brethren
        # check if this border hit has density larger than the current rho_b and update
        if(iNode.isBorder and rho_b[ci] < iNode.rho):
            rho_b[ci] = iNode.rho

    # flag points in cluster with density < rho_b as halo points, then fill the cluster vector
    for iNode in nd:
        ci = iNode.clusterIndex
        if(ci!=-1 and iNode.rho < rho_b[ci]):
            #iNode.isHalo = True
            pass # temporarly disabled until debugged (it seems that it does not work for eta<0)
        if(ci!=-1):
            current_clusters[ci].append(iNode)
            if (verbosityLevel>=2):
                print "Pushing hit ", iNode, " into cluster with index ", ci
                print "   rho_b[ci]: ", rho_b[ci], ", iNode.rho: ", iNode.rho, " iNode.isHalo: ", iNode.isHalo

    return current_clusters

# make 2D clusters out of rechists (need to introduce class with input params: delta_c, kappa, ecut, ...)
def makeClusters(rHitsCollection, ecut = ecut):
    # init 2D hexels lists
    points = [[] for i in range(0,det_layers)] # initialise list of per-layer-lists of hexels
    clusters = [[] for i in range(0,det_layers)] # initialise list of per-layer-clusters

    # loop over all hits and create the Hexel structure, skip energies below ecut
    for rHit in rHitsCollection:
        if (rHit.layer >= det_layers): continue # current protection
        if(rHit.energy < ecut): continue
        points[rHit.layer].append(Hexel(rHit))

    # loop over all layers, and for each layer create a list of clusters
    for layer in range(0, det_layers):
        if (len(points[layer]) == 0): continue # protection
        points_0 = [hex.x for hex in points[layer]] # list of hexels'coordinate 0 for current layer
        points_1 = [hex.y for hex in points[layer]] # list of hexels'coordinate 1 for current layer
        hit_kdtree = spatial.KDTree(zip(points_0, points_1), leafsize=1000) # create KDTree
        maxdensity = calculateLocalDensity(points[layer], hit_kdtree) # get the max density
        #print "layer: ", layer, ", max density: ", maxdensity, ", total hits: ", len(points[layer])
        calculateDistanceToHigher(points[layer], hit_kdtree) # get distances to the nearest higher density
        clusters[layer] = findAndAssignClusters(points[layer], points_0, points_1, hit_kdtree, maxdensity) # get clusters per layer
        #print "found: ", len(clusters[layer]), " clusters."

    # return the clusters list
    return clusters

# get basic clusters from the list of 2D clusters
def getClusters(clusters):
    # init the lists
    thisCluster = []
    clusters_v = []
    # loop over all layers and all clusters in each layer
    layer = 0
    for clist_per_layer in clusters:
        index = 0
        for cluster in clist_per_layer:
            energy = 0
            position = calculatePosition(cluster)
            for iNode in cluster:
                if (not iNode.isHalo):
                    energy += iNode.weight
            if (verbosityLevel>=1):
                print "Layer: ", layer, "| 2D-cluster index: ", index, ", No. of cells = ", len(cluster), ", Energy  = ", energy, ", Phi = ", position.phi(), ", Eta = ", position.eta(), ", z = ", position.z()
                for iNode in cluster:
                    if (not iNode.isHalo):
                        print "Layer: ", layer, "|                    ",       "  detid = ", iNode.detid, ", weight  = ", iNode.weight, ", phi = ", iNode.phi, ", eta = ", iNode.eta

            clusters_v.append(BasicCluster(energy = energy, position = position, thisCluster = cluster))
            index += 1
        layer += 1
    return clusters_v

# get position of the multi-cluster, based on the positions of its 2D clusters weighted by the energy
def getMultiClusterPosition(multi_clu, vz):
    if(len(multi_clu) == 0): return ROOT.Math.XYZPoint()
    acc_rho = 0.0
    acc_eta = 0.0
    acc_phi = 0.0
    totweight = 0.
    for layer_clu in multi_clu:
        x = layer_clu.x
        y = layer_clu.y
        point_r2 = (x*x + y*y)
        point_z = layer_clu.z-vz
        point_h = pow(point_r2 + point_z*point_z,0.5)
        weight = layer_clu.energy * len(layer_clu.thisCluster) # need to check this (is it weight = energy * size ?)
        if not (y != 0. or x != 0.): print "Cluster position somehow in beampipe."
        if not (point_z != 0.): print "Layer-cluster position given as reference point."
        point_r = pow(point_r2,0.5)
        acc_rho += point_r * weight
        acc_phi += math.atan2(y,x) * weight
        acc_eta += -1. * math.log(point_r/(point_z + point_h)) * weight
        totweight += weight
    invweight = 1.0/totweight
    temp = ROOT.Math.RhoEtaPhiPoint(acc_rho*invweight,acc_eta*invweight,acc_phi*invweight)
    return ROOT.Math.XYZPoint(temp.x(),temp.y(),temp.z())

# get energy of the multi-cluster, based on its 2D clusters
def getMultiClusterEnergy(multi_clu):
    acc = 0.
    for layer_clu in multi_clu:
        acc += layer_clu.energy
    return acc

# make multi-clusters stasrting from the 2D clusters
def makePreClusters(clusters, multiclusterRadius = multiclusterRadius, minClusters = minClusters):
    # get clusters in one list (just following original approach)
    thecls = getClusters(clusters)

    # init lists and vars
    thePreClusters = []
    vused = [0.]*len(thecls)
    used = 0
    # indices sorted by decreasing energy
    es = sorted(range(len(thecls)), key=lambda k: thecls[k].energy, reverse=True)
    # loop over all clusters
    index = 0
    for i in range(0,len(thecls)):
        if(vused[i]==0):
            temp = [thecls[es[i]]]
            if (thecls[es[i]].z>0): vused[i] = 1
            else: vused[i] = -1
            used += 1
            for j in range(i+1,len(thecls)):
                if(vused[j]==0):
                    distanceCheck = 9999.
                    if(realSpaceCone):
                        distanceCheck = distanceReal2(thecls[es[i]],thecls[es[j]])
                    else:
                        distanceCheck = distanceDR2(thecls[es[i]],thecls[es[j]])
                    if( distanceCheck < multiclusterRadius*multiclusterRadius and int(thecls[es[i]].z*vused[i])>0 ):
                        temp.append(thecls[es[j]])
                        vused[j] = vused[i]
                        used += 1
            if(len(temp) > minClusters):
                position = getMultiClusterPosition(temp, 0)
                energy = getMultiClusterEnergy(temp)
                thePreClusters.append(BasicCluster(energy = energy, position = position, thisCluster = temp))
                print "Multi-cluster index: ", index, ", No. of 2D-clusters = ", len(temp), ", Energy  = ", energy, ", Phi = ", position.phi(), ", Eta = ", position.eta(), ", z = ", position.z()
                index += 1
    return thePreClusters

