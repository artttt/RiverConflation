
# coding: utf-8

# Paired river network conflation
# 
# Works on the theory that the defining property of any location in a river network is its catchment. Thus in another representation of the same river network the same 'river location' can be found by finding the most similar catchment.
# 
# Testing has shown that this method works well except where the river network line work has a level of detail that is not well supported well supported by the catchment data. Even in these cases the method may still have some value.
# 
# Please reference this publication whenever this code is used
# http://geomorphometry.org/system/files/Read2011geomorphometry.pdf
# Thanks
# 
# Testing has shown that it runs in a reasonable time/RAM usage with large datasets (approx 100000 subcatchments). However testing has been limited so YMMV.

# In[ ]:

# note for development
# I start ipython notebook with the  --script switch so that i get an importable .py version of this notebook on every save


# In[ ]:

from shapely.prepared import prep
from shapely.geometry import shape
from shapely.geometry import mapping,Polygon,Point,LineString
from shapely.ops import transform
import fiona
from rtree import index
import networkx as nx
import collections
import os


# first step is to produce a networkx multiDiGraph from your source data.
# 
# The code below does that specifically for the Australian geofabric dataset in ESRI geodatabase format
# 
# (note the geodatabase has to be upgraded in ESRI software first so that gdal can read it).
# 
# For other datasets it would need to be changed to produce the same output from the different inputs.
# 
# This needs to be done for both the conflation source and destination networks
# 
# the key thing is some sort of from and to node attribution to build the network structure between subcatchments. if this isnt available it can be generated from the coords of the stream links. Code for this not here but available on request.
# 
# subcatchments and streams need a one to one relationship. streams without a catchment are (i think) fully handled. catchments without streams are ignored.

# In[ ]:

def remove_geofabric_catch_duplicates(DG):
    '''deal with special case of duplicate catchment polygons in the geofabric
    
    Im not sure why the geofabric allows multiple stream features to 
    have a relationship to the same subcatchment.
    Maybe it is ok and i just need to change my expected data model..... will have to think about that.
    Its a bit inconvinient though. 
    Is there a good way to weed these out logically and consistently?
    Here is my attempt. Based on the fact that they seem to occur only at flow splits.
    '''
    ts = nx.topological_sort(DG)
    for n in ts:
        new_set = set()
        for f,t,k,data in DG.in_edges_iter(n,data=True,keys=True):
            new_set.update([(data['cid'])])
        for f,t,k,data in DG.out_edges_iter(n,data=True,keys=True):
            if data['cid'] in new_set:
                data['cid'] = None
                data['subCatch'] = Polygon()
    
    #unfortunatly many cases still remain
    #simple way of handling them is to make an arbitary choice
    #first in keeps the catchment
    #not ideal
    ss = set()
    for f,t,k,data in DG.edges_iter(data=True,keys=True):
        if data['cid'] is not None:
            if data['cid'] in ss:
                print 'WARNING: duplicate catchment removed catchment id = ' + str(data['cid'])
                data['cid'] = None
                data['subCatch'] = Polygon()
            #assert data['cid'] not in ss, 'duplicate catchment id ' + str(data['cid'])
            ss.add(data['cid'])

def read_geofabric_data(netGDB,bbox=None):
    catch = {}
    with fiona.open(netGDB, layer='AHGFCatchment') as c:
        for feat in c.items(bbox=bbox):
            geom = shape(feat[1]['geometry'])
            cid = feat[1]['properties']['HydroID']
            assert cid not in catch #shouldnt be duplicates 
            catch[cid] = geom
    
    DG=nx.MultiDiGraph()
    with fiona.open(netGDB, layer='AHGFNetworkStream') as c:
        for feat in c.items(bbox=bbox):
            streamLink = shape(feat[1]['geometry'])
             #for some reason these are coming in as multipart features with only one part - no need for this
            assert streamLink.type == 'MultiLineString'
            assert len(streamLink.geoms) == 1
            streamLink = streamLink.geoms[0]
            
            ##remove this - just here for testing
            #if streamLink.representative_point().y > -40.5: #tasmania for testing
            #    continue
            
            sid = feat[1]['properties']['HydroID']
            cid = feat[1]['properties']['DrainID']
            fid = feat[1]['properties']['From_Node']
            tid = feat[1]['properties']['To_Node']
            subCatch = catch.get(cid,Polygon())
            DG.add_edge(fid, tid, id=sid,cid=cid,subCatch=subCatch,stream=streamLink)
            
    return DG


# In[ ]:

def network_copy_offset(DG,distance=0.0):
    '''copies a network and shifts the subcatchment and stream coordinates by a distance
    
    used for testing conflation simpily using only one network
    Can be used to provide an understanding of the sensitivity of the network to positional error in features'''
    DG2 = DG.copy()
    if distance != 0:
        for f,t,k,data in DG2.edges_iter(data=True,keys=True):
            data['subCatch'] = transform(lambda x, y, z=None: (x+distance, y+distance), data['subCatch'])
            data['stream'] =  transform(lambda x, y, z=None: (x+distance, y+distance), data['stream'])
    return DG2


# In[ ]:

# remember rtrees are not picklable - doh
# TODO: work out the rtree bulk loading method - its quicker apparently
def build_index(DG):
    '''build spatial and other indexes'''
    
    index_generator = ((0, data['subCatch'].bounds,(f,t,k)) for f,t,k,data in DG.edges_iter(data=True,keys=True) if not data['subCatch'].is_empty)
    return index.Index(index_generator)

def build_index_slow(DG):
    DG_idx = index.Index()
    for f,t,k,data in DG.edges_iter(data=True,keys=True):
        if not data['subCatch'].is_empty:
            DG_idx.insert(0, data['subCatch'].bounds,obj=(f,t,k))
    return DG_idx


# In[ ]:

def upstream_edge_set_old_method(DG):
    '''build up a list of upstream edge ids as an attribute in the network
    
    WARNING: these can have a large memory footprint'''
    ts = nx.topological_sort(DG)
    for n in ts:
        new_set = set()
        for f,t,k,data in DG.in_edges_iter(n,data=True,keys=True):
            new_set.update((data['ids']))
        for f,t,k,data in DG.out_edges_iter(n,data=True,keys=True):
            data['ids'] = new_set.union([(f,t,k)])

def upstream_edge_set(DG):
    '''build up a list of upstream edge ids as an attribute in the network
    
    WARNING: these can have a large memory footprint'''
    for e,ids in upstream_edge_set_iter(DG):
        DG.get_edge_data(*e)['ids'] = ids
            
def upstream_edge_set_iter(DG):
    '''build up a list of upstream edge ids in the network
    
    returns ((f,t,k),ids)
    where (f,t,k) identifies the edge
    and ids is a set of edge identifiers in the the same style.
    
    avoids the large memory footprint of storing these
    and is more efficent then generating them at random for edges
    because it uses the network topology to avoid multiple network traversals.'''
    temp_dict = {}
    ts = nx.topological_sort(DG)
    for n in ts:
        new_set = set()
        for f,t,k,data in DG.in_edges_iter(n,data=True,keys=True):
            new_set.update(temp_dict[(f,t,k)])
            del temp_dict[(f,t,k)]
        for f,t,k,data in DG.out_edges_iter(n,data=True,keys=True):
            ids = new_set.union([(f,t,k)])
            temp_dict[(f,t,k)] = ids
            yield ((f,t,k),ids)
            
def node_upstream_edge_set(DG):
    '''build up a list of upstream edge ids in the destination network
    
    WARNING: these can have a large memory footprint'''
    for n in DG:
        new_set = set()
        for f,t,k,data in DG.in_edges_iter(n,data=True,keys=True):
            new_set.update((data['ids']))
        DG.node[n]['ids'] = new_set
        

#not used at present but might come in handy
def upstream_node_set(DG):
    '''build up a list of upstream nodes in the destination network
    
    WARNING: these can have a large memory footprint'''
    ts = nx.topological_sort(DG)
    for n in ts:
        new_set = set([n])
        for n2 in DG.predecessors_iter(n):
            new_set.update(DG.node[n2]['nids'])
        DG.node[n]['nids'] = new_set


# In[ ]:

def sub_catch_area(DG):
    '''calc subcatchment area'''
    for f,t,k,data in DG.edges_iter(data=True,keys=True):
        data['area'] = data['subCatch'].area
        
def catch_area_needs_ids(DG):
    '''calc catchment area
    
    works with anabranching networks without double counting'''
    sub_catch_area(DG)
    for f,t,k,data in DG.edges_iter(data=True,keys=True):
        data['catchArea'] = sum(DG.get_edge_data(*e)['area'] for e in data['ids'])
        
def catch_area_no_ids(DG):
    '''calc catchment area
    
    a simpiler but slower variation on catch_area
    doesnt need ids precomputed
    works with anabranching networks without double counting'''
    sub_catch_area(DG)
    for e,ids in upstream_edge_set_iter(DG):
        DG.get_edge_data(*e)['catchArea'] = sum(DG.get_edge_data(*e2)['area'] for e2 in ids)
        
def node_catch_area(DG):
    '''calc catchment area to nodes
    
    works with anabranching networks without double counting'''
    for n in DG:
        DG.node[n]['catchArea'] = sum(DG.get_edge_data(*e)['area'] for e in DG.node[n]['ids'])



def catch_area(DG):
    '''calc catchment area
    
    works with anabranching networks without double counting'''
    # does catch area for both nodes and edges
    # about 3 times faster than other area calculations
    #although a bit more complex to understand the code
    #doesnt need 'ids' from upstream_edge_set
    temp_dict = {}
    ts = nx.topological_sort(DG)
    for n in ts:
        new_set = set()
        for f,t,k,data in DG.in_edges_iter(n,data=True,keys=True):
            new_set.update(temp_dict[(f,t,k)])
            del temp_dict[(f,t,k)]
            DG.node[n]['catchArea'] = sum(a for a,f,t,k in new_set)
        for f,t,k,data in DG.out_edges_iter(n,data=True,keys=True):
            e_set = new_set.union([(data['subCatch'].area,f,t,k)])
            temp_dict[(f,t,k)] = e_set
            data['catchArea'] = sum(a for a,f,t,k in e_set)
            


# In[ ]:

def build_overlaps(DG1,DG2,DG2_idx):
    '''build up dictionary of overlapping areas between 2 graphs'''

    for f,t,data in DG1.edges_iter(data=True):
        geom = data['subCatch']
        data['overlaps'] = {}
        if not geom.is_empty:
            prepGeom = prep(geom)
            for e in DG2_idx.intersection(geom.bounds,objects='raw'):
                nGeom = DG2.get_edge_data(*e)['subCatch']
                #prepGeom = prep(geom) #im a little wary of prep i have had it grind things to a halt when reused a lot
                if prepGeom.intersects(nGeom):
                    area = geom.intersection(nGeom).area
                    if area > 0.0:
                        data['overlaps'][e] = area


# In[ ]:

def build_full_overlaps(DG):
    '''stores full overlaps for all edges
    
    use is not recomended due to excessive memory foot print
    use full_overlaps to produce what you need when you need it
    '''
    for f,t,k,data in DG.edges_iter(data=True,keys=True):
        data['fullOverlaps'] = full_overlaps(DG,data['ids'])
        
def full_overlaps(DG,edges):
    '''aggregate a table of overlapping area for set of edges.
    
    Usually done for the entire catchment above and including an edge.
    assumes build_overlaps has been run with DG as the first network
    works with anabranching river networks (MultiGraphs) without double counting
    '''
    fullOverlaps = collections.defaultdict(float)
    for e in edges:
        for e2, overlapArea in DG.get_edge_data(*e)['overlaps'].iteritems():
            fullOverlaps[e2] += overlapArea
    return fullOverlaps


# If speed or memory usage is a problem then there are a few changes that could be implemented
# 
# The first speed up might be to split DG1 up into seperate basins if possible. These could be run in parallel without any code changes.
# 
# While find all matches trys to be very efficent in a large network if you are only interested in a few catchment then you could write a similar routine that only ran for those catchments to save some time.
# 
# The 2 slow points are build_overlaps and find_all_matches (including full_overlaps).
# - could look at using cython for these. Although in build_overlaps the 85% of the time is spent in shapely/geos doing intersections
# - The next speed step is that both build_overlaps and find_all_matches could be easily implemented in parallel. would be best to prune the memory foot print of the required objects before copying for other processes (or use shared mem).
# - One memory saving is that there is no need to keep the spatial features in memory after build_overlaps is finished (need to extract node locations from streams in DG1 for find_matches).
# 

# In[ ]:

def find_all_matches(DG1,DG2,DG2_idx,searchRadius,sizeRatio=0,maxMatchKeep=1):
    '''for each catchment in DG1 find list the best matches

    see find_matches for more info
    '''
    matches = {}
    #build ids on the fly with the iterator to save a lot of memory
    for e,eids in upstream_edge_set_iter(DG1):
        m = find_matches(e,eids,DG1,DG2,DG2_idx,searchRadius,sizeRatio,maxMatchKeep)
        if m:
            matches[e] = m
    return matches

def find_matches(e,eids,DG1,DG2,DG2_idx,searchRadius,sizeRatio=0,maxMatchKeep=1):
    '''for a catchment in DG1 find list the best matches

    searchRadius
    Limits the search to subcatchments within this radius
    
    sizeRatio
    limits search and results to pairs that are similar in size
    ie a sizeRatio of 0.5 wont test a pair where one catchment is double the area of the other
    the closer to 1 the more the closer in area the 2 catchments need to be.
    Mostly a time saving tweak Used in conjunction with searchRadius allows a bigger search area without a big time penalty
        
    maxMatchKeep
    limits the number of potential matches that are kept to a short list of n items.
    The best items are chosen using a simple test of the quality of the match
    This is used to reduce the memory footprint of the returned dictionary
    defaults to 1, set to None to keep all
    '''
    data = DG1.get_edge_data(*e)
    # remember dont keep fullOverlaps - they use a lot of RAM
    fullOverlaps = full_overlaps(DG1,eids)
    if len(fullOverlaps) == 0:
        # not going to find anything so move on and save nothing
        return
    match = []
    searchBounds = Point(data['stream'].coords[-1]).buffer(searchRadius).bounds
    for e in DG2_idx.intersection(searchBounds,objects='raw'):
        data2 = DG2.get_edge_data(*e)
        if sizeRatio > min(data['catchArea'],data2['catchArea'])/max(data['catchArea'],data2['catchArea']):
            continue
        # for this pair of edges work out the area of the overlap (its a subset of full overlaps)
        overlap = sum( fullOverlaps.get(e2,0) for e2 in data2['ids'])
        if overlap > 0.0:
            qual = (overlap+overlap)/(data['catchArea']+data2['catchArea'])
            #quality score goes first so it is easy to sort by
            match.append((qual,overlap,e))
    if len(match) == 0:
        # nothing overlapping within the search bounds. This bit of the network is vastly different
        # dont save any result
        return
    #just keep the best few matches
    match.sort(reverse=True)
    return match[0:maxMatchKeep]
    


# In[ ]:

#untested - 

def find_all_node_matches(DG1,DG2,DG2_idx,searchRadius,maxMatchKeep=1):
    '''same as find_all_matches but working to a node not a subcatchment
    '''
    matches = {}
    for n in DG1:
        nIDs = DG1.node[n]['ids']
        fullOverlaps = full_overlaps(DG1,nIDs)
        if len(fullOverlaps) == 0:
            # not going to find anything so move on and save nothing
            continue
        match = []
        for f,t,k,data in DG1.in_edges_iter(n,data=True,keys=True):
            searchBounds = Point(data['stream'].coords[-1]).buffer(searchRadius).bounds
        for f2,t2,k2 in DG2_idx.intersection(searchBounds,objects='raw'):
            data2 = DG2.get_edge_data(f2,t2,k2)
            # for this pair of edges work out the area of the overlap (its a subset of full overlaps)
            overlap = sum( fullOverlaps.get(e2,0) for e2 in DG2.node[t2]['ids'])
            if overlap > 0.0:
                qual = (overlap+overlap)/(DG1.node[n]['catchArea']+DG2.node[t2]['catchArea'])
                #quality score goes first so it is easy to sort by
                match.append((qual,overlap,e))
        if len(match) == 0:
            # nothing overlapping within the search bounds. This bit of the network is vastly different
            # dont save any result
            continue
        #just keep the best few matches
        match.sort(reverse=True)
        matches[n] = match[0:maxMatchKeep]
    return matches


# In[ ]:

def best_match(matches):
    return [(k,l[0][2],l[0][0]) for k,l in matches.iteritems()]


# In[ ]:

def write_debug_lines(DG1,DG2,best,fileName):
    '''a simple output to show how each to node matches up
    
    handy for looking for errors in the conflation'''
    schema = {
        'geometry': 'LineString',
        'properties': {'qual': 'float'},
    }
    with fiona.open(fileName, 'w', 'ESRI Shapefile', schema) as c:
        for e1,e2,qual in best:
            p1 = DG1.get_edge_data(*e1)['stream'].coords[-1]
            p2 = DG2.get_edge_data(*e2)['stream'].coords[-1]
            geom = LineString(LineString([p1, p2]))
            c.write({
                'geometry': mapping(geom),
                'properties': {'qual': qual},
            })


# Note: not finished the reimplementation
# 
# TODO
# 
# unify the matches at a confluence to one match. difference between these and subcatchment match is also useful.
# 
# transfer/rewrite code that maps the conflation out along reaches and subcatchments
# 
# transfer/rewrite code that highlights topology differences
# 
# 
# 

# In[ ]:

#code to test the idea of bringing the conflation of catchments at confluence being conflated together.
#one way of doing this is done with find_all_node_matches.
#The way being tested here is to use the existing list of good conflations to look at all inflowing catchments
#to a confluence together and pick the best one.
#each catchment is given equal weighting (not area weighted!!)
#seems to work quite well.
def confluence_matches(DG1,matches):
    #TODO: tidy up this code
    node_matches = []
    for n in DG1.nodes_iter():
        in_matches = []
        for f,t,k,data in DG1.in_edges_iter(n,data=True,keys=True):
            if (f,t,k) in matches:
                in_matches.append (matches[(f,t,k)])
        if len(in_matches) ==0:
            continue

        #get a list of the potential toNodes in the second network
        uniqueNode = set()
        for m in in_matches:
            for qual,overlap,e in m:
                f2,t2,k2 = e
                uniqueNode.add(t2)
        n_dict = {}
        nScore = []
        for n2 in uniqueNode:
            n_dict[n2]= []
            for m in in_matches:
                #sometimes a toNode appears in 2 potential conflation matches - pick the best one
                ll = [(qual,overlap,e) for qual,overlap,e in m if e[1] == n2]
                if len(ll) > 0:
                    n_dict[n2].append (max(ll))

            s = sum(qual for qual,overlap,e in n_dict[n2])
            nScore.append((s/len(in_matches),n,n2))

        #print n_dict[max(nScore)[1]]
        node_matches.append(max(nScore))
        #print max(nScore),n
    return node_matches

def write_debug_lines_confluence_matches(DG1,DG2,node_matches,fileName):
    '''a simple output to show how each to node matches up
    
    handy for looking for errors in the conflation'''
    schema = {
        'geometry': 'LineString',
        'properties': {'qual': 'float'},
    }
    with fiona.open(fileName, 'w', 'ESRI Shapefile', schema) as c:
        for qual,n1,n2 in node_matches:
            for f,t,k,data in DG1.in_edges(n1,data=True,keys=True)[0:1]:
                p1 = data['stream'].coords[-1]
            for f,t,k,data in DG2.in_edges(n2,data=True,keys=True)[0:1]:
                p2 = data['stream'].coords[-1]
            geom = LineString(LineString([p1, p2]))
            c.write({
                'geometry': mapping(geom),
                'properties': {'qual': qual},
            })


# In[ ]:

def write_catch(DG,fileName):
    schema = {
        'geometry': 'Polygon',
        'properties': {'id': 'int',
                      'area': 'float'},
    }
    with fiona.open(fileName, 'w', 'ESRI Shapefile', schema) as c:
        for f,t,data in DG.edges_iter(data=True):
            if not data['subCatch'].is_empty:
                c.write({
                    'geometry': mapping(data['subCatch']),
                    'properties': {'id': data['id'],
                                    'area':data['catchArea']}
                    })

def write_stream(DG,fileName):
    schema = {
        'geometry': 'LineString',
        'properties': {'id': 'int',
                      'area': 'float'},
    }
    with fiona.open(fileName, 'w', 'ESRI Shapefile', schema) as c:
        for f,t,data in DG.edges_iter(data=True):
            if not data['stream'].is_empty:
                c.write({
                    'geometry': mapping(data['stream']),
                    'properties': {'id': data['id'],
                                   'area':data['catchArea']}
                })


# In[ ]:

#some code used for error testing in development
#assumes full overlap
#will get some show up due to floating point issues
def test_stuff(DG):
    for f,t,data in DG.edges_iter(data=True):
        area = data['subCatch'].area
        area2 = 0.0
        for k, v in data['overlaps'].iteritems():
                area2 += v
        if (abs(area - area2) > 0.000000000000000001):
            print area,area2

    #a test for errors
    for f,t,data in DG.edges_iter(data=True):
        area = data['catchArea']
        area2 = 0.0
        for k, v in data['fullOverlaps'].iteritems():
                area2 += v
        if (abs(area - area2) > 0.00000000000001):
            print area,area2

