
# coding: utf-8

# Please be aware that this example needs 5 - 6 GB of RAM and approx. 40 minutes to run

# In[ ]:

get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')


# In[ ]:

import riverConflation as rc
import os
import networkx as nx


# In[ ]:

#both datasets used for testing can be downloaded here
#http://www.bom.gov.au/water/geofabric/download.shtml
#note the geodatabases had to be upgraded in ArcCatalog to V10+ so that gdal/fiona can read them.
netGDB1 = r'C:\data\temp\SH_Network_GDB_V2_1_1\SH_Network_GDB\SH_Network.gdb'
netGDB2 =r'C:\data\temp\SH_Network_GDB_V3_0_PG\SH_Network_GDB\SH_Network.gdb'
pkl_path = r'C:\data\temp'
#checks all nodes within at least this distance for a conflation
#its a time saving optimisation and should be set on the large size to avoid missing the best result
#A guide to setting this is: set it a bit bigger then the largest expected separation of conflated nodes
#for small networks networks just set it very large
#in projection units (1degree approx 100km)
# a few kms is probably reasonable usually. Set low here for faster demo.
searchRadius = 0.01
#number of potential matches to keep for each feature
maxMatchKeep = 10


# First thing is to read in your data. This example uses the Australian Geofabric V2 and V3 data. Other datasets would need their own customised data prep code.
# 
# In the next 2 steps ignore the duplicate catchment warnings for the case of testing the code.
# I havent dealt with all the minor details of the geofabric quite right.

# In[ ]:

DG2 = rc.read_geofabric_data(netGDB2)
rc.remove_geofabric_catch_duplicates(DG2)
nx.write_gpickle(DG2, os.path.join(pkl_path, 'PG_conflation2.p'))
DG2_idx = rc.build_index(DG2)


# In[ ]:

DG1 = rc.read_geofabric_data(netGDB1,DG2_idx.bounds)
rc.remove_geofabric_catch_duplicates(DG1)
nx.write_gpickle(DG1, os.path.join(pkl_path, 'PG_conflation.p'))


# you can start from here by loading in the pickles that were created with the above code earlier.
# Run the imports and global variables code at the top first though

# In[ ]:

DG1 = nx.read_gpickle(os.path.join(pkl_path, 'PG_conflation.p'))
DG2 = nx.read_gpickle(os.path.join(pkl_path, 'PG_conflation2.p'))
DG2_idx = rc.build_index(DG2)
# starting from pickles = 1 minute 2GB
#starting from scratch = 12 minutes


# In[ ]:

#This is done seperate to finding matches because it takes a while so its nice to split it out for debugging
#%%timeit -r1 -n1
rc.build_overlaps(DG1,DG2,DG2_idx)
# 15 minutes 1GB


# In[ ]:

rc.catch_area(DG1)
rc.catch_area(DG2)
# 1 min


# The next step is to sum up areas to find the catchment overlap for every combination. We are only interested in the best or a short list of the overlaps that match well
# The simple approach is a brute force exhustive test of all combinations. This works well for a few thousand (75 minutes for 17k x 17k) sub catchments in each graph however it would not scale well as network sizes increase.
# 
# There are a few ways to reduce the set of catchments to test for a match. One issue to keep in mind is to not make assumptions about how similar the two networks topology might be. The approach taken is to use a tunable spatial proximity limit that should be set to a size that is expected to ensure finding the best matches within that radius, setting it too small would cause missed matches, too large would just take longer.
# There is also a limit on the pairs of catchments searched based on area similarity.this works well as good matches have to, by defininition, be of a similar size.
# 
# Generally the results from this method are not very sensitive to these parameters, if set conservatively, other then faster processing time.
# 

# In[ ]:




# In[ ]:

#%%timeit -n1 -r1
rc.upstream_edge_set(DG2)
#<10sec approx 2GB


# In[ ]:

sizeRatio=0.5
matches = rc.find_all_matches(DG1,DG2,DG2_idx,searchRadius,sizeRatio,maxMatchKeep)
# depends on the search radius used.
# 8 minutes with searchRadius=0.005 (500m) and a sizeRatio=0
# 7.5 minutes with searchRadius=0.01 (1km) and a sizeRatio=0.5


# In[ ]:

best = rc.best_match(matches)


# In[ ]:

# simple outputs
# more complete outputs still to be re implemented.... stay tuned.


# In[ ]:

rc.write_debug_lines(DG1,DG2,best,os.path.join(pkl_path, 'debug_lines.shp'))


# In[ ]:

# a more refined match of nodes considering each inflow to a confluence.
#find_all_matches needs to have saved a shortlist of matches by setting maxMatchKeep to something like 10
node_matches = rc.confluence_matches(DG1,matches)
rc.write_debug_lines_confluence_matches(DG1,DG2,node_matches,os.path.join(pkl_path, 'debug_lines_nodes.shp'))


# In[ ]:

rc.write_catch(DG1,os.path.join(pkl_path, 'catch1.shp'))
rc.write_catch(DG2,os.path.join(pkl_path, 'catch2.shp'))
rc.write_stream(DG1,os.path.join(pkl_path, 'stream1.shp'))
rc.write_stream(DG2,os.path.join(pkl_path, 'stream2.shp'))

