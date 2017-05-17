__author__ = 'Erwin Chen'
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.table import Table, join, Column
import numpy as np
import scipy.stats as stats
from scipy.sparse import lil_matrix
import scipy.spatial.distance as distance
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.neighbors import NearestNeighbors, DistanceMetric
from sklearn import mixture, preprocessing
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.cluster import v_measure_score, homogeneity_completeness_v_measure
import cPickle as pickle


def jaccard(a,b):
    """
    Calculate Jaccard distance between two arrays
    Parameters:
    -----------
    a: an array
    b: an array
    """
    A = np.array(a, dtype='int')
    B = np.array(b, dtype='int')
    A = A[np.where(A > -1)[0]]
    B = B[np.where(B > -1)[0]]
    union = np.union1d(A,B)
    intersection = np.intersect1d(A,B)
    return 1.0 - len(intersection)*1.0 / len(union)


def get_friends(data, num_element, n_chem, n_rv):
    """
    Get neartest neighbors in both chemical and radial velocity spaces for every star.
    Parameters:
    -----------
    data: a matrix that contains chemical abundances and radial velocities
    num_element: number of elements in the matrix
    n_chem: number of nearest neighbors in chemical space 
    n_rv: number of nearest neighbors in radial velocity space
    """
    data = np.array(data)
    index_chem = np.arange(0, num_element, 1)
    nbrs_chem = NearestNeighbors(n_neighbors=n_chem, algorithm='ball_tree', metric='manhattan').fit(data[:,index_chem])
    distances_chem, indices_chem = nbrs_chem.kneighbors(data[:,index_chem])
    index_rv = np.arange(num_element, len(data[0]), 1)
    rv_data = np.copy(data[:,index_rv])
    if len(rv_data[0]) < 2:
        rv_data = rv_data.reshape(-1, 1)
    nbrs_rv = NearestNeighbors(n_neighbors=n_rv, algorithm='ball_tree').fit(rv_data)
    distances_rv, indices_rv = nbrs_rv.kneighbors(rv_data)
    indices = []
    for n in range(len(indices_chem)):
        indices.append(np.intersect1d(indices_chem[n], indices_rv[n]))
    indices = np.array(indices)
    return indices


def iterator_dist(indices):
    """
    An iterator that calculates ans stores the Jaccard distance between every two stars
    Parameters:
    -----------
    indices: a list of indices of neighbors for every star
    """
    for n in range(len(indices)):
        for m in range(n+1, len(indices)):
            dist = jaccard(indices[n], indices[m])
            if dist < 1:
                yield (n, m, dist)


EPS = 0.35
MIN_SAMPLES = 8
N_ELEM = 9
N_CHEM = 500
N_RV = 300


# load data from APOGEE
ap_file = fits.open('results-unregularized-matched.fits')
ap_data = ap_file[1].data
feature_names = ['APOGEE_ID', 'GLON', 'GLAT', 'RA', 'DEC', 'VHELIO_AVG', 'LOGG', 'TEFF', 'PMRA', 'PMDEC', 
                 'AL_H', 'NA_H', 'O_H', 'MG_H','C_H', 'N_H', 'V_H', 'TI_H', 'CA_H','FE_H', 'K_H', 'MN_H', 'NI_H', 'SI_H', 'S_H', 
                 'SNR']
feature_names = np.array(feature_names)
element_names = ['AL_H', 'NA_H', 'O_H', 'MG_H','C_H', 'N_H', 'V_H', 'TI_H', 'CA_H','FE_H', 'K_H', 'MN_H', 'NI_H', 'SI_H', 'S_H']
element_names = np.array(element_names)
elements = np.array([name.replace('_H', '').title() for name in element_names])
print elements

# append data into columns
ap_cols = []
for name in feature_names:
    ap_cols.append(ap_data.field(name))
ap_cols = np.array(ap_cols)
ap_cols = ap_cols.T
# create a table with the data columns
dtype = ['float' for n in range(len(feature_names))]
dtype[0] = 'string'
ap_table = Table(data=ap_cols, names=feature_names, dtype=dtype)

VSCATTER = ap_data.field('VSCATTER')
VERR = ap_data.field('VERR')
VERR_MED = ap_data.field('VERR_MED')

print stats.describe(VSCATTER)
print stats.describe(VERR[np.where(VERR < 999999.0)[0]])
print stats.describe(VERR_MED[np.where(VERR_MED < 999999.0)[0]])


# In[4]:

# add membership and number labels for clusters
known_clusters = np.loadtxt('table4.dat', usecols=(0, 1), dtype=('S', 'S'), unpack=True)
member_IDs = known_clusters[0]
member_names = known_clusters[1]
labels = np.zeros(len(member_IDs))-1
cluster_names = list(set(member_names))
print cluster_names
k = 0
for name in cluster_names:
    index = np.where(member_names == name)[0]
    labels[index] = k
    k += 1
names = ['APOGEE_ID', 'cluster_name', 'label']
dtype=['string', 'string', 'int']
member_table = Table(data=[member_IDs, member_names, labels], names=names, dtype=dtype)
ap_table = join(ap_table, member_table, keys='APOGEE_ID', join_type='left')

# fill missing values
ap_table['cluster_name'].fill_value = 'background'
ap_table['label'].fill_value = -1
for element in element_names:
    ap_table[element].mask = np.isnan(ap_table[element])
    ap_table[element].fill_value = -9999.
ap_table = ap_table.filled()

# get stars with 15 elements
ap_stars_15 = np.arange(len(ap_table))
for element in element_names:
    ap_stars_15 = np.intersect1d(ap_stars_15, np.where(ap_table[element] > -999.)[0])
ap_table_15 = ap_table[ap_stars_15]
print "There are %i stars with 15 elements."%len(ap_table_15)
halo_index = np.where((ap_table_15['GLAT'] < -10.) | (ap_table_15['GLAT'] > 10.))[0]
ap_table_15 = ap_table_15[halo_index]
print "In the halo, there are %i stars with 15 elements."%len(ap_table_15)

# get globular cluster members with 15 elements
globular_names = np.array(['M107', 'M53', 'M92', 'M67', 'M5', 'M13', 'M3', 'M2', 'M15', 'N5466'])
globular_members = np.array([], dtype='int')
globular_labels = np.array([], dtype='int')
k = 0
save_list = []
for name in globular_names:
    cluster_members = np.where(ap_table_15['cluster_name'] == name)[0]
    if len(cluster_members) <= 0:
        print "%s has %i members"%(name, 0)
    else:
        cluster_labels = ap_table_15['label'][cluster_members][0]
        globular_members = np.append(globular_members, cluster_members)
        globular_labels = np.append(globular_labels, cluster_labels)
        save_list.append(k)
        print "%s has %i members"%(name, len(cluster_members))
    k += 1
globular_names = np.array([globular_names[i] for i in save_list])
print globular_names
print globular_labels

Fe_index = np.where(element_names == 'FE_H')[0][0]
chem = [ap_table_15[element]-ap_table_15['FE_H'] for element in element_names]
chem[Fe_index] = ap_table_15['FE_H']
chem.append(ap_table_15['VHELIO_AVG'])
use_chem_RV = np.array(chem[6:]).T
print use_chem_RV.shape


indices = get_friends(use_chem_RV, 9, N_CHEM, N_RV)


lengths = np.array([len(indices[n]) for n in range(len(indices))])
H, edges = np.histogram(lengths)
for n in range(len(H)):
    print H[n], edges[n+1]

for name in globular_names:
    c_members = np.where(ap_table_15['cluster_name'] == name)[0]
    print name
    print lengths[c_members]

non_noise = np.where(lengths > 1)[0]
print len(non_noise)
# show remaining globular clusters
for name in globular_names:
    members_gc = np.where(ap_table_15['cluster_name'] == name)[0]
    remain = np.intersect1d(members_gc, non_noise)
    if len(remain) <= 0:
        print "%s has %i members"%(name, 0)
    else:
        print "%s has %i members"%(name, len(remain)), len(remain)*1.0/len(members_gc)



S = lil_matrix((len(non_noise), len(non_noise)))
for (n, m, dist) in iterator_dist(indices[non_noise]):    
    S[n,m] = dist   
    S[m,n] = dist

for name in globular_names:
    members = np.where(ap_table_15["cluster_name"][non_noise] == name)[0]
    print name
    d = S[members][:,members].toarray()
    for row in d:
        d_nz = np.sort(row[row != 0])
        if len(d_nz) > 4:
            print d_nz[:3]

db = DBSCAN(eps=0.35, min_samples=8, metric='precomputed', n_jobs=-1).fit(S, lengths[non_noise])
labels = db.labels_
n_clumps = np.amax(labels)
true_labels = ap_table_15["label"][non_noise]
print n_clumps
print len(np.where(labels != -1)[0]), len(labels), len(np.where(labels != -1)[0])*1.0/len(labels)

# get recovered members
recovered = np.where([])
for n in range(n_clumps):
    group = np.where(labels == n)[0]
    group_labels = true_labels[group]
    members_in_group = np.where(group_labels != -1)[0]
    numbers = list(set(group_labels[members_in_group]))
    for num in numbers:
        len_cluster = len(np.where(group_labels[members_in_group] == num)[0])
        len_clusters = len(members_in_group)
        if (len_cluster == len_clusters):
            recovered = np.append(recovered, non_noise[group[np.where(group_labels == num)[0]]])
recovered = np.intersect1d(recovered, globular_members)
print ap_table_15['cluster_name'][recovered]

for name in globular_names:
    members = np.where(ap_table_15["cluster_name"] == name)[0]
    rec_member = np.intersect1d(members, recovered)
    print name, len(rec_member)*1.0/len(members), len(members)

for name in globular_names:
    members = np.where(ap_table_15["cluster_name"] == name)[0]
    print name, np.std(ap_table_15['VHELIO_AVG'][members])
#     plt.hist(ap_table_15['VHELIO_AVG'][members])
#     plt.show()


pickle.dump(labels, open("ISNN_500_300_5_8_labels.p", "wb"))
pickle.dump(non_noise, open("ISNN_500_300_5_8_stars.p", "wb"))
pickle.dump(S, open("ISNN_500_300_5_8_S.p", "wb"))




