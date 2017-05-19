__author__ = 'Erwin Chen'
import math
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

def run():
    EPS = 0.35
    MIN_SAMPLES = 8
    N_ELEM = 9
    N_CHEM = 500
    N_RV = 300
    N_CUT = 1
    DATAFILE_NAME = 'results-unregularized-matched.fits'
    FEATURE_NAMES = ['APOGEE_ID', 'VHELIO_AVG', 'V_H', 'TI_H', 'CA_H', 'FE_H', 'K_H', 'MN_H', 'NI_H', 'SI_H', 'S_H']
    ELEMENT_NAMES = ['V_H', 'TI_H', 'CA_H','FE_H', 'K_H', 'MN_H', 'NI_H', 'SI_H', 'S_H']
    MEMBERFILE_NAME = 'table4.dat'

    ## load data from APOGEE
    ap_file = fits.open(DATAFILE_NAME)
    ap_data = ap_file[1].data
    feature_names = np.array(FEATURE_NAMES)
    element_names = np.array(ELEMENT_NAMES)
    elements = np.array([name.replace('_H', '').title() for name in element_names])
    print "The following elements are used for clustering: "
    print elements

    ## append data into columns
    ap_cols = []
    for name in feature_names:
        ap_cols.append(ap_data.field(name))
    ap_cols = np.array(ap_cols)
    ap_cols = ap_cols.T

    ## create a table with the columns
    dtype = ['float' for n in range(len(feature_names))]
    dtype[0] = 'string'
    ap_table = Table(data=ap_cols, names=feature_names, dtype=dtype)

    ## load membership file
    known_clusters = np.loadtxt(MEMBERFILE_NAME, usecols=(0, 1), dtype=('S', 'S'), unpack=True)
    member_IDs = known_clusters[0]
    member_names = known_clusters[1]
    labels = np.zeros(len(member_IDs))-1
    cluster_names = list(set(member_names))
    print "The following clusters are in the dataset: "
    print cluster_names

    ## add membership and numerical label to table
    k = 0
    for name in cluster_names:
        index = np.where(member_names == name)[0]
        labels[index] = k
        k += 1
    names = ['APOGEE_ID', 'cluster_name', 'label']
    dtype=['string', 'string', 'int']
    member_table = Table(data=[member_IDs, member_names, labels], names=names, dtype=dtype)
    ap_table = join(ap_table, member_table, keys='APOGEE_ID', join_type='left')

    ## fill missing values
    ap_table['cluster_name'].fill_value = 'background'
    ap_table['label'].fill_value = -1
    for element in element_names:
        ap_table[element].mask = np.isnan(ap_table[element])
        ap_table[element].fill_value = -9999.
    ap_table = ap_table.filled()

    ## get stars with valid values for all elements
    ap_stars = np.arange(len(ap_table))
    for element in element_names:
        ap_stars = np.intersect1d(ap_stars, np.where(ap_table[element] > -9999.)[0])
    ap_table = ap_table[ap_stars]
    print "There are %i stars with valid values for all elements."%len(ap_table)
    halo_index = np.where((ap_table['GLAT'] < -10.) | (ap_table['GLAT'] > 10.))[0]
    ap_table_halo = ap_table[halo_index]
    print "In the halo, there are %i stars with 15 elements."%len(ap_table_halo)

    ## get globular cluster members with valid values for all elements
    globular_names = np.array(['M107', 'M53', 'M92', 'M67', 'M5', 'M13', 'M3', 'M2', 'M15', 'N5466'])
    globular_members = np.array([], dtype='int')
    globular_labels = np.array([], dtype='int')
    k = 0
    save_list = []
    for name in globular_names:
        cluster_members = np.where(ap_table_halo['cluster_name'] == name)[0]
        if len(cluster_members) <= 0:
            print "%s has %i members"%(name, 0)
        else:
            cluster_labels = ap_table_halo['label'][cluster_members][0]
            globular_members = np.append(globular_members, cluster_members)
            globular_labels = np.append(globular_labels, cluster_labels)
            save_list.append(k)
            print "%s has %i members"%(name, len(cluster_members))
        k += 1
    globular_names = np.array([globular_names[i] for i in save_list])
    print "The following globular clusters are in the dataset: "
    print globular_names
    print "The numerical labels for globular clusters are as follows: "
    print globular_labels

    ## compose a matrix that contains chemical abundances and radial velocity
    Fe_index = np.where(element_names == 'FE_H')[0][0]
    chem = [ap_table[element]-ap_table['FE_H'] for element in element_names]
    chem[Fe_index] = ap_table['FE_H']
    chem.append(ap_table['VHELIO_AVG'])
    chem_RV = np.array(chem).T
    chem = np.delete(chem_RV,-1,1)
    print "The shape of the matrix is ",
    print chem_RV.shape
    
    ## get the nearest neighbors in chemical-velocity space
    indices = get_friends(chem_RV, len(elements), N_CHEM, N_RV)

    ## histogram over numbers of neighbors to find min_samples
    lengths = np.array([len(indices[n]) for n in range(len(indices))])
    H, edges = np.histogram(lengths)
    print "Number of Stars/ Threshold"
    for n in range(len(H)):
        print H[n], edges[n+1]

    ## select stars with more than N_cut neighbors
    non_noise = np.where(lengths > N_CUT)[0]
    print "%i stars will be used for clustering."%len(non_noise)

    ## show remaining globular clusters
    for name in globular_names:
        members_gc = np.where(ap_table_halo['cluster_name'] == name)[0]
        remain = np.intersect1d(members_gc, non_noise)
        if len(remain) <= 0:
            print "%s has %i members left"%(name, 0)
        else:
            print "%s has %i members, %.2f percent remaining"%(name, len(remain), len(remain)*100.0/len(members_gc))

    ## compose distance matrix
    S = lil_matrix((len(non_noise), len(non_noise)))
    for (n, m, dist) in iterator_dist(indices[non_noise]):    
        S[n,m] = dist   
        S[m,n] = dist

    ## DBSCAN clustering
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES, metric='precomputed', n_jobs=-1).fit(S, lengths[non_noise])
    labels = db.labels_
    n_clumps = np.amax(labels) + 1
    print "%i clusters found"%n_clumps
    print "#Categorized as Member/ Ratio of Member"
    print len(np.where(labels != -1)[0]), len(np.where(labels != -1)[0])*1.0/len(labels)

    pickle.dump(ap_table_halo, open("ap_table_halo.p", "wb"))
    pickle.dump(labels, open("SNN_DBSCAN_labels.p", "wb"))
    pickle.dump(non_noise, open("SNN_non_noise_stars.p", "wb"))
    pickle.dump(S, open("SNN_distance_matrix.p", "wb"))

if __name__ == '__main__':
    run()
