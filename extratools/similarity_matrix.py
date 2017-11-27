#~language: Python 2.7.6
#~libraries:
#numpy 1.8.2
#scikit-learn (sklearn) 0.19.0
#~~~~> author: Alexander Ramos Diaz
#~description: this small code is an example of how to compute a similarity
#matrix, it requires numpy and scikit-learn, which is a python library for
#machine learning, it uses the following equation to compute the euclidean
#distances: dist(x, y) = sqrt(dot(x, x) - 2 * dot(x, y) + dot(y, y)) between.
#for more information go to:
#http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.euclidean_distances.html
#For specific and more precise computations, scipy.spatial.distance is recommended.
# The objective of the code is to obtain the similarity matrix from samples:
# Sm = I - (1/d)E
# Where Sm is the similarity matrix, E is the euclidian distances matrix obtained
# from a sample matrix D.

import numpy as np
import time

print "This is an example of how to compute a similarity matrix:"
print

#~step0~ this is only an example matrix (you also can load your matrix from a file):
samples = np.array([[0.2, 3.4, 5.6, 7.8,3.1], [9.1, 11.2, 13.4,1.6,2.72],[1.5,2.3,3.3,7.6,8.9]], np.float64)
print "Example input matrix:"
print samples
print

#~step1~ computing similarity matrix using scikit-learn:
M1 = np.asmatrix(samples)
from sklearn.metrics.pairwise import euclidean_distances
EucMatrix=euclidean_distances(M1, M1)
OneMatrix=np.ones_like(EucMatrix)
distmax = np.amax(EucMatrix)
SiMatrix = OneMatrix - EucMatrix/(distmax)
SiMatrix = np.asarray(SiMatrix)

#~step2~ shows the new matrix under similarity transformation:
print "Output similarity matrix using scikit-learn:"
print SiMatrix
print
print "run time =" , time.clock() #shows runtime


#~step3~ computing similarity matrix using scipy
M2 = np.asmatrix(samples)
from scipy.spatial.distance import pdist
M2distances = pdist(M2,metric='euclidean')
print type(M2distances), M2distances.shape
onesmatrix = np.ones_like(M2distances)
print type(onesmatrix), onesmatrix.shape
maxdistance = np.amax(M2distances)
M2similarity = onesmatrix - M2distances/(maxdistance)
from scipy.spatial.distance import squareform
M2similarity = squareform(M2similarity)
M2similarity = np.asarray(M2similarity)

#~step4~ shows the new matrix under similarity transformation:
print "Output similarity matrix using scipy:"
print M2similarity
print
print "run time =" , time.clock() #shows runtime
