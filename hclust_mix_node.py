"""
Hopfield clustering

Reads expression data, creates a Hopfield network and creates various plots
such as the weight matrix, the energy surface and the relaxation of the
network.

Version : 1.00
Author  : Stefan Maetschke

hclust_mix_node is only a version of hclust that produces:
- output text files for further analysis
- images as .png files

This version works in clusters, by using loop_hclustmix_node.sh (available in the
extratools folder) you can run this code in a cluster.

.st (state matrices) and .trk (plotted tracks in the PCA) files are outputs for
further analysis such as attractors state identification and state comparison.

Alexander Ramos Diaz
"""
#Python 2.7.6 and:
#numpy 1.8.2
#scipy 0.13.3
#matplotlib 0.19.0
#pandas 0.17.1
#scikit-learn 0.19.0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:
#The order of how libraries are listed from 1 to 4 allows image output files.
import matplotlib #1
matplotlib.use('Agg') #2
from pylab import * #3
import matplotlib.pyplot as plt #4
import numpy as np
import time
import pandas as pd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def load_data(filepath, do_norm, do_log='AUTO'):
    """Load tab separated expression data.
       First row with sample ids is ignored.
       Second row needs to contains subtype labels
       Columns will be sorted according to subtype label.
       filepath -- expression data in tab format
       do_norm -- normalizes expression data if true
       do_log  -- log2 transformation of expression data if true.
                  if do_log=='AUTO' log2 transformation is performed if
                  skew of data is > 1.
       returns a tuple (sample_labels, gene_names, expression_matrix)
       where expression matrix contains samples in rows and genes in columns
    """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #defining global variable genenames and specieslabel
    global genenames
    global specieslabel

   #this label will be used for output file (specific of the data set):
    specieslabel = filepath[9:11]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def trim(label): return label.replace('"','')
    #print "Loading data ..."
    with open(filepath) as f:
        f.next() #skip sample id
        mat = [line.rstrip().split('\t') for line in f]
    mat = zip(*sorted(zip(*mat), key = lambda c: c[0]))
    genes = [row[0] for row in mat[1:]]
    labels = map(trim,mat[0][1:])
    data = array([map(float,row[1:]) for row in mat[1:]], dtype=float).T
    if do_norm: data = normalize(data, do_log)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:

    #assigns the list genes to the global list genenames for further use
    genenames = genes

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return labels, genes, data


def log2_transform(samples):
    """gene-wise log2 transformation of data. Typically required for
       single-channel data.
       samples -- expression matrix (samples in rows)
       returns log2 transformed expression matrix
    """
    #print "Log2 transformation ..."
    def save_median(sample):
        m = median(g)
        return 1.0 if abs(m)<1e-8 else m
    save_log2 = vectorize(lambda x: x if x<=0 else log2(x))
    return array([save_log2(g/save_median(g)) for g in samples.T]).T


def z_score_transform(samples):
    """sample-wise z-score normalization of data.
       samples -- expression matrix (samples in rows)
       returns z-score normalized expression matrix
    """
    def save_std(sample):
        s = std(sample)
        return 1.0 if s<1e-5 else s
    return array([(s-mean(s))/save_std(s) for s in samples])


def normalize(samples, do_log=False):
    """Normalize expressiond ata
       samples -- expression matrix (samples in rows)
       do_log -- log2 transform data
                 True, False, 'AUTO'
                 if do_log=='AUTO' log2 transformation is performed if
                 skew of data is > 1.
       returns normalized expression matrix
    """
    from scipy.stats import skew
    #print "Normalizing data ..."
    if do_log or (do_log=='AUTO' and skew(samples.T, axis=None) > 1.0):
        samples = log2_transform(samples)
    samples = z_score_transform(samples)
    return samples


def class2idx(annotations):
    """Convert string labels to numerical indices (=class ids)
       annotations -- list of (string) annotations
       returns list of numerical annotations [0, len(set(a))[
    """
    idx = {l:i for i,l in enumerate(set(annotations))}
    return array([idx[l] for l in annotations])


def idx2color(idx):
    """Convert numerical value to color code.
       idx -- integer value
       returns color code
    """
    colors = {0:'b', 1:'r', 2:'g', 3:'c', 4:'m', 5:'k', 6:'y'}
    return colors[idx%7]


def states2clusterids(states):
    """Convert state vectors to cluster class labels.
       states -- list of network states
       returns map from state vector to cluster class label.
    """
    def statestr(state):
        return "".join(['0' if v<0 else '1' for v in state])
    statestrs = [statestr(s) for s in states]
    n = len(set(statestrs))
    idmap = dict(zip(set(statestrs),range(n)))
    return array([idmap[s] for s in statestrs])


def feature_selection(data, n=None):
    """Selects genes with largest variance.
       data -- tuple (labels, genes, samples) from load_data()
       n -- number of features to select
            if None the ellbow of the variance over #features plot is used.
       returns tuple (labels, genes, samples) with selected features
    """
    def find_best_n(scores):
        from numpy.linalg import norm
        pts = array(list(enumerate(sorted(scores,reverse=True))))
        s,e = pts[0,:],pts[-1,:]
        b = e-s
        bn = b/norm(b)
        ds = [(p,norm((p-s)-dot((p-s),bn)*bn)) for p in pts]
        p,d = max(ds, key=lambda (p,d): d)
        return int(p[0])
    #print "Selecting features ..."
    labels, genes, samples = data
    c = class2idx(labels)
    nc = len(set(labels))
    scores = samples.std(axis=0)
    tuples = zip(genes,samples.T, scores)
    tuples.sort(key=lambda(g,s,v): -v)
    if not n: n = find_best_n(scores)
    genes, samples, _ = zip(*(tuples[:n]))
    return labels, genes, vstack(samples).T


"""vectorized signum function"""
signum = vectorize(lambda x: -1 if x<0 else (0 if x==0 else +1))


def hopfield_train(states, prune=None):
    """Create hopfield network
       states -- list of discretized samples
       prune -- pruning threshold. If None no pruning is performed
       returns symmetric weight matrix with zero-diagonal
    """
    n = len(states[0])
    W = zeros((n,n))
    for state in states:
        v = signum(state)
        W = W + outer(v,v)
    W[diag_indices(n)] = 0
    W = W/len(states)
    if prune: W[abs(W)<prune] = 0
    return W


def hopfield_ask(W, states, n=30):
    """Recall of patterns from Hopfield network
       W -- weight matrix created by hopfield_train()
       states -- list of discretized samples
       n -- number of recall steps to perform
       returns states after recall
    """

    states = signum(states)
    for _ in xrange(n):
        new_states = signum(dot(states,W))
        states = new_states
    return states


def hopfield_energy(W, samples):
    """Energy function of the Hopfield network.
       W -- weight matrix created by hopfield_train()
       samples -- expression matrix with sample vectors
       returns vector with energies of sample vectors
    """
    return array([-0.5*dot(dot(signum(s).T,W),signum(s)) for s in samples])


"PLOTS-----------------------------------------------------------------------------------------------------------------"
def plot_relaxation(data, n=10, prune=None):
    """Plot the relaxation of the state matrix over n steps.
       data -- tuple (labels, genes, samples) from load_data()
       n -- number of relaxation steps
       prune -- pruning threshold. If None no pruning is performed
    """
    labels, genes, samples = data
    W = hopfield_train(samples, prune=prune)
    figure()
    subplot(1, n+1, 1)
    imshow(samples, origin="lower", interpolation="nearest", cmap=cm.RdYlGn_r)
    axis('off')
    states = signum(samples)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #defines empty list
    relaxation=[]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in xrange(n):
        subplot(1, n+1, i+2)
        axis('off')
        imshow(states, origin="lower", interpolation="nearest", cmap=cm.RdYlGn_r)
        new_states = signum(dot(states,W))
        states = new_states
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #generates output files for each state matrix and a concatenated matrix:
        statefilename =  specieslabel + '_S' + str(i+1) + '.st'
        statematrix = states.transpose()
        statematrix = pd.DataFrame(statematrix , index=genenames, columns=None)
        statematrix.to_csv(statefilename, sep="\t")
        relaxation.append(states)
        relaxed_states = np.asarray(relaxation)
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #finding relaxed states
    #print
    #print "--------------RELAXATION INFORMATION--------------"
    #print
    #print "relaxed_states is ", type(relaxed_states)
    #print "relaxed_states shape is", relaxed_states.shape
    #print "stage number = ", relaxed_states.shape[1]
    #print "gene number = ", relaxed_states.shape[2]
    #print "state matrices (relaxation steps) =", relaxed_states.shape[0]
    #print

    #computes the new dimension (shape) components for reshaping
    dim1 = (relaxed_states.shape[0])*(relaxed_states.shape[1])
    dim2 = relaxed_states.shape[2]

    #generates a two-dimensional array(matrix)
    relaxed_states = np.reshape(relaxed_states,newshape=(dim1,dim2))
    #print "relaxed_states matrix contains all the state matrices."
    #print "matrix dimension = ", relaxed_states.shape
    #print

    #generates data frame to format an output file
    relaxed_states = pd.DataFrame(relaxed_states, index=None, columns=genenames)

    #creates a .st output file with all the relaxation matrices
    relaxedname = specieslabel + '_allstates' + '.st'
    relaxed_states.to_csv(relaxedname, sep="\t")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for relaxed states
    image_1 = '1_relaxation_' + specieslabel + '.png'
    plt.savefig(image_1, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def density(W):
    """Compute fraction of non-zero entries in weigth matrix.
       W -- weight matrix created by hopfield_train()
       returns density
    """
    r,c = W.shape
    return 1.0*count_nonzero(W)/r/c


def plot_weight_matrix(data, prune=None, bin=True):
    """Plot the weight matrix of the network.
       data -- tuple (labels, genes, samples) from load_data()
       prune -- pruning threshold. If None no pruning is performed
       bin -- plot binarized weights if true
    """
    from numpy import count_nonzero, array
    #print "Plotting weight matrix ..."
    labels, genes, samples = data
    W = hopfield_train(samples, prune=prune)
    if bin: W = array([signum(row) for row in W])
    w = max(abs(W.min()),abs(W.max()))
    figure()
    imshow(W, vmin=-w, vmax=+w, origin="upper", interpolation="nearest", cmap=cm.RdYlBu)
    axis('off')
    colorbar()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for weight matrix
    image_2 = '2    _weight_matrix_' + specieslabel + '.png'
    plt.savefig(image_2, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_pruning(data, alpha=0.2):
    """Plot adjusted rand index (ARI), estimated rand index (ERI) and
       density of the weight matrix for different pruning thresholds,
       and select the best threshold.
       data -- tuple (labels, genes, samples) from load_data()
       alpha -- allowed decrease of ERI
       returns pruning threshold
    """
    from numpy import nditer
    from sklearn import metrics
    #print "Pruning network ..."
    labels, genes, samples = data
    W = hopfield_train(samples, prune=None)
    labels_t = class2idx(labels)
    labels_0 = states2clusterids(hopfield_ask(W, samples))
    ws = [abs(w) for w in nditer(W)]
    ts = sorted(list(set(ws)))
    aris, eris, ds = [],[],[]
    n_clusters = len(set(labels_0))
    for t in ts:
        W[abs(W)<t] = 0
        labels_pred = states2clusterids(hopfield_ask(W, samples))
        eri = metrics.adjusted_rand_score(labels_0, labels_pred)
        ari = metrics.adjusted_rand_score(labels_t, labels_pred)
        eris.append(eri)
        aris.append(ari)
        ds.append(density(W))
    tes = [(t,s) for t,s in zip(ts,eris) if s > 1-alpha]
    t = max(tes, key=lambda (t,s):t)[0]
    figure()
    plot(ts,aris,'--', color='r')
    plot(ts,eris,'+-', color='b')
    plot(ts,ds,'.-', color='g')
    vlines(t, -0.05, 1.05, colors='k', linestyles='dotted')
    xlabel('Threshold')
    ylabel('TRI, ERI, Density')
    legend(['TRI','ERI','Density'])
    gca().set_ylim(ymin=-0.05, ymax=1.05)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for pruning threesholds
    image_3 = '3_pruning_threesholds_' + specieslabel + '.png'
    plt.savefig(image_3, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return t


def scatter_plot(points, cidx, marker='o', size=50, alpha=1.0):
    """2D scatter plot
        points -- matrix with point coordinates
        cidx -- class label (numeric) for each point
        marker -- point markers
        size -- diameter of point
        alpha -- transparency
    """
    for i in xrange(len(set(cidx))):
        scatter(points[cidx==i,0], points[cidx==i,1], color=idx2color(i),
                marker=marker, s=size, alpha=alpha)


def scatter_plot3d(ax, points, energies, cidx, marker='o', size=50, alpha=1.0):
    """3D scatter plot
       ax -- 3D axis
       points -- matrix with point coordinates
       energies -- Hopfield energies of points
       cidx -- class label (numeric) for each point
       marker -- point markers
       size -- diameter of point
       alpha -- transparency
    """
    for i in xrange(len(set(cidx))):
        ax.scatter(points[cidx==i,0], points[cidx==i,1], energies[cidx==i],
                   c=idx2color(i), marker=marker, s=size, alpha=alpha)


def create_mesh_points(model2D, W, res):
    """Create high-dimensional mesh points
       model2D -- PCA model
       W -- weight matrix of Hopfield network
       res -- resolution of mesh
       returns high-dimensional mesh points
    """
    w = max(W.min(),W.max())*2.0
    x = y = linspace(-w, +w, res)
    X, Y = meshgrid(x, y)
    points2d = array([[x, y] for x, y in zip(ravel(X), ravel(Y))])
    return model2D.inverse_transform(points2d)


def plot_landscape(data, res=50, prune=None):
    """Plot the energy landscape of the Hopfield network in different ways.
       data -- tuple (labels, genes, samples) from load_data()
       res -- resolution of mesh
       prune -- pruning threshold. If None no pruning is performed
    """
    from sklearn.decomposition import PCA
    #print "Plotting landscape ..."
    labels, genes, samples = data
    cidx = class2idx(labels)

    model2D = PCA(n_components=2)
    model2D.fit(samples)

    W = hopfield_train(samples, prune)
    states = samples.copy()
    trajectories2d = [model2D.transform(samples)]
    trajectories = [states]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #this list will storage the states:
    track=[]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in xrange(10):
        states = hopfield_ask(W, states, n=1)
        trajectories.append(states)
        track.append(states)
        newtrack = np.asarray(track)
        trajectories2d.append(model2D.transform(states))
    samples_a = trajectories[-1]
    samples_a_2d = trajectories2d[-1]
    samples2d = trajectories2d[0]
    trajectories = vstack(trajectories)
    trajectories2d = hstack(trajectories2d)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #print
    #print "------------LANDSCAPE INFORMATION--------------"
    #print
    #print "newtrack", type(newtrack), newtrack.shape
    #print "dim1(newtrack) = ", newtrack.shape[0]
    #print "dim2(newtrack) = ", newtrack.shape[1]
    #print "dim3(newtrack) = ", newtrack.shape[2]
    #print

    #computes the new dimension (shape) components for reshaping
    component1 = (newtrack.shape[0])*(newtrack.shape[1])
    component2 = newtrack.shape[2]
    #print "component1 is ", type(component1)
    #print "component2 is ", type(component2)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    newtrack = np.reshape(newtrack,newshape=(component1,component2))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #print
    #print "Data types for output file (samples):"

    #shows objects' names and dimensions
    #print "samples_a", type(samples_a), samples_a.shape
    #print "samples_a_2d", type(samples_a_2d), samples_a_2d.shape
    #print "samples2d", type(samples2d), samples2d.shape
    #print "trajectories", type(trajectories), trajectories.shape
    #print "trajectories2d", type(trajectories2d), trajectories2d.shape
    #print "track", type(track)
    #print "newtrack", type(newtrack), newtrack.shape
    #print "-----------------------------------------------"
    #print

    #generates data frames to build output tables
    output_samples_a =  pd.DataFrame(samples_a, index=None, columns=None)
    output_samples_a_2d =  pd.DataFrame(samples_a_2d, index=None, columns=None)
    output_samples2d =  pd.DataFrame(samples2d, index=None, columns=None)
    output_trajectories = pd.DataFrame(trajectories, index=None, columns=None)
    output_trajectories2d = pd.DataFrame(trajectories2d, index=None, columns=None)
    output_track = pd.DataFrame(newtrack, index=None, columns=genenames)


    #generates output files in .trk format
    output_samples_a.to_csv(specieslabel + '_samples_a.trk', sep="\t")
    output_samples_a_2d.to_csv(specieslabel + '_samples_a_2d.trk', sep="\t")
    output_samples2d.to_csv(specieslabel + '_samples_2d.trk', sep="\t")
    output_trajectories.to_csv(specieslabel + '_trajectories.trk', sep="\t")
    output_trajectories2d.to_csv(specieslabel + '_trajectories2d.trk', sep="\t")
    output_track.to_csv(specieslabel + '_output_track.trk', sep="\t")
    #print
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mesh = create_mesh_points(model2D, W, res)
    points = vstack((mesh,trajectories))
    energies = hopfield_energy(W, points)
    points = model2D.transform(points)
    w = max(points.min(),points.max())*1.5
    x,y = mgrid[-w:w:complex(0,res),-w:w:complex(0,res)]
    from scipy.interpolate import griddata
    z = griddata(points, energies, (x,y), method='nearest')

    from mpl_toolkits.mplot3d import Axes3D
    ax = figure().gca(projection='3d')
    ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, alpha=0.9)
    ax.set_xlabel('1st pc'); ax.set_ylabel('2nd pc'); ax.set_zlabel('Energy')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for energy(pseudo-potential) landscape
    image_4 = '4_energy_landscape_' + specieslabel + '.png'
    plt.savefig(image_4, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    ax = figure().gca(projection='3d')
    ax.plot_wireframe(x,y,z, linestyles='solid', colors=(0.5,0.5,0.5,0.5), alpha=0.5)
    ax.set_xlabel('1st pc'); ax.set_ylabel('2nd pc'); ax.set_zlabel('Energy')
    samples_e = hopfield_energy(W, samples)
    scatter_plot3d(ax, samples2d, samples_e, cidx)
    samples_a_e = hopfield_energy(W, samples_a)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for energy(pseudo-potential) landscape
    image_5 = '5_PCA_landscape_' + specieslabel + '.png'
    plt.savefig(image_5, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    fig = figure()
    imshow(z.T, extent=(-w,w,-w,w), origin='lower', cmap=cm.binary)
    xlabel("1st pc"); ylabel("2nd pc")
    scatter_plot(samples2d, cidx, marker='o', alpha=0.5)
    legend(set(labels), prop={'size':10})
    for i,t in zip(cidx,trajectories2d):
        xt,yt = t[::2],t[1::2]
        plot(xt,yt,'-'+idx2color(i), alpha=0.5)
    plot(samples_a_2d[:,0], samples_a_2d[:,1], 'og', markersize=10, alpha=0.5)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generates image for energy(pseudo-potential) landscape
    image_6 = '6_PCA_contour_plot_' + specieslabel + '.png'
    plt.savefig(image_6, format='png')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~ENDOFPLOTS--------------------------------------------------------------------


def print_usage():
    """#print usage info"""
    print "hclust.py <filename> [-p -n -f]"
    print "-p : enable pruning"
    print "-n : enable normalization"
    print "-f : enable feature selection"


def main(args):
    """Load data and create all plots.
       args -- command line arguments
    """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #print
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #print "running..."
    data = load_data(args[1], '-n' in args)
    if '-f' in args: data = feature_selection(data)
    t = plot_pruning(data) if '-p' in args else None
    plot_relaxation(data, prune = t, n=8)
    plot_weight_matrix(data, prune = t, bin=False)
    plot_landscape(data, prune = t)
    #print
    #print "Finished."
    #print "Run time =", time.clock()
    #print
    #show()



if __name__ == "__main__":
    import sys
    args = sys.argv
    if len(args) < 2:
        print_usage()
    else:
        main(args)
