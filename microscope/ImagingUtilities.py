import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as optimize
from sklearn.neighbors import NearestNeighbors
from skimage import filters, draw, exposure, feature, util, transform, io, color
from scipy.special import jn

#from astropy.modeling.functional_models import AiryDisk2D
tolerance = 1
def gaussian2D(height, center_x, center_y, width_x, width_y,offset):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)+offset
def homemadeAiry(height,center_x,center_y,radius,c):
    """Homemade airy function."""
    Airyparam = np.pi*1.2196
    center_x = float(center_x)
    center_y = float(center_y)
    return lambda x,y: height* ((2*jn(1,np.sqrt(((x-center_x)**2+(y-center_y)**2))*Airyparam/radius)/(0.0001+np.sqrt(((x-center_x)**2+(y-center_y)**2))*Airyparam/radius)))**2+c
def homemadeAirydecoupled(height,center_x,center_y,radiusx,radiusy,c):
    """Homemade airy function. #TODO Pls check!!"""
    Airyparam = np.pi*1.2196
    center_x = float(center_x)
    center_y = float(center_y)
    radiusx = radiusx/Airyparam
    radiusy = radiusy/Airyparam
    return lambda x,y: height* ((2*jn(1,np.sqrt((((x-center_x)/radiusx)**2+((y-center_y)/radiusy)**2))))/(np.sqrt((((x-center_x)/radiusx)**2+((y-center_y)/radiusy)**2))))**2+c


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    data =data.astype(float)
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    c = data.min()
    return height, x, y, width_x, width_y,c

def fit2D(data,model):
    """Returns (height, x, y, width_x, width_y)
    the model parameters of a 2D distribution found by a fit
    Model must take parameters (height, x, y, width_x, width_y). 
    """
    #print(type(data))
    params = moments(data)
    length = len(data[:,0])
    errorfunction = lambda p:  np.ravel(model(*p)(*np.indices(data.shape)) -
                                 data)
    
    if model == homemadeAiry: # or model == AiryDisk2D
        params = list(params)
        params.pop(4)
        tuple(params)
        #TODO sort out bounds being sensible/responsive to windowsize
        bounds = ([0,0,0,0,0],[512,length,length,length,100])
    else: 
        bounds = ([0,0,0,0,0,0],[512,length,length,length,length,100])
        #print("starting params",params)
    #p, pcov,infodict, mesg, ier = optimize.least_squares(errorfunction, params,bounds = ([0,0,0,0,0],[512,60,60,60,60,100]))#,full_output=1,)
    #print(errorfunction(params))
    #print(params)
    result = optimize.least_squares(errorfunction, params,bounds  = bounds,max_nfev = 5)#,full_output=1,)
    try:
        J = result['jac']
        cov = np.linalg.inv((J.T).dot(J))
        errs = np.sqrt(np.diagonal(cov))
    except:
        errs = [100,100,100,100,100,100]
    p = result.x

    #TODO: Reilaible error estimate 
    return p, errs
def quickfit2D(data, model):
    errorfunction = lambda p:  np.ravel(model(*p)(*np.indices(data.shape)) -
                                 data)
    params = (100,data.shape[0]/2,data.shape[1]/2)
    bounds = ([0,0,0],[512,data.shape[0],data.shape[1]])
    length = len(data[:,0])
    result = optimize.least_squares(errorfunction, params,bounds  = bounds,max_nfev = 5)#,full_output=1,)
    try:
        J = result['jac']
        cov = np.linalg.inv((J.T).dot(J))
        errs = np.sqrt(np.diagonal(cov))
    except:
        errs = [100,100,100]
    p = result.x

    #TODO: Reilaible error estimate 
    return p, errs


def preProcess(filepath):
    """ Takes a tiff image path and ensures it is in the correct grayscale format and is normalised sensibly. 
    """
    if filepath.endswith(".npy"):
        image = np.load(filepath) 
    else:
        image = io.imread(filepath) 
        image = color.rgb2gray(color.rgba2rgb(image))
    return image

def BlobDetection(image,scaling = 10):
    """ Function to detect Blobs. Takes an image and scaling arg, which controls the sensitiviy of blob detector.
        Performs a 5x5 px blur first to help with blob detection
        Returns an open cv keypoints object. 
    """
    image = exposure.rescale_intensity(image,out_range = (0,scaling)) #NB this sets the sensitiviy of the blob detection!
    blobs = feature.blob_log(image,min_sigma=7,num_sigma = 3)

    windowsize = 50
    # Mesh for taking a radial slice
    meshx,meshy = np.meshgrid(np.arange(windowsize*2),np.arange(windowsize*2))
    R = np.sqrt((meshx-windowsize)**2+(meshy-windowsize)**2)
    nanoparticles = []
    otherblobs = []
    for blob in blobs:
            y, x, r = blob
            if  image[int(y),int(x)]>0.1: #filter weak blobs
                blobimage = image[int(y)-windowsize:int(y)+windowsize,int(x)-windowsize:int(x)+windowsize]
                if blobimage.shape == (windowsize*2,windowsize*2): #filter blobs on edge of image
                    radial = blobimage[(R>= 15-0.5)&(R<= 15+0.5)]
                    std = radial.std()
                    midpoint = blobimage[windowsize,windowsize]
                    diff = (midpoint-radial.max())/midpoint
                    print(x,y)
                    print(midpoint,np.max(radial),diff)
                    if diff>0.7: #Filter oblong blobs TODO this seems very image dependant...
                        nanoparticles.append((x,y))
                    else:
                        otherblobs.append((x,y))
    print(nanoparticles)
    nanoparticles = np.array(nanoparticles)
    otherblobs = np.array(otherblobs)
    fig,axes = plt.subplots(1,1)
    try:
        axes.scatter(nanoparticles[:,0],nanoparticles[:,1],color = "red",s = 2)
    except:
        pass
    axes.scatter(otherblobs[:,0],otherblobs[:,1],color = "cyan",s = 2)
    axes.imshow(image)
    plt.show()
    return nanoparticles



def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    squareedge = int((data.shape[0])/2)-10
    return np.array(radialprofile)[:squareedge], np.linspace(0,squareedge,squareedge)+0.001 
def AiryDisk1D(r,I_0,radius,c):
    """
        Airy disk function, in format required by curve_fit
    """
    Airyparam = np.pi*1.2196
    return I_0*(2*jn(1,r*Airyparam/radius)/(r*Airyparam/radius))**2+c
def fit1D(data,model):
    xpoints = np.linspace(0,len(data),len(data))+0.0001
    fit, pcov  = optimize.curve_fit(model,xpoints,data,p0 = (100,170,10))
    errs = np.sqrt(np.diag(pcov))
    return fit,errs
def GetMagnification(loc1,loc2,image1,image2):
    pass #TODO use this as a wrapper 
def match_row(array,row):
    """
        Match rows of a 2D array based on norm. 
        Returns index of best match and a bool to indicate not being able to find a row.
        Assumes that nearest neighbours are similar in both images, i.e. not too many particles have moved out of FoV  
    """
    best_match = 0
    best_difference = 100000
    sucess = True
    for i in range(len(array[:,0])):
        #Loop over rows to find best match
        difference = np.sum((array[i,:]-row)**2)
        if difference < best_difference:
            best_difference = difference
            best_match = i
    if best_difference> 100: #arbitrairy cutoff for poor fit
        print("NO MATCH FOUND")
        sucess = False
    return int(best_match), sucess 
def sort_locations(loc1,loc2):
    """
        Takes 2 sets of locations and sorts them into nice format for affine magnifier using the scikit-lean kNN implementation
        This also sorts the issue of different numbers of features being detected.
        Will return 2 lists, possibly empty if no good matches are found...
    """
    distancelist = [] #List for both arrays
    for locations in [loc1,loc2]:
        nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(np.array(locations))
        distances, indices = nbrs.kneighbors(np.array(locations))
        distancelist.append(distances)
    newloc1 = []
    newloc2 = []
    for i in range(len(loc1)):
        index, sucess = match_row(distancelist[1],distancelist[0][i,:])
        if sucess:
            newloc1.append(loc1[i])
            newloc2.append(loc2[index])
    return newloc1, newloc2
def combine_windows(image,nanoparticles,windowsize):
    """
        Takes img and locations, 
 
        Extracts windows around nanoparticles wholly contained within field of view
        windowsize = square box
    """
    oversamplesdshape = (windowsize*10,windowsize*10)

    total = np.zeros(oversamplesdshape)

    counter = 0
    for point in nanoparticles:
        print(point)
        
        blobimage = image[int(point[1])-windowsize:int(point[1])+windowsize,int(point[0])-windowsize:int(point[0])+windowsize]
        print(blobimage.shape)
        blobimage = transform.resize(blobimage,oversamplesdshape,order = 0)
        fit,err = fit2D(blobimage,gaussian2D)
        print("Fit, errors", fit,err)
        c_x, c_y = fit[1], fit[2]
        # plt.matshow(blobimage)
        # plt.show()
        matrix = np.float32([[1,0,(oversamplesdshape[0]/2 -c_x)],[0,1,(oversamplesdshape[1]/2-c_y)]])
        tfrom = transform.AffineTransform(translation=((oversamplesdshape[0]/2 -c_x),(oversamplesdshape[1]/2-c_y)))
        shifted = transform.warp(blobimage, tfrom.inverse,mode = "edge",order = 0)
        # plt.matshow(shifted)
        total += shifted
        # plt.show(0)
        counter +=1
    if counter == 0:
        counter = 1
        print("Warning: no particles accepted, outputting ones")
        total = np.ones(oversamplesdshape)
    else:
        total = (total/counter)
        print("Average of ",counter,"Windows")
 
    return total
def read_params(folderpath):
    """
        Reads parameters from a readme file in comma seperated format
        image, param \n 
        Ignores lines which cant be read
    """
    f = open(folderpath+"readme.txt", 'r')
    x = f.readlines()
    f.close()
    imagelist, paramlist = [], []
    for line in x:
        line = line.strip()
        seperated = line.split(",")
        try:
            imagelist.append(int(seperated[0]))
            paramlist.append(float(seperated[1]))
        except:
            pass

    imagelist = np.array(imagelist)
    paramlist = np.array(paramlist)
    return imagelist, paramlist
