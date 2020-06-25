import numpy as np
import microscope.ImagingUtilities as imu
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import convolve
from skimage import restoration, filters, draw, exposure, feature, util, transform
import time
from scipy.spatial import Voronoi
from scipy.signal import find_peaks
import os
import time
from sklearn.metrics import pairwise_distances
import scipy.optimize as optimize
import cProfile , pstats

"""
Code to generate a sample lattice image

Inputs: Psf radius; lattice points

"""


def generate_lattice(image_shape, lattice_vectors, offset) :
    """
    Generate points from a lattice which lie withing a given shape.
    
    https://stackoverflow.com/questions/6141955/efficiently-generate-a-lattice-of-points-in-python
    NB: this generates a lattice with as many points within the shape as possible, meaning that different lattices have different numbers of sites for the same shape...
    """
    center_pix = np.array(image_shape) // 2
    # Get the lower limit on the cell size.
    dx_cell = max(abs(lattice_vectors[0,0]), abs(lattice_vectors[1,0]))
    dy_cell = max(abs(lattice_vectors[0,1]), abs(lattice_vectors[1,1]))
    # Get an over estimate of how many cells across and up.
    nx = image_shape[0]//dx_cell
    ny = image_shape[1]//dy_cell
    # Generate a square lattice, with too many points.
    # Here I generate a factor of 4 more points than I need, which ensures 
    # coverage for highly sheared lattices.  If your lattice is not highly
    # sheared, than you can generate fewer points.
    x_sq = np.arange(-nx, nx, dtype=float)
    y_sq = np.arange(-ny, nx, dtype=float)
    x_sq.shape = x_sq.shape + (1,)
    y_sq.shape = (1,) + y_sq.shape
    # Now shear the whole thing using the lattice vectors
    x_lattice = lattice_vectors[0,0]*x_sq + lattice_vectors[1,0]*y_sq +offset[0] # off set would come here???
    y_lattice = lattice_vectors[0,1]*x_sq + lattice_vectors[1,1]*y_sq +offset[1]
    # Trim to fit in box.
    mask = ((x_lattice < image_shape[0]/2.0)
             & (x_lattice > -image_shape[0]/2.0))
    mask = mask & ((y_lattice < image_shape[1]/2.0)
                    & (y_lattice > -image_shape[1]/2.0))
    x_lattice = x_lattice[mask]
    y_lattice = y_lattice[mask]
    # Translate to the centre pix.
    x_lattice += center_pix[0]
    y_lattice += center_pix[1]
    # Make output compatible with original version.
    out = np.empty((len(x_lattice), 2), dtype=float)
    out[:, 0] = y_lattice
    out[:, 1] = x_lattice
    return out

def populate_atoms_convolve(latticelist,imageshape,psf_radius,density = 0.3, savefile = "",layers = [1.0],psf_noise = 0.2,exposure_scale = 1,save = True,seed =11546546 ):
    """
        Create an image via a convolution. First populates latttice sites with atoms, by taking samples from a binomial distribution, then adding variation to the atom signal by drawing from a log guassian distribution. After convolution with the psf additive noise is applied following the formula in the Andor manual.
    """
    background = 1000*exposure_scale
    signal = 2000*exposure_scale
    N_r = 10
    gain_factor = 10
    sensitivity = 4.46
    rescale_factor = 5

    image = np.zeros((imageshape[0]*rescale_factor,imageshape[1]*rescale_factor))
    trutharray = np.zeros(len(latticelist))
    j = 0
    for layer in layers:
        targetlayer = False
        if layer == 1.0:
            truthlist = []
            targetlayer = True

        #Generate "Atoms" PseudoRandomly, with noisy binomial distribution
        np.random.seed(seed+j)
        print("layer,seed",layer,seed+j) 
        j+=1
        heights= (layer*np.random.binomial(1,density,len(latticelist)))
        print(heights)
        
        for i,point in enumerate(latticelist):
            r = np.sqrt((point[0]-imageshape[0]/2)**2+(point[1]-imageshape[1]/2)**2)
            if r>(imageshape[0]*0.45): #radial logic
                height = 0# layer*np.random.binomial(1,0.9) +np.random.normal(loc = 0,scale = 0.1)
            else:
                height = heights[i]
                #print(height)
            if height >0: # to allow hights between 0 or 1
                height *=  (signal)*np.random.lognormal(mean = 0,sigma = psf_noise) # Add noise here to remove ambiguity
                image[int(np.rint(point[1]*rescale_factor)),int(np.rint(point[0]*rescale_factor))] += height #TODO speed up by using psf array...
            if targetlayer: 
                truthlist.append(height !=0) # spaces with no psf will be here
                trutharray[i] += int(height !=0)

    image = convolve(image,airy_psf(30*rescale_factor,psf_radius*rescale_factor)/255,mode = "same")
    image = transform.resize(image,imageshape,order=0)
    # fig,axes = plt.subplots(1,1)
    # axes.imshow(image)
    # axes.scatter(latticelist[:,0],latticelist[:,1],s = 1,color = "red")
    # plt.show()
    image += background
    # plt.matshow(image)
    # plt.show()
    image += np.random.normal(loc = 0,size = imageshape)*np.sqrt(2*(gain_factor/sensitivity)*image+N_r**2) #Noise... 
    # plt.matshow(image)
    # plt.show()
    # image = exposure.rescale_intensity(image, out_range = (0,255))
    truthandlattice = np.column_stack((np.array(latticelist),trutharray))
    if save:
        np.save("AtomTestImages/"+savefile+"data.npy",image)
        #np.save("AtomTestImages/"+savefile+"truth.npy",groundtruth)
        #np.save("AtomTestImages/"+savefile+"truthlist.npy",np.array(truthlist))
        np.save("AtomTestImages/"+savefile+"truthandlattice.npy",truthandlattice)
    return np.array(image), np.array(truthandlattice)

def populate_atoms_convolve_plot(latticelist,imageshape,psf_radius,density = 0.3, savefile = "",layers = [1.0],psf_noise = 0.2,exposure_scale = 1,save = True,seed =11546546 ):
    """
        Create an image
    """
    background = 1000*exposure_scale
    signal = 2000*exposure_scale
    N_r = 10
    gain_factor = 10
    sensitivity = 4.46
    rescale_factor = 5

    image = np.zeros((imageshape[0]*rescale_factor,imageshape[1]*rescale_factor))
    trutharray = np.zeros(len(latticelist))
    j = 0
    for layer in layers:
        targetlayer = False
        if layer == 1.0:
            truthlist = []
            targetlayer = True

        #Generate "Atoms" PseudoRandomly, with noisy binomial distribution
        np.random.seed(seed+j)
        print("layer,seed",layer,seed+j) 
        j+=1
        heights= (layer*np.random.binomial(1,density,len(latticelist)))
        print(heights)
        
        for i,point in enumerate(latticelist):
            # r = np.sqrt((point[0]-imageshape[0]/2)**2+(point[1]-imageshape[1]/2)**2)
            # if r>(imageshape[0]*0.45): #radial logic
            #     height = 0# layer*np.random.binomial(1,0.9) +np.random.normal(loc = 0,scale = 0.1)
            # else:
            height = heights[i]
                #print(height)
            if height >0: # to allow hights between 0 or 1
                height *=  (signal)*np.random.lognormal(mean = 0,sigma = psf_noise) # Add noise here to remove ambiguity
                heights[i] = height
                image[int(np.rint(point[1]*rescale_factor)),int(np.rint(point[0]*rescale_factor))] += height #TODO speed up by using psf array...
            if targetlayer: 
                truthlist.append(height !=0) # spaces with no psf will be here
                trutharray[i] += int(height !=0)

    image = convolve(image,airy_psf(30*rescale_factor,psf_radius*rescale_factor)/255,mode = "same")
    image = transform.resize(image,imageshape,order=0)
    fig,axes = plt.subplots(1,1)
    # axes.imshow(image)
    axes.scatter(latticelist[:,0],latticelist[:,1],c = heights  ,cmap = "inferno",marker="s")
    
    axes.set_xlim([300,400])
    axes.set_ylim([200,300])
    axes.set_aspect("equal")
    plt.figure()
    plt.hist(heights,bins = 100,fc = "orange",ec = "purple")
    
    # plt.show()

    image += background
    # plt.matshow(image)
    # plt.show()
    image += np.random.normal(loc = 0,size = imageshape)*np.sqrt(2*(gain_factor/sensitivity)*image+N_r**2) #Noise... 
    # plt.matshow(image)
    # plt.show()
    # image = exposure.rescale_intensity(image, out_range = (0,255))
    truthandlattice = np.column_stack((np.array(latticelist),trutharray))
    if save:
        np.save("AtomTestImages/"+savefile+"data.npy",image)
        #np.save("AtomTestImages/"+savefile+"truth.npy",groundtruth)
        #np.save("AtomTestImages/"+savefile+"truthlist.npy",np.array(truthlist))
        np.save("AtomTestImages/"+savefile+"truthandlattice.npy",truthandlattice)
    fig,axes = plt.subplots(1,1)
    axes.imshow(image,cmap = "inferno")
    # axes.scatter(latticelist[:,0],latticelist[:,1],c = heights/1000,cmap = "inferno",marker="s")
    
    axes.set_xlim([300,400])
    axes.set_ylim([200,300])
    axes.set_aspect("equal")
    # plt.show()
    return np.array(image), np.array(truthandlattice)


#   
def bin_lattice(image, latticesites, latticevectors):
    """
        Histogram an image based on defined lattice. 
        Inputs: image ndarray
                latticesites ndarray
                latticevectors tuple?
        Outputs
            totals list, with same order as latticesites
    """
    #scale = 5
    #image = transform.resize(image,(image.shape[0]*scale,image.shape[1]*scale),order = 0,mode = "edge")
    
    halfvectors =latticevectors/2

    totals = []
    #print(halfvectors)
    firstpoint = latticesites[0]
    pixelsize = 1
    i = 0
    for i in range(latticesites.shape[0]):
        point = latticesites[i,:]
        square = np.array([[point[1]-(halfvectors[0][1]+halfvectors[1][1]),point[0]-(halfvectors[0][0]+halfvectors[1][0]),],[point[1]-(-halfvectors[0][1]+halfvectors[1][1]),point[0]-(-halfvectors[0][0]+halfvectors[1][0])],[point[1]-(-halfvectors[0][1]-halfvectors[1][1]),point[0]-(-halfvectors[0][0]-halfvectors[1][0])],[point[1]-(halfvectors[0][1]-halfvectors[1][1]),point[0]-(halfvectors[0][0]-halfvectors[1][0])]])
        # if np.all(square>0) and np.all(square<image.shape[0]):
            #NB the way the lattice is returned means point[0] is 2nd index of the mask
        mask = draw.polygon2mask(image.shape, square)

        totals.append(np.sum(image[mask]))
    
    return np.array(totals)



def find_blobs(image,Plot = False):
    """
        Detects blobs in an image file using blob_log, then filters according to height and circularity via a guassian fit. 

        Plots the remaining blobs for a visual inspection. 
    """
    if isinstance(image,str):
        image = np.load(image)
    windowsize = 7
    print("looking for blobs")
    blobs = feature.blob_log(image,min_sigma = 3, max_sigma= 6, num_sigma=3,threshold=0.9,overlap = 0.5)
    print("blobs found")

    atoms = []
    meshx,meshy = np.meshgrid(np.arange(windowsize*2),np.arange(windowsize*2))
    #print(meshx)
    #print(meshy)
    R = np.sqrt((meshx-windowsize)**2+(meshy-windowsize)**2)
    for blob in blobs:
        y, x, r = blob
        if  image[int(y),int(x)]>100: #filter weak blobs
            blobimage = image[int(y)-windowsize:int(y)+windowsize,int(x)-windowsize:int(x)+windowsize]
            if blobimage.shape == (windowsize*2,windowsize*2): #filter blobs on edge of image
                #try: #mostlikely a value error on fit2D 
                try:
                        
                        # R  = np.sqrt((meshx-x)**2+(meshy-y)**2)
                        radial = blobimage[(R>= 4-0.5)&(R>= 4+0.5)] #TODO check this!! seems to be wrong logic
                        std = radial.std()
                        diff = (np.max(blobimage)-radial.max())/np.max(blobimage)
                        #print(std,diff)
                        #print(x,y)
                        #print(fit[1],fit[2])
                        if diff>0.5: #Filter oblong blobs TODO this seems very image dependant...
                            atoms.append((x,y))
                            # fit, err = imu.fit2D(blobimage,imu.gaussian2D)
                            
                            # if abs(fit[3]<windowsize*0.7):
                                # print("blob plotting")
                                # c = plt.Circle((x,y),fit[3],color  = "red",linewidth = "2", fill = False)
                                # axes.add_patch(c)
                                #print(x,y)
                                #print(fit[1],fit[2])
                                # atoms.append((x-(windowsize-fit[1]),y-(windowsize-fit[2])))
                                #print(x,y)
                except:
                    #print('Fail on blob fitting')
                    pass
    atoms = np.array(atoms)
    if Plot:
        fig, axes = plt.subplots(1,1)
        axes.imshow(image)
        axes.scatter(atoms[:,0],atoms[:,1],s = 2,color = 'red')

    return atoms

def find_lattice_vectors_new(atoms,min_angle,max_angle,steps,spacing_bounds = (3,5),Plot = False):
    """ Search angles while allowing for some variation in lattice spacing over set bounds
        For reliablity use a bruteforce search method.    
    """ 
    # pr = cProfile.Profile()
    # pr.enable()

    qualities = []
    firstangles = np.linspace(min_angle,max_angle,steps)
    for angle in firstangles:
        qualities.append(angle_quality_free_spacing(angle,atoms))
    qualities = np.array(qualities)
    best1 = firstangles[np.argmax(qualities)]
    spacing1 = angle_quality_free_spacing(best1,atoms,ReturnSpacing=True)

    #NB assume approximately orthogonal to speed up
    secondangles = np.linspace((best1+np.pi/2-0.1),(best1+np.pi/2+0.1),int(steps/10))
    secondqualities = []
    for angle in secondangles:
        secondqualities.append(angle_quality_free_spacing(angle,atoms))
    best2 = secondangles[np.argmax(np.array(secondqualities))]
    spacing2 = angle_quality_free_spacing(best2,atoms,ReturnSpacing=True)

    if Plot:
        plt.figure()
        plt.plot(firstangles,qualities,color = "purple")
        plt.xlabel("Angle/Radian")
        plt.ylabel("Integrated FFT peak")
        # plt.axvline(best1,color = "red",ls = "--")
        #plt.show()
    latticevectors = np.array(((np.cos(best1)*spacing1,np.sin(best1)*spacing1),(np.cos(best2)*spacing2,np.sin(best2)*spacing2)))
    latticevectors = np.squeeze(latticevectors)
    # pr.disable()
    # pr.print_stats()
    print("best anglees",best1,best2)
    return latticevectors



def angle_quality(angle,points, Plot = False):
    """
    Take a Nx2 array of points and and a lattice angle determine the quality of the guess of lattice angle

    To do this we analyse a histogram of projected mutual distances, as in Wietenburg's thesis. 

    Then we look for a peak in the FFT of this histogram around the spatial frequency of the lattice. 

    The function returns a float corresponding to the integrated FFT PSD near the lattice peak. 

    """
    spacing = 4.655 #TODO this would be a sensible class variable...

    #Get projected mutual distances
    normal = np.array([np.sin(angle), np.cos(angle)])
    projected = np.inner(points,normal.T)
    mutualdistance = pairwise_distances(projected.reshape(-1,1)) 

    #Histogram 
    n_bins = int(18*np.max(mutualdistance)/spacing) #Such that fft contains the signal and a bit more
    binned, edges = np.array(np.histogram(np.ravel(mutualdistance),n_bins))

    # Take FFT of histogram
    spectrum = np.fft.fft(binned)
    freq = np.fft.fftfreq(binned.shape[-1],d = np.max(np.ravel(mutualdistance))/n_bins)

    #Integrate over peak TODO how robust is the hardcoded window size?
    window = np.argwhere(np.logical_and(freq>(1/spacing-0.1),freq<(1/spacing+0.1)))
    
    #Plot histogram and fft for a visual feedback
    if Plot:
        plt.figure("Samplehistogram")
        plt.hist(np.ravel(mutualdistance),bins = n_bins,lw = 2,fc = (np.random.rand(1)[0],0,0,0.5))
        plt.xlim(10,50)
        plt.figure("FFT")
        plt.plot(freq,np.abs(spectrum))
        plt.show()
    return np.abs(spectrum[window]).sum()
def angle_quality_free_spacing(angle,points, Plot = False,ReturnSpacing = False):
    """
    Take a Nx2 array of points and and a lattice angle determine the quality of the guess of lattice angle

    To do this we analyse a histogram of projected mutual distances, as in Wietenburg's thesis. 

    Then we look for a peak in the FFT of this histogram around the spatial frequency of the lattice. 

    The function returns a float corresponding to the integrated FFT PSD near the lattice peak. 

    """
    #Get projected mutual distances
    normal = np.array([np.sin(angle), np.cos(angle)])
    projected = np.inner(points,normal.T)
    mutualdistance = pairwise_distances(projected.reshape(-1,1)) 

    #Histogram 
    n_bins = int(18*np.max(mutualdistance)/4) #Such that fft contains the signal and a bit more
    binned, edges = np.array(np.histogram(np.ravel(mutualdistance),n_bins))

    # Take FFT of histogram
    spectrum = np.abs(np.fft.fft(binned))**2
    freq = np.fft.fftfreq(binned.shape[-1],d = np.max(np.ravel(mutualdistance))/n_bins)

    # Find spacing by looking for largest peaks
    peaks, properties = find_peaks(spectrum,height=np.mean(spectrum))
    #print(freq[peaks])
    #print(properties)
    if len(peaks)==0: #Sometimes no peaks are to be found...
        return 0
    heights = properties['peak_heights']
    sortedpeaks = peaks[np.argsort(heights)]
    
    #print(heights[np.argsort(heights)])
    
    spacing = 1/freq[sortedpeaks[-1]]
    #Integrate over peak TODO how robust is the hardcoded window size?
    window = np.argwhere(np.logical_and(freq>(1/spacing-0.5),freq<(1/spacing+0.5)))
    
    #Plot histogram and fft for a visual feedback
    if Plot:
        plt.figure("Samplehistogram")
        plt.hist(np.ravel(mutualdistance),bins = n_bins,lw = 2,fc = (np.random.rand(1)[0],0,0,0.5))
        plt.xlim(10,50)
        plt.figure("FFT")
        plt.plot(freq,np.abs(spectrum))
        plt.scatter(freq[sortedpeaks[-2:]],spectrum[sortedpeaks[-2:]],color = 'red')
        plt.show()
    if abs(spacing)>10 or abs(spacing)<2: # set a bound on acceptable fft peak positions... This may confuse the optimiser??
        return 0 
    if ReturnSpacing:
        return spacing
    return np.abs(spectrum[window]).sum()

def search_range(anglerange,points,angletruth = None):
    """
        Wrapper for producing a plot to show optimal lattice angle
    """
    plt.figure("Quality vs angle")
    qualities = []
    for angle in anglerange:
        qualities.append(angle_quality(angle,points,False))
    plt.scatter(anglerange,qualities)
    plt.ylabel("FFT peak at 1/ lattice spacing")
    plt.xlabel("Radian lattice angle. Image generated with 0.1")
    if angletruth:
        plt.axvline(angletruth,color="red")
    plt.show()
    return 0 
def rich_lucy_homebrew(image, psf, iterations=50, clip=True,boundary = "reflect"):
    """
        Copy of the scikit image RL algo
    """
    
    print("shape intial",image.shape)
    image = image.astype(np.float)
    # image = np.pad(image,psf.shape[0],mode= "reflect")
    psf = psf.astype(np.float)
    im_deconv = np.full(image.shape, 0.5)
    psf_mirror = psf[::-1, ::-1]

    for _ in range(iterations):
        relative_blur = image / convolve(im_deconv, psf,mode = "same")#,boundary= boundary)
        im_deconv *= convolve(relative_blur, psf_mirror,mode = "same")#,boundary= boundary)

   #  if clip:
   #      im_deconv[im_deconv > 1] = 1
   #      im_deconv[im_deconv < -1] = -1
    # print("shape before",im_deconv.shape)
    # im_deconv = util.crop(im_deconv,psf.shape[0])
    # print("shape after",im_deconv.shape)
    return im_deconv
 



def find_lattice_vectors(atoms, guess_angle):
    """
        Estimate lattice angles and offset via analysis of sparse points

        find points => filter based on aloneness => estimate lattice. 
        TODO fix this documentation

    """

    # atoms = find_blobs(image)
    print("Finished blob")
    
    minfunc = lambda x,*args: (10-angle_quality(x/100,*args))/1e4
    result = optimize.basinhopping(minfunc, 100*(guess_angle+0.01),T = 1,minimizer_kwargs={"method":"Nelder-Mead","args":(atoms,False)},stepsize = 1,niter=100)
    print(result)
    #angle_quality(result["x"]/100,atoms,True)
    orthogonal_result = optimize.basinhopping(minfunc, 100*(result["x"]/100+1.59),T = 1,minimizer_kwargs={"method":"Nelder-Mead","args":(atoms,False)},stepsize = 0.5,niter=10)
    print(orthogonal_result)
    #angle_quality(orthogonal_result["x"]/100,atoms,True)
    angle1,angle2 = result["x"]/100,orthogonal_result["x"]/100
    # angles = (0.78,0.78+np.pi/2)
    spacing = 4.655
    latticevectors = np.array(((np.cos(angle1)*spacing,np.sin(angle1)*spacing),(np.cos(angle2)*spacing,np.sin(angle2)*spacing)))
    latticevectors = np.squeeze(latticevectors)
    # angle = 0.1
    # latticevectors = ((np.cos(angle)*spacing,np.sin(angle)*spacing),(-np.sin(angle)*spacing,np.cos(angle)*spacing))
    
    
   
    # fig,ax = plt.subplots(1,1)
    # ax.imshow(image)
    # # for point in generate_lattice(image.shape,latticevectors,bestoffset):
    # #     c = plt.Circle(point,radius=0.2,color = "red")
    # #     ax.add_patch(c)
    # # plt.show()
    return latticevectors


def find_offset(image,atoms, latticevectors,Plot= False):
    offsetminfunc = lambda offset,*args: offset_quality(offset,*args)
    # offset_result = optimize.basinhopping(offsetminfunc, (0,0),T = 5,minimizer_kwargs={"method":"Nelder-Mead","args":(atoms,latticevectors,image.shape)},stepsize = 1,niter=50)
    # print(offset_result)
    # bestoffset = offset_result["x"]
    offset_result = optimize.minimize(offsetminfunc,(0,0),args = (atoms,latticevectors,image.shape))#,bounds = [[-spacing/2,spacing/2],[-spacing/2,spacing/2]])
    print(offset_result)
    bestoffset = offset_result["x"]
    print(offset_quality((0,0),atoms,latticevectors,image.shape))
    if Plot:
        fig,ax = plt.subplots(1,1)
        ax.imshow(image)
        lattice = generate_lattice(image.shape,latticevectors,bestoffset)
        ax.scatter(lattice[:,0],lattice[:,1],s = 1,color = "red")
    return bestoffset

def offset_quality(offset,atoms,lattice_vectors,shape = (512,512)):
    lattice = generate_lattice(shape,lattice_vectors,offset)
    differences = []
    i = 0
    if len(atoms)>40:
        N = 40
    else:
        N = len(atoms)

    for i in range(N): #TODO make this a random sample
        differences.append(np.sum(np.abs(lattice-np.array(atoms[i]))**2,axis = 1).min())
     #   print(differences[i])
    total_deviation = np.sum(differences)
    #print("total",total_deviation)
    return total_deviation
def find_atoms(image,psf,lattice_vectors,offset,iterations =50,plot = None):

    """ Take an image and using a known psf and lattice, return a boolean list of atom occupations
    Uses deconvolution algos
            Image: ndarray
            psf: ndarray
            lattic_vectors: tuple/array
            offset: tuple
            iterations: int, optional
        outputs
            atoms; ndarray of bools
    """ 
    t1 = time.time()
    lattice_array = generate_lattice(image.shape,lattice_vectors,offset)
    print("time for lattice",time.time()-t1)
    
    t2 = time.time()

    image = np.pad(image,int(psf.shape[0]/2),mode= "reflect")
    deconvolved = restoration.richardson_lucy(image,psf,iterations,False)
    deconvolved = util.crop(deconvolved,int(psf.shape[0]/2))
    
    print("time for devonvolve",time.time()-t2)
    t3 = time.time()
    totals = bin_lattice(deconvolved,lattice_array,lattice_vectors)
    #totals = bin_image(deconvolved,lattice_array,4.655) #WARNING, OLD METHOD
    print("time for binning",time.time()-t3)
    t4 = time.time()
    threshold = filters.threshold_minimum(totals)
    atoms = (np.array(totals)> threshold)
    print("time for threshing", time.time()-t4)
    if plot:
        plt.figure()
        plt.hist(totals, 100)
        plt.axvline(threshold,0,10000,color="red")
        plt.xlabel("Counts")
        plt.ylabel("Occurances")
    atomsandlattice = np.column_stack((lattice_array,atoms))
    return atoms, deconvolved, atomsandlattice

def hist_and_thresh(deconvolved,lattice,lattice_vectors):
    """ Take a deconvolved image and produce histogram..."""
    totals = bin_lattice(deconvolved,lattice,lattice_vectors)
    threshold = filters.threshold_minimum(totals)
    atoms = (np.array(totals)> threshold)
    return totals, threshold, atoms  


def airy_psf(window,psf_radius):
    """ Retrun a numpy array of a psf normalised to max height 255 

    """

    psf  = np.zeros((window,window))
    for i in range(window):
        for j in range(window):
            psf[i,j] = imu.homemadeAiry(255,window/2,window/2,psf_radius,0)(i,j)
            if i == j == int(window/2): #Catch for divide by zero
                psf[i,j] = 255
    return psf
def find_nearest(row,a):
    #print(a-row)
    #print(np.sum(abs(a-row),axis = 1))
    return np.sum(abs(a-row),axis = 1).argmin()

def check_fidelity(truthlist,atomlist,axes = None):
    """
        check fidelity when lattice lists are not nesesarrily the same shape
        assumes thruthlist and atomlist contain data on lattice positions.
    """
    truthlattice = truthlist[:,0:2]
    atomlattice = atomlist[:,0:2] 

    #best = find_nearest(atomlattice[34],truthlattice)
    best_array = np.apply_along_axis(find_nearest,1,atomlattice,truthlattice)

    truthlattice = truthlattice[best_array]
    truthlist = truthlist[best_array]
    error = (np.sum(np.abs(truthlattice-atomlattice),axis = 1))
    matches = 0
    errors = 0
    for i in range(atomlattice.shape[0]):
        if error[i]<1:
            if truthlist[i,2] == atomlist[i,2]:
                matches +=1
            else:
                errors +=1
    
    if axes:
        axes.cla()
        matchpoints = []
        errorpoints = []
        for i in range(atomlattice.shape[0]):
            if error[i]<2:
                if truthlist[i,2] == atomlist[i,2]:
                    matchpoints.append(atomlattice[i])
                else:
                    errorpoints.append(atomlattice[i])
        matchpoints = np.array(matchpoints)
        errorpoints = np.array(errorpoints)
        axes.scatter(matchpoints[:,0],matchpoints[:,1],color = 'white',s = 1)
        try:# this will fail if there are no errorpoints...
            axes.scatter(errorpoints[:,0],errorpoints[:,1],color = 'red',marker = 'x',s = 10)
        except:
            pass
    # print(matches/atomlattice.shape[0])
    # print(atomlattice[errors])
    
    # fig,ax =plt.subplots(1,1)
    # ax.scatter(atomlist[:,0],atomlist[:,1],s = atomlist[:,2]*10+1,alpha = 0.5)
    # ax.scatter(truthlist[:,0],truthlist[:,1],s = truthlist[:,2]*10+1,alpha = 0.5)

    # plt.show()
    # print(np.argwhere(error>1))
    # plt.plot(error)
    # plt.show()
    # print(atomlattice[34])
    #print(1-errors/matches)
    if matches>0:
        fidelity = 1-errors/matches
    else:
        fidelity =0
    #fidelity = matches/atomlattice.shape[0] # This includes errors at edges from misalignment of lattices. Which is a bit unfair?
    return fidelity
# def hist_and_threshKmeans(deconvolved,lattice,lattice_vectors):
#     """ Take a deconvolved image and produce histogram..."""
#     totals = bin_lattice(deconvolved,lattice,lattice_vectors)
#     threshold = filters.threshold_minimum(totals)
#     atoms = (np.array(totals)> threshold)
#     return totals, threshold, atoms  
