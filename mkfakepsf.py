#! /usr/bin/env python2.5
"""
2010.09.18
S.Rodney

Module for constructing fake SN point spread functions.

This module is called by plantfake to build a set of fake SN psf
images.  For each fake SN we define a separate psf profile for each of
the flt frames that it will be planted in.  The psf profile is
constructed by TinyTim to incorporate geometric distortion effects.
Each psf is scaled to a specified flux value, mimicking real SN
colors.

TODO: add an option for using measured psfs in place of TinyTim.
"""
import exceptions
import os
import sys

# WFC3-IR inter-pixel capacitance kernel
# Outdated kernel reported in TinyTim WFC3 manual: 
# # http://tinytim.stsci.edu/static/TinyTim_WFC3.pdf
# KERNEL_WFC3IR = [ [ 0.002,  0.038,  0.002], 
#           [ 0.038,  0.840,  0.038], 
#           [ 0.002,  0.038,  0.002] ]

# Newer kernel from ISR WFC3-2008-41
KERNEL_WFC3IR = [ [ 0.0007, 0.025, 0.0007 ],
                  [ 0.025,  0.897, 0.025 ], 
                  [ 0.0007, 0.025, 0.0007 ] ]

# diameter of psf (arcsec) 
# (the psf size recommended by tinytim for WFC3-IR, UVIS and 
#  for ACS is about 3.0 arcsec.  We are measuring photometry
# with aperture sizes less than 0.5"
PSFSIZE = 2.0


def mkTinyTimPSF( x, y, fltfile, ext=1,
                  fileroot='tinytim', psfdir='tinytim',
                  specfile='flat_flam_tinytim.dat',
                  verbose=False, clobber=False ):
    """ run tinytim to construct a model psf 
    in the distorted (flt) frame
    """
    # TODO : use a redshifted SN spectrum !!!
    #   (will have to build each psf separately)

    # TinyTim generates a psf image centered on 
    # the middle of a pixel
    #  we use "tiny3 SUB=5" for 5x sub-sampling,
    #  then interpolate and shift the sub-sampled psf
    #  to recenter it away from the center of the pixel
    #  Re-bin into normal pixel sampling, and then 
    #  convolve with the charge diffusion kernel from
    #  the fits header. 

    import time
    from numpy import iterable, zeros
    from scipy.ndimage import zoom
    # from util.fitsio import tofits
    import cntrd
    import pyfits

    imhdr0 = pyfits.getheader( fltfile, ext=0 )
    imhdr = pyfits.getheader( fltfile, ext=ext )
    instrument = imhdr0['instrume']
    detector = imhdr0['detector']
    if ( instrument=='WFC3' and detector=='IR' ) : 
        camera='ir'
        filt = imhdr0['filter']
        ccdchip = None
    elif ( instrument=='WFC3' and detector=='UVIS' ) : 
        camera='uvis'
        filt = imhdr0['filter']
        ccdchip = imhdr['CCDCHIP']
    elif ( instrument=='ACS' and detector=='WFC' ) : 
        camera='acs'
        filter1 = imhdr0['filter1']
        filter2 = imhdr0['filter2']
        if filter1.startswith('F') : filt=filter1
        else : filt=filter2
        ccdchip = imhdr['CCDCHIP']


    if not os.path.isfile( specfile ) : 
        if 'TINYTIM' in os.environ :
            tinytimdir = os.environ['TINYTIM']
            specfile = os.path.join( tinytimdir, specfile )
    if not os.path.isfile( specfile ) :
        thisfile = sys.argv[0]
        if thisfile.startswith('ipython'): thisfile = __file__
        thisdir = os.path.dirname( thisfile )
        specfile = os.path.join( thisdir, specfile )
    if not os.path.isfile( specfile ) :
        raise exceptions.RuntimeError("Can't find TinyTim spec file %s"%
                                      os.path.basename(specfile) )
    if verbose :
        print( "Using TinyTim spectrum file : %s"%specfile )

    if not iterable(x) : x = [x]
    if not iterable(y) : y = [y]

    if iterable( ext ) :
        extname = ''.join( str(extbit).lower() for extbit in ext )
    else :
        extname = str(ext).lower()
    coordfile = "%s.%s.coord"%(os.path.join(psfdir,fileroot),extname)

    if not os.path.isdir(psfdir) : os.mkdir( psfdir )
    newcoords=True
    if os.path.isfile( coordfile ) and not clobber :
        print( "%s exists. Not clobbering."%coordfile )
        newcoords=False

    psfstamplist = []
    if newcoords: fout = open( coordfile ,'w')
    allexist = True

    i=0
    for xx,yy in zip(x,y):
        # here we give the integer component of the x,y coordinates.
        # after running tiny3, we will shift the psf to account for the 
        # fractional coordinate components
        if newcoords: print >>fout,"%6i %6i"%(int(xx),int(yy))
        if len(x)<=100:
            psfstamp = os.path.join(psfdir,"%s.%s.%02i.fits"%(fileroot,extname,i))
        else : 
            psfstamp = os.path.join(psfdir,"%s.%s.%03i.fits"%(fileroot,extname,i))
        if not os.path.isfile( psfstamp ) : allexist = False
        psfstamplist.append( psfstamp)
        i+=1
    if newcoords: fout.close()

    if allexist and not clobber : 
        print( "All necessary tinytim psf files exist. Not clobbering.")
        return( psfstamplist )

    queryfile = "%s.%s.query"%(os.path.join(psfdir,fileroot),extname)
    fout = open( queryfile ,'w')
    despace = 0.0 # as of tinytim v7.1 : must provide 2ndary mirror despace
    if camera == 'ir' : 
        print >> fout, """23\n @%s\n %s\n 5\n %s\n %.1f\n %.1f\n %s.%s."""%(
            coordfile, filt.lower(), specfile, 
            PSFSIZE, despace, os.path.join(psfdir,fileroot), extname )
    elif camera == 'uvis' : 
        print >> fout,"""22\n %i\n @%s\n %s\n 5\n %s\n %.1f\n %.1f\n %s.%s."""%(
            ccdchip, coordfile, filt.lower(), specfile, 
            PSFSIZE, despace, os.path.join(psfdir,fileroot), extname )
    elif camera == 'acs' : 
        print >> fout,"""15\n %i\n @%s\n %s\n 5\n %s\n %.1f\n %.1f\n %s.%s."""%(
            ccdchip, coordfile, filt.lower(), specfile, 
            PSFSIZE, despace, os.path.join(psfdir,fileroot), extname )
    fout.close()
    
    # run tiny1 to generate the tinytim paramater file
    command1 =  "cat %s | %s %s.%s.in"%(
        queryfile, 'tiny1',
        os.path.join(psfdir,fileroot), extname )
    if verbose : print command1
    os.system( command1 )

    # run tiny2 to generate the distortion-free psfs
    command2 =  "%s %s.%s.in"%(
        'tiny2', os.path.join(psfdir,fileroot), extname )
    if verbose : print command2
    os.system( command2 ) 

    xgeo_offsets = []
    ygeo_offsets = []
    # fluxcorrs = []  # 2014.07.18 Flux correction disabled by Steve
    #run tiny3 and measure the how much offset the geometric distortion adds
    for Npsf in range(len(x)): 
        command3 =  "%s %s.%s.in POS=%i"%(
            'tiny3', os.path.join(psfdir,fileroot), extname, Npsf )
        if verbose : 
            print time.asctime() 
            print command3
        os.system( command3 )
        #Calculate the expected center of the image
        #Get the dimensions of the stamp. 
        if len(x) <= 100 :
            this_stamp = '%s.%s.%02i.fits' %(os.path.join(psfdir,fileroot),extname,Npsf)
        else : 
            this_stamp = '%s.%s.%03i.fits' %(os.path.join(psfdir,fileroot),extname,Npsf)
        xdim = int(pyfits.getval(this_stamp,'NAXIS1'))
        ydim = int(pyfits.getval(this_stamp,'NAXIS2'))

        #The center will be in dimension/2 + 1
        xcen = float(xdim/2 + 1)
        ycen = float(ydim/2 + 1)
        
        #run phot to measure the true center position
        if instrument =='WFC3': instrument = instrument +'_'+detector
        pixscale = getpixscale( this_stamp )
        fwhmpix = 0.13 / pixscale  # approximate HST psf size, in pixels
        meas_xcen, meas_ycen = cntrd.cntrd(this_stamp,xcen,ycen,fwhmpix)
        #Subtract the expected center from the measured center
        #Save the offsets
        xgeo_offsets.append(meas_xcen - xcen)
        ygeo_offsets.append(meas_ycen - ycen)
        # fluxcorrs.append(meas_fluxcorr)
        #Move this stamp so that it doesn't get overwritten.
        os.rename(this_stamp,os.path.splitext(this_stamp)[0]+'_tiny3.fits')
        
        
    # run tiny3 to add in geometric distortion and 5x sub-sampling
    for Npsf in range(len(x)): 
        command3 =  "%s %s.%s.in POS=%i SUB=5"%(
            'tiny3', os.path.join(psfdir,fileroot), extname, Npsf )
        if verbose : 
            print time.asctime() 
            print command3
        os.system( command3 ) 

    outstamplist = []
    for xx,yy,psfstamp,xgeo,ygeo in zip( x,y,psfstamplist,xgeo_offsets,ygeo_offsets):
        if verbose : print("sub-sampling psf at %.2f %.2f to 0.01 pix"%(xx,yy))

        # read in tiny3 output psf (sub-sampled to a 5th of a pixel)
        psfim = pyfits.open( psfstamp )
        psfdat = psfim[0].data.copy()
        hdr = psfim[0].header.copy()
        psfim.close()
        #If the number of pixels is even then the psf is centered at pixel n/2 + 1
        # if you are one indexed or n/2 if you are zero indexed. 
        #If the psf image is even, we need to pad the right side with a row (or column) of zeros
        if psfdat.shape[0] % 2 == 0: 
            tmpdat = zeros([psfdat.shape[0]+1,psfdat.shape[1]])
            tmpdat[:-1,:] = psfdat[:,:]
            psfdat = tmpdat
        if psfdat.shape[1] % 2 == 0: 
            tmpdat = zeros([psfdat.shape[0],psfdat.shape[1]+1])
            tmpdat[:,:-1] = psfdat[:,:]
            psfdat = tmpdat
       
        #Now the center of the psf is exactly in the center of the image and the psf image has
        #odd dimensions
        
        #TinyTim returns the psf subsampled at a 5th of a pixel 
        #but not necessarily psfim.shape % 5 == 0. 
        
        #Now we need to pad the array with zeros so that center of the image will be in 
        #the center of both the 1/5th pixel image and the psf at native scale. 
        #As the center of the psf is at the center, all we need to do is add pixels to both
        #sides evenly until we have an integer number of native pixels. 
        # All of the rules assume that the dimensions are odd which makes the rules
        #a little confusing.
         
        xpad,ypad = psfdat.shape[1] % 5, psfdat.shape[0] % 5
        
        if xpad == 2:
            tmpdat = zeros([psfdat.shape[0],psfdat.shape[1]+8])
            tmpdat[:,4:-4] = psfdat[:,:]
            psfdat = tmpdat
        elif xpad == 4:
            tmpdat = zeros([psfdat.shape[0],psfdat.shape[1]+6])
            tmpdat[:,3:-3] = psfdat[:,:]
            psfdat = tmpdat
        elif xpad == 1:
            tmpdat = zeros([psfdat.shape[0],psfdat.shape[1]+4])
            tmpdat[:,2:-2] = psfdat[:,:]
            psfdat = tmpdat
        elif xpad == 3:
            tmpdat = zeros([psfdat.shape[0],psfdat.shape[1]+2])
            tmpdat[:,1:-1] = psfdat[:,:]
            psfdat = tmpdat
            
        if ypad == 2:
            tmpdat = zeros([psfdat.shape[0]+8,psfdat.shape[1]])
            tmpdat[4:-4,:] = psfdat[:,:]
            psfdat = tmpdat
        elif ypad == 4:
            tmpdat = zeros([psfdat.shape[0]+6,psfdat.shape[1]])
            tmpdat[3:-3,:] = psfdat[:,:]
            psfdat = tmpdat
        elif ypad == 1:
            tmpdat = zeros([psfdat.shape[0]+4,psfdat.shape[1]])
            tmpdat[2:-2,:] = psfdat[:,:]
            psfdat = tmpdat
        elif ypad == 3:
            tmpdat = zeros([psfdat.shape[0]+2,psfdat.shape[1]])
            tmpdat[1:-1,:] = psfdat[:,:]
            psfdat = tmpdat

        #Add 2 extra pixels on both sides (+ 400 in each dimension) to account for the fractional shift
        #and the geometric distortion

        psfdat100 = zeros([psfdat.shape[0]*20 + 400, psfdat.shape[1]*20 + 400])
        
        #Calculate the fractional shifts
        xfrac,yfrac = int(round(xx % 1 * 100)), int(round(yy % 1 * 100))
        if yfrac == 100: yfrac = 0
        if xfrac == 100: xfrac = 0
        
        #Add the geometric distorition offsets of the centroid into xfrac and yfrac. 
        #This makes the assumption that the distortion centroid offsets are less than 1 pixel
        #This has been the case for all of my tests.
        if verbose: print('Adding %0.2f, %0.2f to correct the center of the psf for geometric distortion' % (xgeo,ygeo))
        xfrac -= int(100*xgeo)
        yfrac -= int(100*ygeo)
        
        if verbose : print("    Interpolating and re-sampling with sub-pixel shift")

        # interpolate at a 20x smaller grid to get 
        #  sub-sampling at the 100th of a pixel level
        #Right now we use the ndimage zoom function which does a spline interpolation, but is fast
        try : psfdat100[200+yfrac:-(200-yfrac),200+xfrac:-(200-xfrac)] = zoom(psfdat,20)
        except ValueError as e:
            print( e )
            import pdb; pdb.set_trace()
        
        # re-bin on a new grid to get the psf at the full-pixel scale          
        psfdat1 = rebin(psfdat100, 100)
        #remove any reference to psfdat100 in hopes of it getting garbage collected 
        #as it is by far the biggest thing we have in memory
        del psfdat100
        psfdat1 = psfdat1 / psfdat1.sum()
    
        # Blur the re-binned psf to account for detector effects:
        # For  UVIS and ACS, read in the charge diffusion kernel
        # For IR, we use a fixed IR inter-pixel capacitance kernel, defined above
        if verbose : print("    convolving with charged diffusion or inter-pixel capacitance kernel")
        if camera == 'ir': kernel=KERNEL_WFC3IR
        else : kernel = getCDkernel( psfim[0].header ) 

        psfdat2 = convolvepsf( psfdat1, kernel )

        # 2014.07.18 : Disabled by Steve
        # Rescale the TinyTim psf to match the measured aperture corrections.
        # psfdat2 *= fluxcorr
        # if verbose : print('Applying a %f flux correction to the TinyTim psf.' % fluxcorr)

        # write out the new recentered psf stamp
        outstamp = psfstamp.replace('.fits','_final.fits')
        hdr['naxis1']=psfdat2.shape[1]
        hdr['naxis2']=psfdat2.shape[0]
        psfim.writeto( outstamp, psfdat2, header=hdr, clobber=True )
        outstamplist.append( outstamp )
        if verbose : print("    Shifted, resampled psf written to %s"%outstamp)

    # return a list of  psf stamps 
    return( outstamplist )    
    

def rebin(a, factor):
    '''
    '''
    #Some fancy indexing to make this sum easier, hopefully it is still fast
    return a.reshape(a.shape[0]/factor,factor,a.shape[1]/factor,factor).sum(axis=3).sum(axis=1)
    
def convolvepsf( psfdat, kernel ):
    """ convolve with the inter-pixel capacitance kernel
    (see ISR WFC3-2008-41:
    http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2008-41.pdf

    TODO : include ACS and UVIS charge diffusion kernels
    """
    from scipy import ndimage
    blurpsf = ndimage.convolve( psfdat, kernel, mode='constant', cval=0.0 )

    # rescale output pixels to fix integral at unity
    blurpsf = blurpsf / blurpsf.sum()
    return( blurpsf )

def getCDkernel( imheader ) :
    """ read the Charge Diffusion kernel planted by tinytim 
    out of the comments section of the given fits header.
    NOTE: we're assuming that the kernel lives in the 
    last three header cards.  Is that safe?
    """
    kernel = []
    for i in range(-3,0) : 
        kernel.append( [ float(val) for val in imheader[i].split() ] )
    return( kernel )

def get_fake_centroid(filename,x,y,instrument,filt):
    """
    Locate the center of a fake psf

    INPUTS: The fake-SN psf image in filename, the expected x,y position
    of the center of the psf, the instrument and filter being modeled.
    RETURNS: xcentroid, ycentroid, fluxcorr
    """
    from pyraf import iraf

    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.unlearn(iraf.apphot.phot)
    iraf.unlearn(iraf.datapars)
    iraf.unlearn(iraf.centerpars)
    #Use the centroid algorithm right now as it seems more robust to geometric distortion.
    iraf.centerpars.calgorithm = 'centroid'
    iraf.centerpars.cbox = 5.0

    iraf.unlearn(iraf.fitskypars)
    iraf.unlearn(iraf.photpars)
    photparams = {
        'interac':False,
        'radplot':False,
        }
    iraf.datapars.readnoise = 0.0
    iraf.datapars.itime = 1.0
    iraf.datapars.epadu = 1.0

    # iraf.digiphot.apphot.fitskypars :
    iraf.unlearn(iraf.fitskypars)
    iraf.fitskypars.salgorithm = 'constant'
    iraf.fitskypars.skyvalue = 0.0

    # iraf.digiphot.apphot.photpars :
    iraf.unlearn(iraf.photpars)
    iraf.photpars.weighting = 'constant'
    iraf.photpars.apertures = 20   # TODO : set this more intelligently !
    iraf.photpars.zmag = 25
    iraf.photpars.mkapert = False

    #Write the coordinate file starting as position x and y
    coxyfile = 'centroid.xycoo'
    coxy = open(coxyfile, 'w')
    print >> coxy, "%10.2f  %10.2f" % (x,y)
    coxy.close()
    if os.path.exists('centroid.mag'): os.remove('centroid.mag')
    iraf.phot(image=filename, skyfile='', coords=coxyfile, output='centroid.mag',
              verify=False, verbose=True, Stdout=1, **photparams)
    f = open('centroid.mag', 'r')
    maglines = f.readlines()
    f.close()
    xcentroid = float(maglines[76].split()[0])
    ycentroid = float(maglines[76].split()[1])

    return xcentroid,ycentroid

def getpixscale( fitsfile, returntuple=False, ext=0 ):
    """ Compute the pixel scale of the reference pixel in arcsec/pix in
    each direction from the fits header cd matrix.
    With returntuple=True, return the two pixel scale values along the x and y
    axes.  For returntuple=False, return the average of the two.
    The input fitsfile may be a string giving a fits filename, a
    pyfits hdulist or hdu.
    """
    from math import sqrt
    import pyfits

    hdr = pyfits.getheader( fitsfile, ext=ext )
    if 'CD1_1' in hdr :
        cd11 = hdr['CD1_1']
        cd12 = hdr['CD1_2']
        cd21 = hdr['CD2_1']
        cd22 = hdr['CD2_2']

        # define the sign based on determinant
        det = cd11*cd22 - cd12*cd21
        if det<0 : sgn = -1
        else : sgn = 1

        if cd12==0 and cd21==0 :
            # no rotation: x=RA, y=Dec
            cdelt1 = cd11
            cdelt2 = cd22
        else :
            cdelt1 = sgn*sqrt(cd11**2 + cd12**2)
            cdelt2 = sqrt(cd22**2 + cd21**2)
    elif 'CDELT1' in hdr.keys() and (hdr['CDELT1']!=1 and hdr['CDELT2']!=1) :
        cdelt1 = hdr['CDELT1']
        cdelt2 = hdr['CDELT2']

    cdelt1 = cdelt1  * 3600.
    cdelt2 = cdelt2  * 3600.

    if returntuple :
        return( cdelt1, cdelt2 )
    else :
        return( (abs(cdelt1)+abs(cdelt2)) / 2. )

