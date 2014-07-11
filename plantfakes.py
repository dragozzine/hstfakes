#! /usr/bin/env python
"""
2010.09.19
S.Rodney

Top-level module for planting fake SN psf profiles into flt images.

This module may be called from the command line to plant fake psf images
into flt images.

The core function is addfltfakes(), which takes as input the name of a single
flt file, and a catalog file specific to that flt image with this format :

 # header comments
 #  X    Y    FLX   SCI
  10.1 20.2   0.3    1
  40.5 50.5   0.6    2

The columns x, y and flux are required.  The 'sci' column is required for
files with multiple SCI extensions for multiple CCDs. It should typically have
values of either 1 or 2 to indicate fits extensions 'SCI,1' and 'SCI,2' (not
necessarily corresponding to 'CCD1' and 'CCD2').

If a header line is absent, we assume the first 3 columns are x,y,flux.  If
no SCI extension column is provided, fakes are planted exclusively in the first
SCI extension.

If additional columns are present, they are recorded in the fits header
for each planted fake -- but only the first 3 characters of the column name
are used for the fits header keyword.  For example, a column labeled 'MAG_F814W'
would be encoded in the header as 'FK000MAG', where the 000 indicates the index
of each planted fake.

It is always assumed that the flux  values correspond to the appropriate
brightness units of the given flt file (i.e. ELECTRONS for ACS and WFC3 UVIS,
ELECTRONS/SEC for WFC3-IR).

"""

#We have removed all of the zeropoint dictionaries and now point a common location for these
#dicts in hstsnphot.
#from tools import hstsnphot

def addfltfakes( fltfile, fltcatfile, psf='TinyTim', verbose=False, clobber=False ):
    """ Add fakes into the given flt.fits file <fltfile> at positions and fluxes
    as specified in the given catalog <fltcatfile>.

    OPTIONS :
      psf    -   The PSF model to use. A string giving either the path to a
                 fits file with the complete PSF model, or the word 'TinyTim'
                 to indicate that TinyTim is to be used for generating
                 position-dependent PSFs.
    """
    import os
    import pyfits
    import numpy as np
    #import util
    import mkfakepsf
    # from numpy import unique, where
    from astropy.io import ascii
    # TODO : check for fake info in header before clobbering (clobber by default?)

    fakefile = fltfile.replace('_flt.fits','_fake_flt.fits')

    # read in the catalog of fake SN locations and magnitudes
    im = pyfits.open( fltfile )

    fltcat = ascii.read( fltcatfile, format='commented_header',  header_start=-1)
    x = fltcat.columns['X'] if 'X' in fltcat.colnames else fltcat.columns[0]
    y = fltcat.columns['Y'] if 'Y' in fltcat.colnames else fltcat.columns[1]
    flux = fltcat.columns['FLX'] if 'FLX' in fltcat.colnames else fltcat.columns[2]
    Nfakes = len(flux)

    sci = np.ones( Nfakes )
    if len(fltcat.colnames)>3 :
        if 'SCI' in fltcat.colnames :
            sci = fltcat.columns['SCI']
        else :
            sci = fltcat.columns[3]

    # In ACS and UVIS we have multiple SCI extensions in the flt
    # files corresponding to the two chips.  Our fake SNe
    # get added to one chip at a time.
    hdr0 = im[0].header
    exptime = hdr0['EXPTIME']
    hdr0.add_blank()
    hdr0.update( 'NFAKES', len(flux) , comment='Number of fake PSFs added')
    scilist = np.unique( sci )
    for isci in scilist :
        ithissci = np.where(sci==isci)[0]

        if psf.lower() == 'tinytim' :
            # make all the TinyTim psf stamps for this chip
            tinytimroot = fltfile.replace( '_flt.fits','_tinytim')
            psfstamplist = mkfakepsf.mkTinyTimPSF(
                x[ithissci], y[ithissci], fltfile, ext=['SCI',isci],
                fileroot=tinytimroot, verbose=verbose, clobber=clobber>1 )
        else :
            # read in a single psf stamp image
            psfimarray = pyfits.getdata( psf )

        # set up for fake PSF data in header
        hdr = im['SCI',isci].header
        hdr.add_blank()
        hdr.update( 'NFAKES', len(flux) , comment='Number of fake PSFs added')

        # set the 'gain' for converting pixel values to electrons.
        if hdr['BUNIT'] == 'ELECTRONS' : gain=1.0
        elif hdr['BUNIT'] == 'ELECTRONS/S' : gain=exptime

        # add fakes to image data array and header
        imdat = im['SCI',isci].data
        for ifake in ithissci :

            # Update the SCI extension header and the primary header
            for col in fltcat.colnames :
                val = fltcat.columns[col][ifake]
                key = 'FK%03i%-3s'%(ifake,col[:3].upper())
                val,fmt = setformat( val )
                hdr.update( key, fmt%val , comment='Fake %03i %s'%(ifake,col) )
                hdr0.update( key, fmt%val , comment='Fake %03i %s'%(ifake,col) )

            # read in the psf model
            if psf.lower() == 'tinytim' :
                psfimarray = pyfits.getdata( psfstamplist[ifake] )

            # scale the psf to the desired flux level, add it to the flt image
            addpsfimage( imdat, x[ifake], y[ifake], psfimarray,
                         scale=flux[ifake], gain=gain, usecopy=False )
            #import pdb; pdb.set_trace()
            # print( imdat.sum()-im['SCI',isci].data.sum() )
        #im['SCI',isci].data = imdat
        #im['SCI',isci].header = hdr

    # write out the faked flt file
    im.writeto( fakefile, clobber=True )
    return( fakefile )


def setformat( val ) :
    """ Define a format string for the given value that is suitable for
     printing to a fits header.
    """
    if isinstance( val, float ) :
        if val>0.01 :
            val = round(val,2)
            fmt = '%.2f'
        elif val>0.001 :
            val = round(val,3)
            fmt = '%.3f'
        else :
            fmt = '%.4e'
    elif isinstance( val, int ) :
        fmt = '%i'
    else :
        val = str( val )
        fmt = '%s'
    return( val, fmt )


def addpsfimage( imdat, x, y, fakestamp, scale=1.0, gain = 1.0,
                 usecopy=False, verbose=False ):
    """ 
    insert a single fake SN (given in the numpy ndarray fakestamp)
    into the image data array imdat at position (x,y), scaled to 
    the given total flux.

    scale - multiply the normalized PSF image by this value to set the total
             flux of the planted star.

    gain -  the (inverse) gain, i.e.  electrons / ADU
            This gives the conversion factor from the image's brightness units
            into units of electrons, for the purpose of adding shot noise.
            If the image data is in electrons, then use gain=1.0
            If its in electrons/sec then use gain = exposure_time
    """
    from numpy import copy, random

    # NOTE : using int instead of round here is important
    #  for getting the fakes to align properly when the 
    # tinytim psf stamps have different sizes for diff't bands
    # Now the psfstamp always has odd dimensions so we can just do int math
    #Since we are zero indexed we don't need to add or subtract 1 to the center position.
    xc = fakestamp.shape[1]/2
    yc = fakestamp.shape[0]/2
    x1 = int(x) - xc - 1 
    x2 = int(x) + (fakestamp.shape[1]-xc) - 1 
    y1 = int(y) - yc - 1 
    y2 = int(y) + (fakestamp.shape[0]-yc) - 1

    if imdat[ y1:y2, x1:x2].shape != fakestamp.shape :
        if verbose : print("not enough room for stamp at %.2f %.2f"%(x,y))
        if usecopy : return( copy(imdat) )
        else: return None

    # # TODO : clean up trimming so that sizes always match
    # # check that stamp fits inside image boundaries
    # if( y1>imdat.shape[0]-1 or x1>imdat.shape[1]-1 or 
    #     y2<1 or x2<1 ) : 
    #     if usecopy : return(copy(imdat))
    #     else : return(0)
    # 
    # # if necessary, trim to fit 
    # y1stamp,y2stamp,x1stamp,x2stamp = 0,fakestamp.shape[0],0,fakestamp.shape[1]
    # if y1<2 : 
    #     y1stamp = 2-int(y1)
    #     y1=2
    # if x1<2 : 
    #     x1stamp = 2-int(x1)
    #     x1=2
    # if y2>imdat.shape[0]-2 : 
    #     y2stamp = fakestamp.shape[0] - int((y2-imdat.shape[0]-2))
    #     y2=imdat.shape[0]-2
    # if x2>imdat.shape[1]-2 : 
    #     x2stamp = fakestamp.shape[1] - int((x2-imdat.shape[1]-2))
    #     x2=imdat.shape[1]-2
    # fakestamp = fakestamp[ y1stamp:y2stamp, x1stamp:x2stamp ]
    
    #Right now we put a floor on the minimum for the psf so the noise calculation makes sense
    fakestamp[fakestamp < 1e-7] = 1e-7

    # scale fake star to the desired net flux
    fakescaled = fakestamp / fakestamp.sum() * scale

    # add in shot noise to each pixel
    fakescaled = random.poisson(fakescaled * float(gain))/float(gain)
        
    # add the stamp to the image data array
    if usecopy :
        # add the fake stamp to a copy of the input
        # image, so as to not tamper with the input image
        imcopy = copy( imdat )
        imcopy[ y1:y2, x1:x2] += fakescaled
        # import pdb; pdb.set_trace()
        return(imcopy)
    else : 
        # adjust the input image itself, saving some memory
        imdat[ y1:y2, x1:x2] += fakescaled
        return None


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser( description= 'Plant fake PSFs into an FLT image.')

    # Required positional argument
    parser.add_argument('fltfile', help='FLT file in which to plant the fakes.')
    parser.add_argument('catalog', help='ASCII Catalog with X,Y,FLUX[,SCI] columns.')

    # optional arguments :
    parser.add_argument('--morehelp', action='store_true', help='Print more detailed help message and exit.')
    parser.add_argument('--psf', metavar='TinyTim', default='TinyTim', help='PSF model to use. (TinyTim or a fits file name)')
    parser.add_argument('--verbose', action='store_true', help='Turn on verbose mode. [False]', default=False)
    parser.add_argument('--clobber', action='store_true', help='Turn on clobber mode. [False]', default=False)

    argv = parser.parse_args()


    if argv.morehelp :
        print( __doc__ )
        parser.print_usage()
        sys.exit(0)

    addfltfakes( argv.fltfile, argv.catalog, psf=argv.psf,
                 verbose=argv.verbose, clobber=argv.clobber )
    return 0

if __name__ == "__main__":
    main()
