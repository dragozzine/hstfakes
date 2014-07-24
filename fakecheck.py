#! /usr/bin/env python
# 2011.09.08 S.Rodney
# 
"""
Use aperture photometry to measure the magnitudes of fake SNe
at known locations in a faked diff image.

Extracts the fake SN info directly from the fits header
"""
from numpy import array, ndarray

# dictionary of AB and Vega zero points for all filters
ZPT_ACS_VEGA = {'F435W': 25.76695,'F475W': 26.16252,'F555W': 25.72747,'F606W': 26.40598,'F625W': 25.74339,'F775W': 25.27728,'F814W': 25.51994,'F850LP':24.32305,}
ZPT_ACS_AB = {'F435W': 25.65777,'F475W': 26.05923,'F555W': 25.7184, 'F606W': 26.49113,'F625W': 25.90667,'F775W': 25.66504,'F814W': 25.94333,'F850LP':24.84245,}
ZPT_WFC3_IR_AB = {'F105W':26.27,'F110W':26.83,'F125W':26.25,'F140W':26.46,'F160W':25.96,'G141':25.96,}
ZPT_WFC3_IR_VEGA = {'F105W':25.63,'F110W':26.07,'F125W':25.35,'F140W':25.39,'F160W':24.70,'G141':24.70 }
ZPT_WFC3_UVIS_AB = {'F218W':22.94,'F225W':24.06,'F275W':24.14,'F336W':24.64,'F350LP':26.94,'F390W':25.37,'F438W':24.83,'F475W':25.69,'F555W':25.78,'F600LP':25.85,'F606W':26.08,'F625W':25.52,'F775W':24.85,'F814W':25.09,}
ZPT_WFC3_UVIS_VEGA = {'F218W':21.25,'F225W':22.40,'F275W':22.65,'F336W':23.46,'F350LP':26.78,'F390W':25.15,'F438W':24.98,'F475W':25.79,'F555W':25.81,'F600LP':25.53,'F606W':25.99,'F625W':25.37,'F775W':24.47,'F814W':24.67,}

ZPT = dict(ZPT_WFC3_IR_VEGA)
ZPT.update( ZPT_ACS_VEGA )
ZPT.update( ZPT_WFC3_UVIS_VEGA )


# aperture corrections : 
RAP_PIX  = array( [2,3,4,5] )
RAP_ASEC  = { 
    'UVIS':RAP_PIX*0.04,
    'ACS':RAP_PIX*0.05,
    'IR':RAP_PIX*0.09 }
APCOR = { 
    'F125W':array( [0.566, 0.270, 0.182, 0.150 ]),
    'F160W':array( [0.685, 0.361, 0.223, 0.168 ]),
    'F606W' :array([0.461,0.247,0.183,0.151,0.088 ]),
    'F625W' :array([0.461,0.247,0.183,0.151,0.088 ]), # JUST A COPY OF 606!!
    'F775W' :array([0.533,0.274,0.197,0.162,0.087 ]),
    'F814W' :array([0.543,0.300,0.203,0.173,0.087 ]),
    'F850LP':array([0.584,0.330,0.217,0.184,0.117 ]),
    }

def getApCorrection( filter, radius, radunit='pix' ) :
    """ returns the aperture correction for the given 
    filter and radius """
    if radunit == 'pix' : 
        rap = RAP_PIX
    else : 
        if filter in ['F125W','F160W'] : rap = RAP_ASEC['IR']
        elif filter in ['F350LP'] : rap = RAP_ASEC['UVIS']
        else : rap = RAP_ASEC['ACS']

    apcorlist = APCOR[filter] 
    irap = rap.tolist().index( radius )
    return( apcorlist[ irap ] )


def plotFakeMags( imfile, coofile=None, magfile=None, 
                  interactive=False) : 
    """ get coords, do photometry, plot results
    if savetype is specified, then the plot is saved with the 
    given extension ('png','pdf','eps', etc.)
    """
    if interactive : 
        from pylab import clf, gca, axes, savefig, ioff
        clf()
        ax = gca()
    else : 
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        #from matplotlib.backends.backend_agg import FigureCanvasAgg
        ax = plt.axes()
    
    inmaglist = getInputMags( imfile )
    coofile = mkFakeCoordFile( imfile, coofile=coofile )
    # outmaglist = runApphot(imfile, coofile=coofile, magfile=magfile)
    photometry = runApphot(imfile, coofile=coofile, magfile=magfile)
    
    if photometry in [ None, [] ] :  return( None )

    ax.plot( [inmaglist.min()-0.2, inmaglist.max()+0.2 ], 
             [inmaglist.min()-0.2, inmaglist.max()+0.2 ], 'k-' )

    colors = ['r','darkorange','g','b','darkorchid','k']
    for i,c in zip( range( len( RAP_PIX ) ), colors ) : 
        ax.errorbar( inmaglist+RAP_PIX[i]*0.01, 
                     photometry.mags(RAP_PIX[i])+RAP_PIX[i]*0.01, 
                     photometry.magerrs(RAP_PIX[i]), 
                     marker='o', color=c, ls=' ',
                     label='%i pix'%RAP_PIX[i])
    
    ax.legend( loc='lower right', numpoints=1, title='Aperture' )
    ax.set_xlabel('input magnitude + 0.01$\cdot$R$_{ap}$')
    ax.set_ylabel('recovered magnitude + 0.01$\cdot$R$_{ap}$')
    ax.set_title('%s Fake Mag Check'%imfile)

    ax.set_xlim( [inmaglist.min()-0.2, inmaglist.max()+0.2 ] )
    ax.set_ylim( [inmaglist.min()-0.2, inmaglist.max()+0.2 ] )

    plt.savefig( imfile.replace( '.fits','_magcheck.pdf' ) )


def getInputMags( imfile ) :
    """ get a list of input fake SN mags """
    import pyfits

    maglist = []
    hdr = pyfits.getheader( imfile ) 
    Nfake  = hdr['NFAKESNE']
    for i in range( Nfake ) : 
        x = hdr['FAKE%03iX'%i]
        y = hdr['FAKE%03iY'%i]
        mag = hdr['FAKE%03iM'%i]
        maglist.append( mag ) 

    return( array(maglist) )

    
def mkFakeCoordFile( imfile, coofile=None ):
    """ extract x,y coordinates 
    of fake SNe from the imfile header """
    import pyfits

    if not coofile : 
        coofile = imfile.replace('.fits', '.coo')

    fout = open(coofile,'w')
    hdr = pyfits.getheader( imfile ) 
    Nfake  = hdr['NFAKESNE']
    for i in range( Nfake ) : 
        x = hdr['FAKE%03iX'%i]
        y = hdr['FAKE%03iY'%i]
        print >> fout, '%12.3f  %12.3f'%(x,y)
    fout.close()

    return( coofile ) 



def runApphot(imfile, coofile=None, magfile=None):
    """ 
    use iraf.digiphot.apphot 
    to collect aperture photometry
    """
    import pyraf
    from pyraf import iraf
    import os
    import pyfits

    if not coofile : 
        coofile = imfile.replace('.fits', '.coo')
    if not magfile : 
        magfile = imfile.replace('.fits', '.mag')

    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    
    if os.path.exists( magfile ) : 
        os.remove( magfile ) 

    hdr = pyfits.getheader( imfile )
    if 'FILTER' in hdr: 
        filt=hdr['FILTER']
        filtkey='FILTER'
    elif 'FILTER1' in hdr :
        filt=hdr['FILTER1']
        filtkey='FILTER1'
    if filt.startswith('CLEAR'):
        filt = hdr['FILTER2']
        filtkey='FILTER2'

    if filt not in APCOR.keys() :
        return( [] )
    
    # iraf.digiphot.apphot.datapars : 
    iraf.datapars.scale = 1.0
    iraf.datapars.fwhmpsf = 2.5
    iraf.datapars.emission = True
    iraf.datapars.sigma = 'INDEF'
    iraf.datapars.datamin = 'INDEF'
    iraf.datapars.datamax = 'INDEF'
    iraf.datapars.noise =  'constant'
    iraf.datapars.ccdread = ''
    iraf.datapars.gain = ''
    iraf.datapars.readnoise = 0.0
    iraf.datapars.epadu = 1.0
    iraf.datapars.exposure = ''
    iraf.datapars.airmass = ''
    iraf.datapars.filter = filt
    iraf.datapars.obstime = ''
    iraf.datapars.itime = 1.0
    iraf.datapars.xairmass = 'INDEF'
    iraf.datapars.ifilter = 'INDEF'
    iraf.datapars.otime = 'INDEF'
    
    # iraf.digiphot.apphot.centerpars : 
    iraf.unlearn( iraf.centerpars )
    iraf.centerpars.calgorithm = 'centroid'
    iraf.centerpars.cbox = 3.0
    iraf.centerpars.cthreshold = 0.0
    iraf.centerpars.minsnratio = 1.0
    iraf.centerpars.cmaxiter = 10.0
    iraf.centerpars.maxshift = 1.0
    iraf.centerpars.clean = False
    iraf.centerpars.rclean = 1.0
    iraf.centerpars.rclip = 2.0
    iraf.centerpars.kclean = 3.0
    iraf.centerpars.mkcenter = False

    # iraf.digiphot.apphot.fitskypars : 
    iraf.unlearn( iraf.fitskypars )
    iraf.fitskypars.salgorithm = 'median'
    iraf.fitskypars.annulus = 25.0
    iraf.fitskypars.dannulus = 40.0
    iraf.fitskypars.skyvalue = 0.0
    iraf.fitskypars.smaxiter = 10.0
    iraf.fitskypars.sloclip = 0.0
    iraf.fitskypars.shiclip = 0.0
    iraf.fitskypars.snreject = 50.0
    iraf.fitskypars.sloreject = 3.0
    iraf.fitskypars.shireject = 3.0
    iraf.fitskypars.khist = 3.0
    iraf.fitskypars.binsize = 0.1
    iraf.fitskypars.smooth = False
    iraf.fitskypars.rgrow = 0.0
    iraf.fitskypars.mksky = False

    # iraf.digiphot.apphot.photpars : 
    iraf.unlearn(iraf.photpars)
    iraf.photpars.weighting = 'constant'
    iraf.photpars.apertures = '2,3,4,5'
    iraf.photpars.zmag = ZPT[filt]
    iraf.photpars.mkapert = False

    iraf.unlearn(iraf.phot)
    iraf.gflush(); iraf.gflush()
    iraf.flpr(); iraf.flpr(); iraf.flpr()

    photparams={
        'interac':False,
        'radplot':False,
        }
    outputstuff = iraf.phot(image=imfile, skyfile='',coords=coofile, output=magfile, 
                            verify=False, verbose=True, Stdout=1,**photparams)
    sourcelist = apphotOutput( magfile ) 
    return( sourcelist )

    #omatrix=[]
    #for line in outputstuff:
    #    rawmags = array( line.split()[4:8], dtype=float )
    #    correctedmags = rawmags - APCOR[filt] 
    #    omatrix.append( correctedmags )
    #return(array(omatrix,dtype=float) )


class apphotOutput( object ) :
    
    def __init__( self, magfile ):
        """ read in the apphot output from magfile """
        from numpy import append 
        self.sources = []

        fin = open( magfile, 'r' )
        linelist = fin.readlines()
        fin.close()

        i = 0
        while i<len(linelist) : 
            line = linelist[i]
            if line.startswith('#') : 
                if line.startswith('#K FILTER') : 
                    self.filter = line.split('=')[1].split()[0].strip('"')
                i+=1
                continue
            if len(line.strip())==0 : 
                i+=1
                continue

            # start of a new measured object
            photobj = source(filter=self.filter)
            photobj.id = getdatum(line[23:].split()[2])

            # read in the block of data lines
            for datline in linelist[i+1:] : 
                if not datline.startswith(' '): break
                photobj.datlines.append( datline )
            i+= len(photobj.datlines) + 1

            # dat line 0: position info
            dat = photobj.datlines[0].split()
            photobj.x = float( dat[0] )
            photobj.y = float( dat[1] )
            photobj.xshift = float( dat[2] )
            photobj.yshift = float( dat[3] )
            #photobj.xerr = float( dat[4] )
            #photobj.yerr = float( dat[5] )

            # dat line 1: sky info
            dat = photobj.datlines[1].split()
            photobj.msky = getdatum( dat[0] )
            photobj.skystdev = getdatum( dat[1] )
            photobj.skyerr = getdatum( dat[5] )
            photobj.skyerrstring = dat[6]

            # dat lines 3+: mag info
            for j in range( 3,len(photobj.datlines) ) : 
                dat = photobj.parsedatline( j )
                if dat[4] == 'INDEF' : continue
                apcorrection = getApCorrection( self.filter, float(dat[0]) )
                photobj.aperture      = append( photobj.aperture, getdatum( dat[0] ) )
                photobj.sum           = append( photobj.sum, getdatum( dat[1] ) )
                photobj.area          = append( photobj.area, getdatum( dat[2] ) )
                photobj.flux          = append( photobj.flux, getdatum( dat[3] ) )
                photobj.rawmag        = append( photobj.rawmag, getdatum( dat[4] ) )
                photobj.mag           = append( photobj.mag, getdatum( dat[4] ) - apcorrection )
                photobj.magerr        = append( photobj.magerr, getdatum( dat[5] ) )
                photobj.photerr       = append( photobj.photerr, getdatum( dat[6] ) )
                photobj.photerrstring = append( photobj.photerrstring, dat[7] )
            
            self.sources.append( photobj )
            continue

    def mags( self, aperture ) : 
        """ collect a vector of mags for the given aperture """
        maglist = []
        for i in range( len( self.sources ) ) : 
            photobj = self.sources[i]
            try : 
                iaper = photobj.aperture.tolist().index( aperture ) 
                maglist.append( photobj.mag[iaper] )
            except : 
                maglist.append( 99.0 )

        return( array( maglist ) )

    def magerrs( self, aperture ) : 
        """ collect a vector of mag errs for the given aperture """
        magerrlist = []
        for i in range( len( self.sources ) ) : 
            photobj = self.sources[i]
            try : 
                iaper = photobj.aperture.tolist().index( aperture ) 
                magerrlist.append( photobj.magerr[iaper] )
            except : 
                magerrlist.append( 99.0 )
        return( array( magerrlist ) )
            

class source( object ) :
    def __init__( self, filter='' ) :
        """ initialize a measured object """
        self.filter=filter
        self.datlines = []
        self.aperture      =  array([],dtype=float)
        self.sum           =  array([],dtype=float)
        self.area          =  array([],dtype=float)
        self.flux          =  array([],dtype=float)
        self.rawmag        =  array([],dtype=float)
        self.mag           =  array([],dtype=float)
        self.magerr        =  array([],dtype=float)
        self.photerr       =  array([],dtype=int)
        self.photerrstring =  array([],dtype=str)

        self.x = 0
        self.y = 0
        self.image = ''
        self.id = None

    def parsedatline( self, index ) :
        """ parse an apphot.phot data line into 
        useful data, guarding against poor iraf formatting"""
        datlist = []
        if len(self.datlines)<index : return( datlist )

        try : 
            datlist.append( getdatum(self.datlines[index][:12]) )    # RAPERT
            datlist.append( getdatum(self.datlines[index][12:26]) )  # SUM
            datlist.append( getdatum(self.datlines[index][26:37]) )  # AREA
            datlist.append( getdatum(self.datlines[index][37:51]) )  # FLUX

            if datlist[-1] > 0 : 
                try : 
                    mag = getdatum(self.datlines[index][51:58])
                    merr = getdatum(self.datlines[index][58:64])
                except : 
                    mag = 99.0
                    merr = 99.0
            else : 
                mag = 99.0
                merr = 99.0

            datlist.append( mag )  # MAG
            datlist.append( merr )  # MERR

            datlist.append( getdatum(self.datlines[index][64:69]) )    # PIERR
            datlist.append( self.datlines[index][69:78] )         # PERROR
        except : 
            return( [] )

        return( datlist ) 


def getdatum( datum ) :
    """ parse an apphot datum
    checks for INDEF, and decides if it is an int or float
    """
    if type(datum)==str : 
        if datum == 'INDEF' : datum = 0
        elif datum.find('.')<0 : 
            try :  datum = int(datum)
            except : pass
        if type(datum) == str : 
            try :  datum = float(datum)
            except : pass

    return( datum ) 
        
