import sys, os, glob
import subprocess
import numpy as np
from . import photometer
from astropy import wcs as pywcs
try:
    import astropy.io.fits as pyfits
except(ImportError):
    import pyfits


def run_command(cmd):
    """Open a child process, and return its exit status and stdout"""
    child = subprocess.Popen(cmd, shell =True, stderr = subprocess.PIPE, 
                             stdin=subprocess.PIPE, stdout = subprocess.PIPE)
    out = [s for s in child.stdout]
    w = child.wait()
    return os.WEXITSTATUS(w), out


class Image(object):

    def __init__(self, imagename):
        self.image_name = imagename
        self.image = pyfits.getdata(imagename)
        self.header = pyfits.getheader(imagename)
        self.wcs = pywcs.WCS(self.header)
        # a rough plate scale in arcseconds
        try:
            self.plate_scale = 3600. *np.sqrt((self.wcs.wcs.cd**2).sum(axis = 1))
        except (AttributeError):
             self.plate_scale = 3600. * np.abs(self.wcs.wcs.cdelt)
             
        self.mask = np.ones_like(self.image, dtype = bool)
        ny, nx = self.mask.shape  
        self.y, self.x = np.meshgrid(np.arange(ny), np.arange(nx))
        
    def sextract(self, catalog_name = 'junk.dat', **kwargs):
        """
        Run SExtractor on the image.  One should already have a default.param
        and a default.sex file laying around with most of the options set.  Howewver, any of the
        options can be set by supplying them as keywords, using a similar format
        to the commandline options.  Note that errors may be raised if the number  of
        apertures does not match the number specified in default.params.

        N.B command line changing of phot_apertures NOT WORKING
        """
        f   = self.image_name
        cmd = ['sex',f,'-CATALOG_NAME ' + catalog_name]
        for k, v in kwargs.iteritems():
            cmd[2] += ' -' + k.upper() + ' ' + v
        self.sex_command = cmd
        stat, out = run_command(' '.join(cmd))
        
        dat = read_sexcat(catalog_name)
        return dat

    def point_source_phot(self, sourcelist, shape, skypars):
        """
        Do simple aperture photomtry at a list of locations.  The list
        is ideally the record array returned by read_sexcat.
        """
        ## Initialize aperture (and background annulus)
        ap = photometer.Circular()
        ap.background = photometer.Annulus(bgtype = skypars['bgtype'])
        #ap.background = photometer.ZeroSky()
        flux = np.zeros(len(sourcelist))
        for i,s in enumerate(sourcelist):
            shape['xcen'] = s['x_image']-0.5
            shape['ycen'] = s['y_image']-0.5
            skypars['xcen'] = shape['xcen']
            skypars['ycen'] = shape['ycen']
            flux[i], dummy = ap.measure_flux(shape, self.image, skypars = skypars)
        return flux
    
    def area_phot(self, regionlist, use_mask = False):
        """
        Add up all the flux in some (large) regions of the sky defined
        by region objects. Fractional pixels not considered, so make
        your polygon big (or magnify the image and WCS).
        """
        if use_mask:
            mask = self.mask
        else:
            mask = np.ones_like(self.image, dtype = bool)
        for region in regionlist:
            mask = mask & self.region_bool_mask(region)
            
        flux = np.nansum(self.image[mask])
        
        return  flux

    def region_bool_mask(self, region):
        
        ny, nx = self.mask.shape        
        x, y = self.x.flatten(), self.y.flatten()
        inds = region.contains(x, y, wcs = self.wcs)
        grid = inds.reshape((ny,nx)).T
        return grid

    
def read_sexcat(filename):
    """Read a SExtractor results catalog, and return a numpy record array.
    Because astropy Tables are annoying.
    """
    f = open(filename,'r')
    colname = [] #parameter header
    skiprows = 0
    oldpos = 0
    for line in f:
        if line[0] !='#':
            ncol = len(line.split())
            break
        else:
            #parse the column descriptor line
            d = line.split()
            pos, name, units = int(d[1]), d[2].lower(),  d[-1][1:-1]
            #deal with vectors
            off = 1
            while pos > oldpos + off:
                if off == 1:
                    colname[-1] += '_1'
                colname.append(oldname + '_' + str(off+1))
                off += 1
            colname.append(name)
            oldpos, oldname = pos, name
            skiprows +=1
        
    f.close()
    dt = len(colname) * ['<f8']
    dtype = np.dtype(zip(colname, dt))
    dat = np.loadtxt(filename, dtype = np.dtype(zip(colname, dt)), skiprows = skiprows)
    return dat


    
def xylist(sourcelist, outname, physical = False, wcs = None):
    cols = ['alpha_j2000', 'delta_j2000']
    xx, yy = sourcelist[cols[0]], sourcelist[cols[1]]
    out = open(outname, 'w')
    out.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    for x,y in zip(xx, yy):
        out.write('fk5;circle({0},{1},5.0")\n'.format(x,y))
    out.close()
    
