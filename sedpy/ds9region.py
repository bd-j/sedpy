### VERY PRELIMINARY!!!! #####

import numpy as np
try:
    from matplotlib import path
except(ImportError):
    pass

region_types = ['circle', 'ellipse', 'polygon', 'box']

def load_regfile(regfile):
    reglist = []
    f = open(regfile,'r')
    for line in f:
        istype = [r in line for r in region_types]
        if True in istype:
            rtype = region_types[np.where(istype)[0][0]]
            numstring = line[line.find("(")+1:line.find(")")].replace('"','')
            if rtype == 'ellipse':
                reglist.append( Ellipse(numstring) )
            elif rtype == 'circle':
                reglist.append( Circle(numstring) )
            elif rtype =='polygon':
                reglist.append( Polygon(numstring) )
    f.close()
    return reglist
    
class Region(object):
    
    def __init__(self, defstring):
        self.parse(defstring)

    def parse(self, defstring):
        raise(NotImplementedError)

    @property
    def unparse(self):
        raise(NotImplementedError)        
            
    def plate_scale(self, wcs):
        try:
            pscale = 3600. *np.sqrt((wcs.wcs.cd**2).sum(axis = 1))
        except (AttributeError):
            pscale = 3600. * np.abs(wcs.wcs.cdelt)
        return pscale.mean()
    
    def print_to(self, fileobj=None, color='green', label=''):
        fstring = 'fk5;{0}({1}) # color={2} text={{{3}}}\n'
        line = fstring.format(self.shape, self.unparse, color, label)
        if fileobj is None:
            print(line)
        else:
            fileobj.write(line)

class Circle(Region):

    shape = 'circle'

    def parse(self, defstring):
        bits = defstring.split(',')
        self.ra = float(bits[0])
        self.dec = float(bits[1])
        self.radius = float(bits[2])

    @property
    def unparse(self):
        return ','.join([str(a) for a in [self.ra, self.dec, self.radius]])
        
    def contains(self, x, y, wcs = None):
        if wcs is not None:
            #should actually do the plate scale separately for x and y
            # and use great circle distances?
            plate_scale = self.plate_scale(wcs)
            r = self.radius/plate_scale
            cx, cy = wcs.wcs_world2pix(self.ra, self.dec, 0)
        else:
            raise ValueError('No WCS given!!!')
        dx, dy = x-cx, y-cy
        return (dx**2 + dy**2) <= r**2
        
class Ellipse(Region):

    shape = 'ellipse'

    def parse(self, defstring):
        bits = defstring.split(',')
        self.ra = float(bits[0])   #degrees 
        self.dec = float(bits[1])  #degrees
        self.a = float(bits[2])    #arcsec
        self.b = float(bits[3])    #arcsec
        self.pa = float(bits[4])   #degrees East of North
        
        #swap a and b if poorly defined
        if self.b > self.a:
            self.a, self.b = self.b, self.a
            self.pa = self.pa + 90.
    @property
    def unparse(self):
        return ','.join([str(a) for a in [self.ra, self.dec, self.a, self.b, self.pa]])

    def contains(self, x, y, wcs = None):
        """
        Determine whether pixel locations x,y are within the ellipse.
        Requires that a WCS be given.  Assumes tangent projection, and
        will fail for large ellipses/images. Assumes N is up and E to
        the left.
        """
        if wcs is not None:
            #should actually do the plate scale separately for x and y
            # and use great circle distances?
            plate_scale = self.plate_scale(wcs)
            a, b = self.a/plate_scale, self.b/plate_scale
            cx, cy = wcs.wcs_world2pix(self.ra, self.dec, 0)
        else:
            raise ValueError('No WCS given!!!')
        dx, dy = x-cx, y-cy
        side1 = ((dx * np.cos(np.deg2rad(self.pa-90)) - dy * np.sin(np.deg2rad(self.pa-90))) / a)**2
        side2 = ((dx * np.sin(np.deg2rad(self.pa-90)) + dy * np.cos(np.deg2rad(self.pa-90))) / b)**2
        return np.sqrt(side1 + side2) <= 1
            
class Polygon(Region):

    shape = 'polygon'

    def parse(self, defstring):
        bits = defstring.split(',')
        bits = np.array([float(b) for b in bits])
        #print( len(bits), bits)
        assert( (len(bits) % 2) == 0)
        self.n_vertices = len(bits)/2
        self.ra = bits[np.arange(self.n_vertices) * 2]       #degrees 
        self.dec = bits[np.arange(self.n_vertices) * 2 + 1]  #degrees 

    @property
    def unparse(self):
        return ','.join([ str(val) for pair in zip(self.ra, self.dec) for val in pair])
        
    def contains(self, x, y, wcs = None):
        """
        Determine whether pixel locations x,y are within the polygon.
        Requires that a WCS be given.  Uses Path objects from matlib,
        should probably use Shapely, or better yet a real spherical
        coordinates implementation
        """
        if wcs is not None:
            vv = wcs.wcs_world2pix(self.ra, self.dec, 1)
            vv[0] -= 1 #convert to numpy indexing
        else:
            vv = (self.ra, self.dec)
            
        vertices = [(r,d) for r,d in zip(vv[0], vv[1])]
        vertices += [vertices[0]]
        poly = path.Path(vertices, closed=True)
        points = np.vstack((x,y)).T
        is_in = poly.contains_points(points)
        #flux = (self.image - background)[is_in].sum()
        return is_in
