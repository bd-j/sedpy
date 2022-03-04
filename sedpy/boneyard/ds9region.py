# -- VERY PRELIMINARY!!!! -----

import numpy as np
try:
    from matplotlib import path
except(ImportError):
    pass


region_types = ['circle', 'ellipse', 'polygon', 'box']


def load_regfile(regfile):
    reglist = []
    f = open(regfile, 'r')
    for line in f:
        istype = [r in line for r in region_types]
        if True in istype:
            rtype = region_types[np.where(istype)[0][0]]
            numstring = line[line.find("(")+1:line.find(")")].replace('"', '')
            if rtype == 'ellipse':
                reglist.append(Ellipse(numstring))
            elif rtype == 'circle':
                reglist.append(Circle(numstring))
            elif rtype == 'polygon':
                reglist.append(Polygon(numstring))
    f.close()
    return reglist


def get_image_angle(wcs, at_coord=None):
    """
    :returns angle:
       In radians, degrees East of North
    """
    if at_coord is None:
        lon, lat = wcs.wcs.crval[:2]
    else:
        lon, lat = at_coord
    x, y = world2pix(wcs, lon, lat)
    lat_offset = lat + 1./3600.
    xoff, yoff = world2pix(wcs, lon, lat_offset)
    angle = np.arctan2(xoff-x, yoff-y)
    try:
        return angle[0]
    except(TypeError, IndexError):
        return angle


def world2pix(wcs, lon, lat, **extras):
    """Wrap the wcs world2pix method to be less finicky and take our defaults.
    """
    try:
        x, y = wcs.wcs_world2pix(lon, lat, 0)
    except(TypeError):
        x, y, _ = wcs.wcs_world2pix(lon, lat, np.array([0]), 0)
    return x, y


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
            pscale = 3600. * np.sqrt((wcs.wcs.cd[:2, :2]**2).sum(axis=1))
        except (AttributeError):
            pscale = 3600. * np.abs(wcs.wcs.cdelt[:2])
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

    def contains(self, points=None, x=None, y=None, wcs=None, **extras):
        if wcs is not None:
            # Should actually do the plate scale separately for x and y
            # and use great circle distances?
            plate_scale = self.plate_scale(wcs)
            r = self.radius / plate_scale
            cx, cy = world2pix(wcs, self.ra, self.dec)
        else:
            raise ValueError('No WCS given!!!')
        if points is not None:
            x, y = points.T

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
        self.pa = float(bits[4])   #degrees East of North for a

        #swap a and b if poorly defined
        if self.b > self.a:
            self.a, self.b = self.b, self.a
            self.pa = self.pa + 90.

    @property
    def unparse(self):
        return ','.join([str(a) for a in [self.ra, self.dec, self.a, self.b, self.pa]])

    def contains(self, points=None, x=None, y=None, wcs=None, **extras):
        """
        Determine whether pixel locations x,y are within the ellipse.
        Requires that a WCS be given.  Assumes tangent projection, and
        will fail for large ellipses/images. Assumes N is up and E to
        the left.

        x, and y are zero indexed, and can be produced by flattened versions of
        `y, x = np.indices(image)`
        """
        if wcs is not None:
            cx, cy = world2pix(wcs, self.ra, self.dec)
            # should actually do the plate scale separately for x and y
            # and use great circle distances?
            plate_scale = self.plate_scale(wcs)
            a, b = self.a / plate_scale, self.b / plate_scale
            theta_image = get_image_angle(wcs)
        else:
            raise ValueError('No WCS given!!!')
        if points is not None:
            x, y = points.T

        # --- Translate ---
        dx, dy = x - cx, y - cy
        # --- Rotate ---
        #  convert to radians from positive x
        theta = np.deg2rad(self.pa) - theta_image + np.pi
        # construct rotation matrix
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        # We flip dx and dy here because why?
        # it seems to work.... something to do with column-major/row-major
        prime = np.dot(R, np.vstack([dy, dx]))

        r = np.hypot(prime[0, :] / a, prime[1, :] / b)
        #side1 = ((dx * np.cos(np.deg2rad(self.pa-90)) - dy * np.sin(np.deg2rad(self.pa-90))) / a)**2
        #side2 = ((dx * np.sin(np.deg2rad(self.pa-90)) + dy * np.cos(np.deg2rad(self.pa-90))) / b)**2
        return r <= 1


class Polygon(Region):

    shape = 'polygon'

    def parse(self, defstring):
        bits = defstring.split(',')
        bits = np.array([float(b) for b in bits])
        #print( len(bits), bits)
        assert (len(bits) % 2) == 0
        self.n_vertices = len(bits)/2
        self.ra = bits[np.arange(self.n_vertices) * 2]       #degrees
        self.dec = bits[np.arange(self.n_vertices) * 2 + 1]  #degrees

    @property
    def unparse(self):
        return ','.join([str(val) for pair in zip(self.ra, self.dec) for val in pair])

    def contains(self, points=None, x=None, y=None, wcs=None,
                 fast=False, pad=[0, 0], **extras):
        """
        Determine whether pixel locations x,y are within the polygon.
        Requires that a WCS be given.  Uses Path objects from matlib,
        should probably use Shapely, or better yet a real spherical
        coordinates implementation
        """
        if wcs is not None:
            vv = world2pix(wcs, self.ra, self.dec)
            #vv[0] -= 1 #convert to numpy indexing
        else:
            vv = (self.ra, self.dec)

        vertices = [(r, d) for r, d in zip(vv[0], vv[1])]
        vertices += [vertices[0]]
        poly = path.Path(vertices, closed=True)
        if points is None:
            points = np.vstack((x, y)).T

        if fast:
            sub = ((points[:, 0] > vv[0].min()-pad[0]) &
                   (points[:, 0] < vv[0].max()+pad[0]) &
                   (points[:, 1] > vv[1].min()-pad[1]) &
                   (points[:, 1] < vv[1].max()+pad[1]))
            sub_is_in = poly.contains_points(points[sub, :])
            is_in = np.zeros(points.shape[0], dtype=bool)
            is_in[sub] = sub_is_in
        else:
            is_in = poly.contains_points(points)
            #flux = (self.image - background)[is_in].sum()
        return is_in
