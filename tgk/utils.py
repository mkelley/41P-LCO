# Licensed under a MIT style license - see LICENSE
"""utils - Various common utilities."""

import numpy as np

class UnableToCenter(Exception):
    pass

def anphot(im, yx, rap, subsample=4, squeeze=True):
    """Simple annular aperture photometry.

    Pixels may be sub-sampled, and sub-sampling may be CPU and memory
    intensive.

    Parameters
    ----------
    im : array or array of arrays
      An image, cube, or array of images on which to measure
      photometry.  For data cubes, the first axis iterates over the
      images.  All images must have the same shape.
    yx : array
      The `y, x` center of the aperture(s), or an Nx2 length array of
      centers. [pixels]
    rap : float or array
      Aperture radii.  The inner-most aperture will be the annulus 0
      to `min(rap)`.  [pixels]
    subsample : int, optional
      The sub-pixel sampling factor.  Set to `<= 1` for no sampling.
      This will sub-sample the entire image.
    squeeze : bool, optional
      Set to `True` to sqeeze single length dimensions out of the
      results.

    Returns
    -------
    n : ndarray
      The number of pixels per annular bin, either shape `(len(rap),)`
      or `(len(yx), len(rap))`.
    f : ndarray
      The annular photometry.  The shape will be one of:
        `(len(rap),)`
        `(len(yx), len(rap))`
        `(len(im), len(yx), len(rap))`

    """

    _im = np.array(im)
    assert _im.ndim in [2, 3], ("Only images, data cubes, or tuples/lists"
                                " of images are allowed.")
    if _im.ndim == 2:
        _im = _im.reshape((1,) + _im.shape)

    yx = np.array(yx, float)
    assert yx.ndim in [1, 2], "yx must be one or two dimensional."
    if yx.ndim == 1:
        assert yx.shape[0] == 2, "yx must have length 2."
        yx = yx.reshape((1, 2))
    else:
        assert yx.shape[1] == 2, "Second axis of yx must have length 2."

    if not np.iterable(rap):
        rap = np.array([rap])

    if subsample > 1:
        _im = np.array([rebin(x, subsample, flux=True) for x in _im])
        yx = yx * subsample + (subsample - 1) / 2.0

    sz = _im.shape[-2:]

    # flatten all arrays for digitize
    N = _im.shape[0]
    M = _im.shape[1] * _im.shape[2]
    _im = _im.flatten().reshape((N, M))

    n = np.zeros((len(yx), len(rap)))
    f = np.zeros((N, len(yx), len(rap)))

    # annular photometry via histograms, histogram via digitize() to
    # save CPU time when mulitiple images are passed
    for i in range(len(yx)):
        r = rarray(sz, yx=yx[i], subsample=10) / float(subsample)
        bins = np.digitize(r.flatten(), rap)
        for j in range(len(rap)):
            ap = np.flatnonzero(bins == j)
            f[:, i, j] = np.sum(_im[:, ap], 1)
            n[i, j] = len(ap)

    n /= float(subsample**2)
    if squeeze:
        return n.squeeze(), f.squeeze()
    else:
        return n, f

def apphot(im, yx, rap, subsample=4, **kwargs):
    """Simple aperture photometry.

    Pixels may be sub-sampled, and sub-sampling may be CPU and memory
    intensive.

    Parameters
    ----------
    im : array or array of arrays
      An image, cube, or tuple of images on which to measure
      photometry.  For data cubes, the first axis iterates over the
      images.  All images must have the same shape.
    yx : array
      The `y, x` center of the aperture(s), or an Nx2 length array of
      centers. [pixels]
    rap : float or array
      Aperture radii.  [pixels]
    subsample : int, optional
      The sub-pixel sampling factor.  Set to `<= 1` for no sampling.
      This will sub-sample the entire image.
    **kwargs
      Any `anphot` keyword argument.

    Returns
    -------
    n : ndarray
      The number of pixels per aperture, either shape `(len(rap),)` or
      `(len(yx), len(rap))`.
    f : ndarray
      The annular photometry.  The shape will be one of:
        `(len(rap),)`
        `(len(yx), len(rap))`
        `(len(im), len(yx), len(rap))`  (for multiple images)

    """
    n, f = anphot(im, yx, rap, subsample=subsample, **kwargs)
    if np.size(rap) > 1:
        return n.cumsum(-1), f.cumsum(-1)
    else:
        return n, f

def cutout(yx, half_size, shape=None):
    """Return a slice to cut out a subarray from an array.

    Parameters
    ----------
    yx : array of ints
      The center of the cutout.
    half_size : int or array of ints
      The half_size of the array.  The cut out will have shape `2 *
      half_size + 1`.
    shape : tuple, optional
      If provided, then the slice will not extend beyond the lengths
      of the axes.
    
    Returns
    -------
    s : slice

    """

    if shape is None:
        shape = (inf, inf)
    if not np.iterable(half_size):
        half_size = (half_size, half_size)

    s = np.s_[max(yx[0] - half_size[0], 0):
              min(yx[0] + half_size[0] + 1, shape[0]),
              max(yx[1] - half_size[1], 0):
              min(yx[1] + half_size[1] + 1, shape[1])]
    return s

def gaussian(x, mu, sigma):
    """A normalized Gaussian curve.

    Parameters
    ----------
    x : array
      Dependent variable.
    mu : float
      Position of the peak.
    sigma : float
      Width of the curve (sqrt(variance)).

    Returns
    -------
    G : ndarray
      The Gaussian function.

    """
    return (np.exp(-(x - mu)**2 / 2.0 / sigma**2) /
            np.sqrt(2.0 * np.pi) / sigma)

def gcentroid(im, yx=None, box=None, niter=1, dim=2, shrink=True, silent=True):
    """Centroid (x-/y-cut Gaussian fit) of an image.

    The best-fit should be bounded within `box`.

    Parameters
    ----------
    im : ndarray
      A 2D image on which to centroid.
    yx : float array, optional
      `(y, x)` guess for the centroid.  The default guess is the image
      peak.
    box : array, optional
      Specify the size of the box over which to compute the centroid.
      This may be an integer, or an array (width, height).  The
      default is to use the whole image.
    niter : int, optional
      When box is not None, iterate niter times.
    dim : int, optional
      Set to 1 to fit two 1D Gaussians, or 2 to fit a 2D Gaussian.
    shrink : bool, optional
      When iterating, decrease the box size by sqrt(2) each time.
    silent : bool, optional
      Suppress any print commands.

    Returns
    -------
    cyx : ndarray
      The computed center.  The lower-left corner of a pixel is -0.5,
      -0.5.

    """

    from photutils.centroids import centroid_1dg, centroid_2dg

    assert dim in [1, 2], '`dim` must be one of [1, 2]'
    if dim == 1:
        centroid_func = centroid_1dg
    else:
        centroid_func = centroid_2dg

    if yx is None:
        yx = np.array(np.unravel_index(np.nanargmax(im), im.shape), float)

    # the array index location of yx
    iyx = np.round(yx).astype(int)

    if box is None:
        box = np.array((im.shape[0], im.shape[1]))
    elif np.size(box) == 1:
        box = np.array((box, box)).reshape((2,))
    else:
        box = np.array(box)

    halfbox = box // 2

    yr = [max(iyx[0] - halfbox[0], 0),
          min(iyx[0] + halfbox[0] + 1, im.shape[0] - 1)]
    xr = [max(iyx[1] - halfbox[1], 0),
          min(iyx[1] + halfbox[1] + 1, im.shape[1] - 1)]
    ap = (slice(*yr), slice(*xr))

    try:
        cyx = centroid_func(im[ap])[::-1]
    except ValueError as e:
        raise UnableToCenter(str(e))

    # convert from aperture coords to image coords
    cyx = cyx + np.r_[yr[0], xr[0]]

    if niter > 1:
        if shrink:
            box = (halfbox * 2 / np.sqrt(2)).round().astype(int)
            box = box + (box % 2 - 1)  # keep it odd

        if not silent:
            print("y, x = {0[0]:.1f}, {0[1]:.1f}, next box size = {1}".format(
                cyx, str(box)))

        if max(box) < 2:
            if not silent:
                print(" - Box size too small -")
            return cyx

        return gcentroid(im, yx=cyx, box=box, niter=niter-1, dim=dim,
                         shrink=shrink, silent=silent)
    else:
        if not silent:
            print("y, x = {0[0]:.1f}, {0[1]:.1f}".format(cyx))

    return cyx

def rarray(shape, yx=None, subsample=0, dtype=float):
    """Array of distances from a point.

    Parameters
    ----------
    shape : array
      The shape of the resulting array `(y, x)`.
    yx : array, optional
      The center of the array `(y, x)`.  If set to None, then the
      center is `(shape - 1.0) / 2.0` (floating point arithmetic).
      Integer values refer to the center of the pixel.
    subsample : int, optional
      Set to `>1` to sub-pixel sample the core of the array.
    dtype : np.dtype or similar, optional
      Set to the data type of the resulting array.

    Returns
    -------
    r : ndarray
      The array of radial values.

    """

    if yx is None:
        yx = (np.array(shape) - 1.0) / 2.0
    else:
        yx = np.array(yx)

    y = yarray(shape, dtype=dtype) - yx[-2]
    x = xarray(shape, dtype=dtype) - yx[-1]
    r = (np.sqrt(x**2 + y**2)).astype(dtype)

    if subsample > 0:
        r = refine_center(rarray, r, yx, 11, subsample, scale=-1, dtype=dtype)

    return r

def rebin(a, factor, flux=False, trim=False):
    """Rebin a 1, 2, or 3 dimensional array by integer amounts.

    Parameters
    ----------
    a : ndarray
      Image to rebin.
    factor : int
      Rebin factor.  Set factor < 0 to shrink the image.  When
      shrinking, all axes of a must be an integer multiple of factor
      (see trim).
    flux : bool
      Set to True to preserve flux, otherwise, surface brightness is
      preserved.
    trim : bool
      Set to True to automatically trim the shape of the input array
      to be an integer multiple of factor.

    Returns
    -------
    b : ndarray
      The rebinned array.

    Notes
    -----
    By default, the function preserves surface brightness, not flux.

    """

    def mini(a, factor):
        b = a[::-factor]
        for i in range(-factor - 1):
            b += a[(i + 1)::-factor]
        if not flux:
            b /= float(-factor)
        return b

    def magni(a, factor):
        s = np.array(a.shape)
        s[0] *= factor
        b = np.zeros(s)
        for i in range(factor):
            b[i::factor] = a
        if flux:
            b /= float(factor)
        return b

    if factor == 1:
        # done!
        return a

    _a = a.copy()
    if factor < 0:
        for i in range(_a.ndim):
            if trim:
                r = _a.shape[i] % abs(factor)
                if r != 0:
                    _a = np.rollaxis(np.rollaxis(_a, i)[:-r], 0, i + 1)

            assert (_a.shape[i] % factor) == 0, (
                "Axis {0} must be an integer multiple of "
                "the minification factor.".format(i))
        f = mini
    else:
        f = magni

    b = f(_a, factor)
    for i in range(len(_a.shape) - 1):
        c = f(np.rollaxis(b, i + 1), factor)
        b = np.rollaxis(c, 0, i + 2)

    return b

def refine_center(func, im, yx, N, subsample, scale=0, **kwargs):
    """Subsample an array generating function near the center.

    Parameters
    ----------
    func : function
      The array generating function, e.g., `rarray`.  The first
      argument of the function is the shape of the result.  Function
      keywords must include `yx` and `subsample`.  See `rarray` for an
      example.
    im : array
      The initially generated array, the center pixels of which will
      be replaced.
    yx : array
      The function's origin.
    N : int
      The center `(N, N)` pixels will be refined.  Best is an odd
      value.
    subsample : int
      The sub-pixel sampling factor.
    scale : float
      Scale the refined center by `subsample**scale`.  The refined
      area will be generated at high resolution, and rebinned
      (averaging) to a resolution of 1 pixel.  When, for example, a
      function has dimensions of length, the rebinned array needs to
      be scaled by `subsample**-1`.
    **kwargs
      Additional keyword arguments for `func`.

    Returns
    -------
    refined : ndarray
      The refined array.

    Notes
    -----
    Numpy rounds to even numbers, which causes an error when working
    with an origin at half-pixel steps (i.e., 50.5 rounds to 50, but
    should be 51).  `refine_center` adds a small value to `yx` to
    mitigate this issue.  If better than 1e-4 pixel precision is
    needed, this function may not work for you.

    """

    yx_f, yx_i = np.modf(yx)  # fractional and integer parts
    yx_i = yx_i.astype(int)

    # where is the center of the NxN region?
    yx_N = (np.ones(2, int) * N - 1.0) / 2.0 + yx_f

    # subsample the NxN region, where is the center in that?
    shape_s = np.ones(2, int) * N * subsample
    yx_s = subsample * (yx_N + 0.5) - 0.5    
    
    # generate the subsampled center
    refined_s = func(shape_s, yx=yx_s, subsample=0, **kwargs)
    refined_s = rebin(refined_s, -subsample, flux=False)
    refined_s *= subsample**scale

    # The region to be refined: xi, yi
    yi, xi = np.indices((N, N)) - N // 2
    yi += yx_i[0]
    xi += yx_i[1]

    # insert into the result
    refined = im.copy()
    i = (yi >= 0) * (xi >= 0) * (yi < im.shape[0]) * (xi < im.shape[1])
    if np.any(i):
        refined[yi[i], xi[i]] = refined_s[i]

    return refined

def xarray(shape, yx=[0, 0], dtype=int):
    """Array of x values.

    Parameters
    ----------
    shape : tuple of int
      The shape of the resulting array, e.g., (y, x).
    yx : tuple of float
      Offset the array to align with this y, x center.
    dtype : numpy.dtype, optional
      The data type of the result.

    Returns
    -------
    x : ndarray
      An array of x values.

    """

    y, x = np.indices(shape, dtype)[-2:]
    y -= yx[0]
    x -= yx[1]

    return x

def yarray(shape, yx=[0, 0], dtype=int):
    """Array of y values.

    Parameters
    ----------
    shape : tuple of int
      The shape of the resulting array, e.g., (y, x).
    yx : tuple of float
      Offset the array to align with this y, x center.
    dtype : numpy.dtype, optional
      The data type of the result.

    Returns
    -------
    y : ndarray
      An array of y values.

    """

    y, x = np.indices(shape, dtype)[-2:]
    y -= yx[0]
    x -= yx[1]

    return y

