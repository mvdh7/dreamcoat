import numpy as np
import vptree
import dreamcoat as dc


lat_range = (50, 55)
lon_range = (-5, 10)
coords = dc.maps.coastline_coords(lat_range, lon_range, resolution="110m")


def test_geodesic_distance():
    """Do we calculate the same distances as http://edwilliams.org/gccalc.htm ?"""
    dist0 = dc.maps.geodesic_distance((0, 1), (0, 0))
    assert isinstance(dist0, float)
    assert np.isclose(dist0, 110.57438855790893)
    dist1 = dc.maps.geodesic_distance((1, 0), (0, 0))
    assert isinstance(dist1, float)
    assert np.isclose(dist1, 111.31949077920639)
    dist2 = dc.maps.geodesic_distance((12, -23), (-135.8, 86.4))
    assert isinstance(dist2, float)
    assert np.isclose(dist2, 12885.15245722158)


def test_coastline_coords():
    """Do we return coordinates of the coastline(s) in the expected format?"""
    assert isinstance(coords, list)
    assert len(coords) > 0
    for coord in coords:
        assert isinstance(coord, tuple)
        assert len(coord) == 2
        assert isinstance(coord[0], float)
        assert isinstance(coord[1], float)


def test_build_vptree():
    """Can we build and use a vantage-point tree for the coastline?"""
    vpt = dc.maps.build_vptree(lat_range, lon_range, resolution="110m")
    assert isinstance(vpt, vptree.VPTree)
    nn0 = vpt.get_nearest_neighbor(coords[0])
    assert isinstance(nn0[0], float)
    assert nn0[0] == 0
    assert np.all(np.isclose(coords[0], nn0[1]))


# test_geodesic_distance()
# test_coastline_coords()
# test_build_vptree()
