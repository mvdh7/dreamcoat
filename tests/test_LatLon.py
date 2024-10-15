import dreamcoat as dc


def test_dd():
    lat = 3.3
    lon = -4.2
    ll = dc.LatLon(lat, lon)
    lon_out, lat_out = ll.to_dd()
    assert lon_out == lon
    assert lat_out == lat


def test_ddm():
    ll_dashes = dc.LatLon("67-12.55-N", "155-8.3-W")
    ll_fancy = dc.LatLon("67°12.55'N", "155°8.3'W")
    assert ll_dashes.latitude_dd == ll_fancy.latitude_dd
    assert ll_dashes.longitude_dd == ll_fancy.longitude_dd
    ll_str = ll_dashes.to_ddm()
    assert isinstance(ll_str, str)


def test_wrap():
    ll_180 = dc.LatLon(longitude=-90)
    ll_360 = dc.LatLon(longitude=270)
    assert ll_180.longitude_dd == ll_360.longitude_dd == -90


# test_dd()
# test_ddm()
# test_wrap()
