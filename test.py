from eos3dmpl import *

# orbit data generation test
eos = eos_core.Eos()
eos.orbit(a=10000,e=0)
data = eos.kepler(data=True,vis=None)

# TLE parsing test
eos.parseTLE()

# pycurl test
c = eos_core.pycurl.Curl()
c.setopt(c.URL, 'https://www.google.com')
with eos_core.BytesIO() as e:
    c.setopt(c.WRITEFUNCTION, e.write)
    c.setopt(eos_core.pycurl.SSL_VERIFYPEER, 0)
    c.perform()
    c.close()

# SGP4 propagation test

l1 = "1 25544U 98067A   18172.92063794 +.00001967 +00000-0 +36850-4 0  9998"
l2 = "2 25544 051.6423 355.1116 0002915 209.8982 300.4127 15.54210114119259"
sat = eos_core.twoline2rv(l1, l2, eos_core.wgs84)
pos, vel = sat.propagate(2018, 6, 22, 21, 0, 0)
