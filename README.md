![AOKOI2](https://github.com/user-attachments/assets/15bcc2a3-0f07-481d-92cc-0dd4b14f913d)

AOKOI is an orbit determination software for minor planets and comets. You provide an obs80 file (astrometry in MPC's 80-column format), it calculates an orbit and plots it for you.
It is currently more aimed at short-arc observations, and multi-opposition orbits with sparse observations may fail to converge to a solution.

If you don't know what the obs80 format is, it is explained here on the MPC website: https://minorplanetcenter.net/iau/info/OpticalObs.html
Here is an example obs80 file from MPC database: https://www.minorplanetcenter.net/tmp2/Ao.txt

To be able to use this program, you will need to put de440.bsp, earth_000101_250316_241218.bpc, naif0012.tls and pck00011.tpc files into the /data folder. You can obtain these
from NAIF website: https://naif.jpl.nasa.gov/naif/data.html
