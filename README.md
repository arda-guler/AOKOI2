![AOKOI2](https://github.com/user-attachments/assets/15bcc2a3-0f07-481d-92cc-0dd4b14f913d)

## About
AOKOI is an orbit determination software for minor planets and comets. You provide an obs80 file (astrometry in MPC's 80-column format), it calculates an orbit and plots it for you.
It is currently more aimed at short-arc observations, and multi-opposition orbits with sparse observations may fail to converge to a solution.

If you don't know what the obs80 format is, it is explained here on the MPC website: https://minorplanetcenter.net/iau/info/OpticalObs.html
Here is an example obs80 file from MPC database: https://www.minorplanetcenter.net/tmp2/Ao.txt

To be able to use this program, you will need to put de440.bsp, earth_000101_250316_241218.bpc, naif0012.tls and pck00011.tpc files into the /data folder. You can obtain these
from NAIF website: https://naif.jpl.nasa.gov/naif/data.html

## Example Result Plots

<p align="center">
  <img src="https://github.com/user-attachments/assets/03f04173-38cc-4fbb-83d5-4a3ff57b6028" />
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/81fb08a7-a1d1-4ebf-87b6-5b3252ee2b2f" />
</p>

## License
AOKOI

Copyright (C) 2025  H. A. Guler

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
<https://www.gnu.org/licenses/>.
