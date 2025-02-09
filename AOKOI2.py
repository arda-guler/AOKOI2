import numpy as np
import spiceypy as spice
from datetime import datetime, timedelta
import re
import sys
import matplotlib.pyplot as plt
import subprocess
import os
import mplcursors
import json

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

title = "AOKOI"
version = "2.0.0"
help_text =\
"""AOKOI is an orbit determination software for minor planets and comets. You provide an obs80 file (astrometry in MPC's 80-column format), it calculates an orbit and plots it for you. """ +\
"""It is currently more aimed at short-arc observations, and multi-opposition orbits with sparse observations may fail to converge to a solution.\n\n"""+\
"""If you don't know what the obs80 format is, it is explained here on the MPC website: https://minorplanetcenter.net/iau/info/OpticalObs.html\n"""+\
"""Here is an example obs80 file from MPC database: https://www.minorplanetcenter.net/tmp2/Ao.txt

To be able to use this program, you will need to put de440.bsp, earth_000101_250316_241218.bpc, naif0012.tls and pck00011.tpc files into the /data folder. You can obtain these """+\
"""from NAIF website: https://naif.jpl.nasa.gov/naif/data.html
"""
authors = "Authors: H. A. Guler"
license_text=\
"""AOKOI
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
"""

# Constants
km_per_AU = 149597870.7
AU = 149597870.7
km_per_mAU = 149597870.7 * 1e-3
s_per_d = 86400
day = 86400
arcsecond = 0.00027778 # deg
observatories = {}

## CLASS DEFINITIONS
class MainBody: # for the Sun and planets
    def __init__(self, name, pos, GM):
        self.name = name
        self.pos = pos
        self.GM = GM # a.k.a. mu, standard gravitational parameter

class MP: # minor planets, comets
    def __init__(self, des, pos, vel):
        self.des = des
        self.pos = pos
        self.vel = vel

class Obs: # MPC 80-column observations
    def __init__(self, obs_str, debug=False):
        if type(obs_str) == str:
            self.obs_str = obs_str
            
            packed_perm = obs_str[0:5]
            packed_prov = obs_str[5:12]
            date_str = obs_str[15:31]
            RA_str = obs_str[32:43]
            DEC_str = obs_str[44:55]
            mag_str = obs_str[65:69]
            obscode_str = obs_str[77:80]

            # date handling
            obs_date_parts = date_str.split()
            year = int(obs_date_parts[0])
            month = int(obs_date_parts[1])
            day = int(float(obs_date_parts[2]))

            decimal_day = float(obs_date_parts[2]) - day
            day_seconds = 86400
            decimal_secs = day_seconds * decimal_day

            obs_datetime = datetime(year, month, day)
            obs_datetime = obs_datetime + timedelta(seconds=decimal_secs)

            # RA - DEC handling
            hour2deg = 360/24
            minute2deg = 360/(24*60)
            second2deg = 360/(24*60*60)

            RA_parts = RA_str.split()
            RA_deg = float(RA_parts[0]) * hour2deg + float(RA_parts[1]) * minute2deg + float(RA_parts[2]) * second2deg

            if obs_str[44] == "-":
                DEC_sign = -1
            else:
                DEC_sign = 1

            DEC_parts = DEC_str.split()
            DEC_deg = (abs(float(DEC_parts[0])) + float(DEC_parts[1]) / 60 + float(DEC_parts[2]) / 3600) * DEC_sign

            # mag handling
            try:
                mag = float(mag_str)
            except ValueError:
                mag = 0

            self.perm = packed_perm
            self.prov = packed_prov
            self.date = obs_datetime
            self.RA = RA_deg
            self.DEC = DEC_deg
            self.mag = mag
            self.obs_code = obscode_str
            
        elif type(obs_str) == list:
            self.RA = obs_str[0]
            self.DEC = obs_str[1]
            self.date = obs_str[2]

    def __str__(self):
        return self.obs_str

    def __repr__(self):
        return str(self)

class PlotPlanet:
    def __init__(self, name, a, e, i, OMG, o, L):
        self.name = name
        self.a = a
        self.e = e
        self.i = i
        self.OMG = OMG # big omega
        self.o = o # small omega with a tilda
        self.L = L

class Observatory:
    def __init__(self, code, lon, cos_phi, sin_phi, name):
        self.code = code
        self.lon = lon
        self.name = name

        # convert MPC's phi and rho to latitude and altitude
        self.lat = np.degrees(np.arcsin(sin_phi))
        self.rho = 6371.0088 * (cos_phi**2 + sin_phi**2)**0.5

    def __repr__(self):
        return (f"Observatory(code={self.code}, name='{self.name}', "
                f"lon={self.lon:.4f}, lat={self.lat:.4f}, rho={self.rho:.4f})")

def performPairAnalysis(knwon_pair, oc, pix, res):
    known_obses = text1.split("\n")
    o1 = Obs(known_pair[0])
    o2 = Obs(known_pair[1])

    pair = ObsPair(o1, o2, pix, res)
    check = pair.check_obs(oc)

def spherical2cartezian(d, RA, DEC):
    x = d * np.cos(DEC) * np.cos(RA)
    y = d * np.cos(DEC) * np.sin(RA)
    z = d * np.sin(DEC)

    return np.array([x, y, z])

def cartezian2spherical(vec):
    d = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
    RA = np.arctan2(vec[1], vec[0])
    DEC = np.arcsin(vec[2] / d)

    RA_deg = np.rad2deg(RA)
    if RA_deg < 0:
        RA_deg += 360

    return d, RA_deg, np.rad2deg(DEC)
    # return d, np.rad2deg(RA), np.rad2deg(DEC)

def sign(x):
    if x >= 0:
        return 1
    return -1

def gravAccel(mp, bodies):
    accel = np.array([0, 0, 0])
    
    for body in bodies:
        dist = np.linalg.norm(body.pos - mp.pos)
        grav_dir = (body.pos - mp.pos) / dist
        grav_mag = body.GM / dist**2
        
        accel = accel + grav_mag * grav_dir

    return accel

def readJSON(filename):
    try:
        with open(filename, "r") as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print("Error: File", filename, "not found!")
    except json.JSONDecodeError as e:
        print("Error: Failed to decode JSON in", filename, "\n", e)

def stepYoshida8(mp, bodies, dt):
    # - - - CONSTANTS - - -
    w1 = 0.311790812418427e0
    w2 = -0.155946803821447e1
    w3 = -0.167896928259640e1
    w4 = 0.166335809963315e1
    w5 = -0.106458714789183e1
    w6 = 0.136934946416871e1
    w7 = 0.629030650210433e0
    w0 = 1.65899088454396

    ds = [w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7]

    cs = [0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165]
    # - - -   - - -   - - -

    for i in range(15):
        mp.pos = mp.pos + mp.vel * cs[i] * dt
        accel = gravAccel(mp, bodies)
        mp.vel = mp.vel + accel * ds[i] * dt

    mp.pos = mp.pos + mp.vel * cs[15] * dt

def computeKepler(r, v, mu=1.3271244004193938e11):
    global km_per_AU
    
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    inclination = np.degrees(np.arccos(h[2] / h_mag))

    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    n_mag = np.linalg.norm(n)

    if n_mag != 0:
        omega = np.degrees(np.arccos(n[0] / n_mag))
        if n[1] < 0:
            omega = 360 - omega
    else:
        omega = 0

    e_vec = (1 / mu) * (np.cross(v, h) - mu * r / r_mag)
    eccentricity = np.linalg.norm(e_vec)

    if n_mag != 0:
        if eccentricity != 0:
            arg_periapsis = np.degrees(np.arccos(np.dot(n, e_vec) / (n_mag * eccentricity)))
            if e_vec[2] < 0:
                arg_periapsis = 360 - arg_periapsis
        else:
            arg_periapsis = 0
    else:
        arg_periapsis = 0

    if eccentricity != 0:
        true_anomaly = np.degrees(np.arccos(np.dot(e_vec, r) / (eccentricity * r_mag)))
        if np.dot(r, v) < 0:
            true_anomaly = 360 - true_anomaly
    else:
        true_anomaly = np.degrees(np.arccos(np.dot(r / r_mag, v / v_mag)))

    specific_energy = v_mag**2 / 2 - mu / r_mag
    if abs(eccentricity - 1) > 1e-8:
        semi_major_axis = -mu / (2 * specific_energy)
    else:
        semi_major_axis = np.inf

    if eccentricity < 1:
        E = 2 * np.arctan(np.tan(np.radians(true_anomaly) / 2) * np.sqrt((1 - eccentricity) / (1 + eccentricity)))
        if E < 0:
            E += 2 * np.pi
        mean_anomaly = np.degrees(E - eccentricity * np.sin(E))
        
    elif eccentricity > 1:
        F = 2 * np.arctanh(np.tan(np.radians(true_anomaly) / 2) * np.sqrt((eccentricity - 1) / (eccentricity + 1)))
        mean_anomaly = np.degrees(eccentricity * np.sinh(F) - F)
        
    else:
        mean_anomaly = None

    return {
        "a": semi_major_axis / AU,
        "e": eccentricity,
        "i": inclination,
        "lon_asc": omega,
        "arg_peri": arg_periapsis,
        "true_anomaly": true_anomaly,
        "mean anomaly": mean_anomaly
    }

def propagate(p0, v0, date_init, date_final, obs_code="Geocentric", mark_date=None, dt=None):
    global AU, day

    # generate bodies
    mp = MP("", p0, v0)
    
    bodies = []
    body_names = ["MERCURY BARYCENTER",
                  "VENUS BARYCENTER",
                  "EARTH BARYCENTER",
                  "MARS BARYCENTER",
                  "JUPITER BARYCENTER",
                  "SATURN BARYCENTER",
                  "URANUS BARYCENTER",
                  "NEPTUNE BARYCENTER",
                  "SUN"]
    
    body_GMs = [2.2031780000000021E+04,
                3.2485859200000006E+05,
                4.0350323550225981E+05,
                4.2828375214000022E+04,
                1.2671276480000021E+08,
                3.7940585200000003E+07,
                5.7945486000000080E+06,
                6.8365271005800236E+06,
                1.3271244004193938e11]

    for i in range(9):
        new_body = MainBody(body_names[i], np.array([0, 0, 0]), body_GMs[i])
        bodies.append(new_body)

    # set time parameters
    time_interval = (date_final - date_init).total_seconds()

    if mark_date:
        md_time_interval = (mark_date - date_init).total_seconds()

    if not dt:
        if mark_date:
            dt = min([md_time_interval / 10, 10 * day])
        else:
            dt = min([time_interval / 10, 10 * day])
    
    N_cycles = int(time_interval // dt) + 1
    date_final_actual = date_init + timedelta(seconds=N_cycles * dt)

    rhos = []
    RAs = []
    DECs = []
    mark_RA = None
    mark_DEC = None

    # numerically propagate orbit
    for cycle in range(N_cycles):
        cycle_date = date_init + timedelta(seconds=cycle * dt)
        cycle_date_str = cycle_date.strftime('%Y-%m-%dT%H:%M:%S')
        t = spice.str2et(cycle_date_str)
        
        for ib, body in enumerate(bodies):
            state, _ = spice.spkezr(body.name, t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
            body.pos = state[:3]

        stepYoshida8(mp, bodies, dt)

        if N_cycles > 5000 and cycle % 1000 == 0:
            percent_done = round(cycle / N_cycles * 100, 2)
            print(f"Propagating: {percent_done}%")

        if obs_code == "Geocentric":
            rho, RA, DEC = getRADEC(cycle_date, mp.pos)
            rhos.append(rho)
            RAs.append(RA)
            DECs.append(DEC)
        else:
            rho, RA, DEC = getRADECObsCode(cycle_date, obs_code, mp.pos)
            rhos.append(rho)
            RAs.append(RA)
            DECs.append(DEC)

        if mark_date and abs((mark_date - cycle_date).total_seconds()) < 2 * dt:
            mark_RA = RA
            mark_DEC = DEC

    return mp.pos, mp.vel, date_final_actual, rhos, RAs, DECs, mark_RA, mark_DEC

def getEarthPos(obs_date, coord='ecliptic'):
    t = spice.str2et(obs_date.strftime('%Y-%m-%dT%H:%M:%S'))
    if coord == 'ecliptic':
        earth_state, _ = spice.spkezr("EARTH", t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    else:
        earth_state, _ = spice.spkezr("EARTH", t, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    earth_pos = earth_state[:3]
    return earth_pos

def readObservatoryData(file_path):
    global observatories
    print("Reading observatory data...")
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("Code"):
                try:
                    code = line[0:3].strip()
                    longitude = float(line[3:13])
                    cos_phi = float(line[13:21])
                    sin_phi = float(line[21:30])
                    name = ''.join(line[30:len(line)-1])
                    observatories[code] = Observatory(code, longitude, cos_phi, sin_phi, name)
                except ValueError as e: # space telescopes etc. obviously don't have these
                    pass

    return observatories

def getObsCodePos(obs_date, obs_code, coord='ecliptic'):
    lat = observatories[obs_code].lat
    lon = observatories[obs_code].lon
    rho = observatories[obs_code].rho
    
    t = spice.str2et(obs_date.strftime('%Y-%m-%dT%H:%M:%S'))

    if coord == 'ecliptic':
        earth_pos, _ = spice.spkpos('EARTH', t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    else:
        earth_pos, _ = spice.spkpos('EARTH', t, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    
    R_EARTH = 6371.0088
    
    x_geo = rho * np.cos(lat_rad) * np.cos(lon_rad)
    y_geo = rho * np.cos(lat_rad) * np.sin(lon_rad)
    z_geo = rho * np.sin(lat_rad)
    geo_pos = np.array([x_geo, y_geo, z_geo])
    
    geo_to_eci = spice.pxform('ITRF93', 'J2000', t)
    eci_pos = geo_to_eci @ geo_pos

    if coord == 'ecliptic':
        eci_to_ecliptic = spice.pxform('J2000', 'ECLIPJ2000', t)
        ecliptic_pos = eci_to_ecliptic @ eci_pos
        obscode_pos = np.array(earth_pos) + ecliptic_pos
    else:
        obscode_pos = np.array(earth_pos) + eci_pos

    return obscode_pos

def getRADEC(obs_date, pos):
    earth_pos = getEarthPos(obs_date, 'equatorial')
    return cartezian2spherical(ecliptic2equatorial(pos) - earth_pos)

def getRADECObsCode(obs_date, obs_code, pos):
    obscode_pos = getObsCodePos(obs_date, obs_code, 'equatorial')
    return cartezian2spherical(ecliptic2equatorial(pos) - obscode_pos)

def equatorial2ecliptic(eq_pos):
    epsilon = np.deg2rad(23.439281)
    
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, np.cos(epsilon), np.sin(epsilon)],
        [0, -np.sin(epsilon), np.cos(epsilon)]
    ])
    
    ecliptic_coords = np.dot(rotation_matrix, eq_pos)
    
    return ecliptic_coords

def ecliptic2equatorial(ec_pos):
    epsilon = np.deg2rad(23.439281)
    
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, np.cos(epsilon), -np.sin(epsilon)],
        [0, np.sin(epsilon), np.cos(epsilon)]
    ])
    
    equatorial_coords = np.dot(rotation_matrix, ec_pos)
    
    return equatorial_coords

# === === === ORBIT DETERMINATION FUNCTIONS === === ===
def placeRelToEarth(obs_date, R, RA, DEC, coord='ecliptic'):
    t = spice.str2et(obs_date.strftime('%Y-%m-%dT%H:%M:%S'))
    if coord == 'ecliptic':
        earth_state, _ = spice.spkezr("EARTH", t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    else:
        earth_state, _ = spice.spkezr("EARTH", t, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    earth_pos = earth_state[:3]
    p = earth_pos + spherical2cartezian(R, np.deg2rad(RA), np.deg2rad(DEC))
    return p

def placeRelToObsCode(obs_date, obs_code, R, RA, DEC, coord='ecliptic'):
    obscode_pos = getObsCodePos(obs_date, obs_code, coord)
    p = obscode_pos + spherical2cartezian(R, np.deg2rad(RA), np.deg2rad(DEC))
    return p

def constructUnitVector(RA, DEC):
    x = np.cos(np.deg2rad(DEC)) * np.cos(np.deg2rad(RA))
    y = np.cos(np.deg2rad(DEC)) * np.sin(np.deg2rad(RA))
    z = np.sin(np.deg2rad(DEC))

    return np.array([x, y, z]) / np.linalg.norm(np.array([x, y, z]))

def get_rho2s(o1, o2, o3):
    mu = 1.3271244004193938e11

    rhat_1 = constructUnitVector(o1.RA, o1.DEC)
    rhat_2 = constructUnitVector(o2.RA, o2.DEC)
    rhat_3 = constructUnitVector(o3.RA, o3.DEC)

    R_1 = getEarthPos(o1.date, 'equatorial')
    R_2 = getEarthPos(o2.date, 'equatorial')
    R_3 = getEarthPos(o3.date, 'equatorial')

    t3mt1 = (o3.date - o1.date).total_seconds()
    rhatprime_2 = (rhat_3 - rhat_1) / t3mt1
    rhatprimeprime_2 = (rhat_3 - 2 * rhat_2 + rhat_1) / (t3mt1**2)

    tau_1 = (o1.date - o2.date).total_seconds()
    tau_3 = (o3.date - o2.date).total_seconds()
    tau = (o3.date - o1.date).total_seconds()

    p_1 = np.cross(rhat_2, rhat_3)
    p_2 = np.cross(rhat_1, rhat_3)
    p_3 = np.cross(rhat_1, rhat_2)

    D_0 = np.dot(rhat_1, p_1)
    D_11 = np.dot(R_1, p_1)
    D_12 = np.dot(R_1, p_2)
    D_13 = np.dot(R_1, p_3)
    D_21 = np.dot(R_2, p_1)
    D_22 = np.dot(R_2, p_2)
    D_23 = np.dot(R_2, p_3)
    D_31 = np.dot(R_3, p_1)
    D_32 = np.dot(R_3, p_2)
    D_33 = np.dot(R_3, p_3)

    A = 1/D_0 * (-D_12 * tau_3 / tau + D_22 + D_32 * tau_1 / tau)
    B = 1 / (6 * D_0) * (D_12 * (tau_3**2 - tau**2) * tau_3 / tau + D_32 * (tau**2 - tau_1**2) * tau_1 / tau)
    E = np.dot(R_2, rhat_2)

    R_2sq = np.dot(R_2, R_2)
    
    a = -(A**2 + 2*A*E + R_2sq)
    b = -2*mu*B*(A+E)
    c = -mu**2 * B**2

    rho_2s = np.roots([1, 0, a, 0, 0, b, 0, 0, c]).real

    return rho_2s

def perfectV0(R1, deltaR, o1, o2):
    R2 = R1 + deltaR

    # p0 = equatorial2ecliptic(placeRelToEarth(o1.date, R1, o1.RA, o1.DEC, 'equatorial'))
    # pf = equatorial2ecliptic(placeRelToEarth(o2.date, R2, o2.RA, o2.DEC, 'equatorial'))
    
    p0 = equatorial2ecliptic(placeRelToObsCode(o1.date, o1.obs_code, R1, o1.RA, o1.DEC, 'equatorial'))
    pf = equatorial2ecliptic(placeRelToObsCode(o2.date, o2.obs_code, R2, o2.RA, o2.DEC, 'equatorial'))

    date_final = o2.date
    date_init = o1.date

    delta_time = (date_final - date_init).total_seconds()
    v0 = (pf - p0) / delta_time # km s-1

    error = float('Inf')
    tol = 1e-4
    while error > tol:
        p_final, v_final, date_final_actual, _, _, _, _, _ = propagate(p0, v0, date_init, date_final)
        R_final, RA_final, DEC_final = getRADECObsCode(o2.date, o2.obs_code, pf)

        p_err = pf - p_final
        v0 = v0 + p_err / delta_time

        error = np.linalg.norm(p_err)

    # print(f"Perfected v0 with {error} km of error.")
    return v0

# Quasi-Herget method
def determineOrbit(obs_all):
    mu = 1.3271244004193938e11
    
    o1 = obs_all[0]
    o3 = obs_all[len(obs_all) - 1]
    o2 = obs_all[int(len(obs_all) / 2)]

    date_init = o1.date
    date_final = o2.date
    date_check = o3.date

    rho_2s = abs(get_rho2s(o1, o2, o3))

    # --- initial observation guess ---
    # initially assuming perfect observation with no errors
    R1 = max(rho_2s) # max. is usually the closest one
    deltaR = 0 * AU
    R2 = R1 + deltaR

    # p0 = equatorial2ecliptic(placeRelToEarth(o1.date, R1, o1.RA, o1.DEC, 'equatorial'))
    # pf = equatorial2ecliptic(placeRelToEarth(o2.date, R2, o2.RA, o2.DEC, 'equatorial'))
    p0 = equatorial2ecliptic(placeRelToObsCode(o1.date, o1.obs_code, R1, o1.RA, o1.DEC, 'equatorial'))
    pf = equatorial2ecliptic(placeRelToObsCode(o2.date, o2.obs_code, R2, o2.RA, o2.DEC, 'equatorial'))

    delta_time = (date_final - date_init).total_seconds()
    # v0 = perfectV0(R1, deltaR, o1, o2)
    v0 = (mu / np.linalg.norm(p0))**0.5 * np.array([-p0[1] / np.linalg.norm(p0), p0[0] / np.linalg.norm(p0), 0])

    orbital_elems = computeKepler(p0, v0)
    print("Initial guess:")
    print(orbital_elems)
    # --- --- --- --- --- --- --- --- ---

    adjust_factor = 1
    adjust_pfactor = 1

    good_fit = False
    retry_count = 0
    max_retry = 50
    while (not good_fit) and retry_count <= max_retry:
        err_val = 0
        for idx_o, o in enumerate(obs_all):
            if o.date != date_init:
                p_check, v_check, date_check_actual, rhos, RAs, DECs, _, _ = propagate(p0, v0, date_init, o.date)
                d_prop, RA_prop, DEC_prop = getRADECObsCode(o.date, o.obs_code, p_check)

                RA_err = o.RA - RA_prop
                DEC_err = o.DEC - DEC_prop

                err_val += RA_err**2 + DEC_err**2

        print(f"Iter: {retry_count}, errScore: {err_val}, errRA: {RA_err}, errDEC: {DEC_err}")
        orbital_elems = computeKepler(p0, v0)
        # print(orbital_elems)

        if abs(RA_err) < 1 * arcsecond and abs(DEC_err) < 1 * arcsecond:
            good_fit = True
        else:
            adjust_vals = [0, 0, 0, 0, 0, 0]
            adjust_vecs = [np.array([0.01, 0, 0]),
                           np.array([-0.01, 0, 0]),
                           np.array([0, 0.01, 0]),
                           np.array([0, -0.01, 0]),
                           np.array([0, 0, 0.01]),
                           np.array([0, 0, -0.01])]

            for i in range(6):
                adjust_vecs[i] *= adjust_factor

            # adjust vel
            for i in range(6):
                v0_1 = v0 + adjust_vecs[i]

                adjust_vals[i] = 0
                for idx_o, o in enumerate(obs_all):
                    if o.date != date_init:
                        p_check_1, v_check_1, _, _, _, _, _, _ = propagate(p0, v0_1, date_init, o.date)
                        _, RA_prop_1, DEC_prop_1 = getRADECObsCode(o.date, o.obs_code, p_check_1)

                        RA_err_1 = o.RA - RA_prop_1
                        DEC_err_1 = o.DEC - DEC_prop_1

                        adjust_vals[i] += RA_err_1**2 + DEC_err_1**2

            idx_min = np.argmin(adjust_vals)
            if adjust_vals[idx_min] < err_val:
                v0 = v0 + adjust_vecs[idx_min]
                adjust_factor *= 2
            else:
                adjust_factor *= 0.1

            # adjust pos
            adjust_pvals = [0, 0, 0, 0, 0, 0]
            adjust_pvecs = [np.array([0.01, 0, 0]),
                           np.array([-0.01, 0, 0]),
                           np.array([0, 0.01, 0]),
                           np.array([0, -0.01, 0]),
                           np.array([0, 0, 0.01]),
                           np.array([0, 0, -0.01])]

            for i in range(6):
                adjust_pvecs[i] *= adjust_pfactor
                
            for i in range(6):
                p0_1 = p0 + adjust_pvecs[i]
                adjust_pvals[i] = 0
                
                for idx_o, o in enumerate(obs_all):
                    if o.date != date_init:
                        p_check_1, v_check_1, _, _, _, _, _, _ = propagate(p0_1, v0, date_init, o.date)
                        _, RA_prop_1, DEC_prop_1 = getRADECObsCode(o.date, o.obs_code, p_check_1)

                        RA_err_1 = o.RA - RA_prop_1
                        DEC_err_1 = o.DEC - DEC_prop_1

                        adjust_pvals[i] += RA_err_1**2 + DEC_err_1**2

            idx_min = np.argmin(adjust_pvals)
            if adjust_pvals[idx_min] < err_val:
                p0 = p0 + adjust_pvecs[idx_min]
                adjust_pfactor *= 2
            else:
                adjust_pfactor *= 0.1

        retry_count += 1

    pf, vf, date_check_actual, rhos, RAs, DECs, _, _ = propagate(p0, v0, date_init, obs_all[-1].date)
    orbital_elems = computeKepler(p0, v0)
    print("Final fit:")
    print(p0, v0)
    print(orbital_elems)

    return p0, v0, date_init, orbital_elems

# === === === === === === === === === === === === === === === ===

def readObsFile(filename='primary.obs'):
    obses = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            new_o = Obs(line)
            obses.append(new_o)

    return obses

class AstrometryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("AOKOI " + str(version))

        current_row = 0
        
        icon_path = "AOKOI2.png"
        if os.path.exists(icon_path):
            self.ico_img = tk.PhotoImage(file=icon_path).subsample(6, 6)
            self.root.iconphoto(False, self.ico_img)
            header_label = tk.Label(root, image=self.ico_img, bg="#212331")
            header_label.grid(row=current_row, column=0, columnspan=2, sticky="n", pady=10)
            current_row += 1
            
        self.root.configure(bg="#212331")

        # Input fields
        tk.Label(root, text="Observations File:", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.primary_file_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.primary_file_entry.insert(0, "primary.obs")
        self.primary_file_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        s = ttk.Style() # silly ttk needs styles to style widgets
        s.configure('AOKOI.TRadiobutton',
                    background='#212331',
                    foreground='#eeeeee')
        
        tk.Label(root, text="Orbit Determination:", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.OD_var = tk.StringVar()
        self.OD_var.set('quasi-herget')
        self.OD_quasiherget = ttk.Radiobutton(
            root,
            text="Quasi-Herget",
            value='quasi-herget',
            variable=self.OD_var,
            style='AOKOI.TRadiobutton')
        self.OD_quasiherget.grid(row=current_row, column=1, padx=10, pady=5, sticky="w")
        current_row += 1

        self.run_button = tk.Button(root, text="RUN", command=self.runAnalysis, bg="#0078d7", fg="#fff", width=20)
        self.run_button.grid(row=current_row, column=0, columnspan=2, pady=20)
        current_row += 1

        self.about_button = tk.Button(root, text="About / Help / License", command=self.showAbout, bg="#aa1100", fg="#eee", width=20)
        self.about_button.grid(row=current_row, column=0, columnspan=2, pady=10)
        current_row += 1

        print(title)
        print("AOKOI", version, "initialized.")

    def showAbout(self):
        about_window = tk.Toplevel(self.root)
        about_window.title("About")
        about_window.geometry("900x650")
        about_window.resizable(False, False)

        tk.Label(about_window, text="About AOKOI", font=("Arial", 14, "bold")).pack(pady=10)
        tk.Label(
            about_window,
            text=(help_text + "\n" + authors),
            wraplength=800,
            justify="left",
        ).pack(pady=10)

        tk.Label(about_window, text="License", font=("Arial", 14, "bold")).pack(pady=10)
        tk.Label(
            about_window,
            text=(license_text),
            wraplength=800,
            justify="left",
        ).pack(pady=10)

        tk.Button(about_window, text="Close", command=about_window.destroy, bg="#0078d7", fg="#fff", width=10).pack(pady=10)

    def runAnalysis(self):
        print("AOKOI is running...")
        
        # Get inputs
        primary_file = self.primary_file_entry.get()
        OD = self.OD_var.get()

        #try:
        # load SPICE kernels
        print("Loading SPICE kernels...")
        spice.furnsh('data/naif0012.tls')
        spice.furnsh('data/de440.bsp')
        spice.furnsh('data/pck00011.tpc')
        spice.furnsh('data/earth_000101_250316_241218.bpc')

        # load observatory data
        readObservatoryData('data/obscodes.txt')
        
        # -- orbit determination --
        # read observation files
        print("Reading input observations from file...")
        obs_all = readObsFile(primary_file)
        
        if OD == 'quasi-herget':
            print("Orbit determination via quasi-Herget method...")
            p0, v0, FO_epoch, orbital_elements = determineOrbit(obs_all)
        # -- -- -- -- -- -- -- -- --

        # check if all observations have same obscode
        same_obs_code = True
        for obs in obs_all[1:]:
            if obs.obs_code != obs_all[0].obs_code:
                same_obs_code = False
                break

        print("Plotting trajectory on celestial sphere...")
        if not same_obs_code:
            _, _, _, _, RAs, DECs, _, _ = propagate(p0, v0, obs_all[0].date, obs_all[-1].date, obs_code="Geocentric", mark_date=None, dt=None)
        else:
            _, _, _, _, RAs, DECs, _, _ = propagate(p0, v0, obs_all[0].date, obs_all[-1].date, obs_code=obs_all[0].obs_code, mark_date=None, dt=None)
            
        self.plotCelestialSphere(obs_all, RAs, DECs)

        print("Plotting 3D orbit map...")
        self.plotOrbit3D(orbital_elements)
        
        print("Done!")

        #except Exception as e:
            #messagebox.showerror("Error", f"An error occurred: {e}")

    def plotCelestialSphere(self, obs_all, RAs, DECs):
        plt.plot(RAs, DECs, label="Estimated Orbit", linestyle="--")
        sc = plt.scatter([o.RA for o in obs_all], [o.DEC for o in obs_all], label="Known Obs.", color="blue")
        plt.grid()
        plt.title("Astrometry Chart")
        plt.legend(loc="upper left")
        plt.xlabel("RA (deg)")
        plt.ylabel("DEC (deg)")
        plt.axis('equal')
        mplcursors.cursor(sc, hover=False).connect("add", lambda sel: sel.annotation.set_text(str(obs_all[sel.index])))
        # plt.show()

    def plotOrbit3D(self, orbital_elements, poly=1000):
        eccentricity = orbital_elements["e"]
        argument_perihelion = orbital_elements["arg_peri"]
        semimajor_axis = orbital_elements["a"]
        ascending_node = orbital_elements["lon_asc"]
        inclination = orbital_elements["i"]
        
        planets = []
        # print("Reading JSON major planet data...")
        planet_data = readJSON("data/planets.json")
        for p in planet_data:
            name = p["name"]
            a = float(p["a"])
            e = float(p["e"])
            i = float(p["i"])
            OMG = float(p["Omega"])
            o = float(p["-omega"])
            L = float(p["L"])

            new_planet = PlotPlanet(name, a, e, i, OMG, o, L)
            planets.append(new_planet)

        planet_colors =\
                      ["gray",
                       "sandybrown",
                       "royalblue",
                       "firebrick",
                       "peru",
                       "moccasin",
                       "powderblue",
                       "slateblue"]
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        argument_perihelion = np.deg2rad(argument_perihelion)
        ascending_node = np.deg2rad(ascending_node)
        inclination = np.deg2rad(inclination)
        
        true_anomaly = np.linspace(0, 2 * np.pi, poly)

        r = (semimajor_axis * (1 - eccentricity**2)) / (1 + eccentricity * np.cos(true_anomaly))

        x = r * np.cos(true_anomaly)
        y = r * np.sin(true_anomaly)
        z = np.zeros_like(true_anomaly)

        x_rot = x * (np.cos(argument_perihelion) * np.cos(ascending_node) - np.sin(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination)) - \
                y * (np.sin(argument_perihelion) * np.cos(ascending_node) + np.cos(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination))

        y_rot = x * (np.cos(argument_perihelion) * np.sin(ascending_node) + np.sin(argument_perihelion) * np.cos(ascending_node) * np.cos(inclination)) + \
                y * (np.cos(argument_perihelion) * np.cos(ascending_node) - np.sin(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination))

        z_rot = x * np.sin(argument_perihelion) * np.sin(inclination) + y * np.sin(argument_perihelion) * np.sin(inclination)

        color = "k"

        ax.plot(x_rot, y_rot, z_rot, color=color)

        for p_idx in range(len(planets)):
            p = planets[p_idx]
            
            eccentricity = p.e
            argument_perihelion = np.radians(p.o)
            semimajor_axis = p.a
            ascending_node = np.radians(p.OMG)
            inclination = np.radians(p.i)
            
            true_anomaly = np.linspace(0, 2 * np.pi, 1000)

            r = (semimajor_axis * (1 - eccentricity**2)) / (1 + eccentricity * np.cos(true_anomaly))

            x = r * np.cos(true_anomaly)
            y = r * np.sin(true_anomaly)
            z = np.zeros_like(true_anomaly)

            x_rot = x * (np.cos(argument_perihelion) * np.cos(ascending_node) - np.sin(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination)) - \
                    y * (np.sin(argument_perihelion) * np.cos(ascending_node) + np.cos(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination))

            y_rot = x * (np.cos(argument_perihelion) * np.sin(ascending_node) + np.sin(argument_perihelion) * np.cos(ascending_node) * np.cos(inclination)) + \
                    y * (np.cos(argument_perihelion) * np.cos(ascending_node) - np.sin(argument_perihelion) * np.sin(ascending_node) * np.cos(inclination))

            z_rot = x * np.sin(argument_perihelion) * np.sin(inclination) + y * np.sin(argument_perihelion) * np.sin(inclination)

            color = planet_colors[p_idx]
            ax.plot(x_rot, y_rot, z_rot, color=color, lw=3, linestyle="--")
        
        ax.scatter(0, 0, 0, color='yellow', label='Sol Barycenter')
        
        ax.set_title('Orbit Plot')
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
        
        # setting map limits
        allmin = abs(min([min(x_rot), min(y_rot), min(z_rot)])) * 1.5
        allmax = abs(max([max(x_rot), max(y_rot), max(z_rot)])) * 1.5
        maplim = max([allmin, allmax])
        
        ax.set_xlim(-maplim, maplim)
        ax.set_ylim(-maplim, maplim)
        ax.set_zlim(-maplim, maplim)
        plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = AstrometryApp(root)
    root.mainloop()
