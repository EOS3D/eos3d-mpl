import math
import numpy as np
import pandas as pd
import ast
import warnings

# Check for pycurl installation
try:
    import pycurl
except:
    warnings.warn("Pycurl not installed: getTLE() not functional", ImportWarning)
    pycurl_present = False
else:
    pycurl_present = True

try:
    from sgp4.earth_gravity import wgs84
    from sgp4.io import twoline2rv
except:
    warnings.warn("SPG4 not installed: getTLE() and SGP4 propagation not functional", ImportWarning)
    sgp4_present = False
else:
    sgp4_present = True

from io import BytesIO
import os

#import cartopy.crs as ccrs
#import cartopy.feature as cfeature



import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib import style
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
from matplotlib import colors as mcolors

from scipy import optimize

import datetime
from dateutil import tz
import time
import calendar

from eos3dmpl.eos_aedata import *
import eos3dmpl.eos_conv as conv
import eos3dmpl.eos_config

class Eos():
    res = 250
    ss = False
    # Available colors
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                    for name, color in colors.items())
    l_colors = [name for hsv, name in by_hsv]
    dCRD = {
        'dX': {
              'Sun'     : 0,
              'Mercury' : 0,
              'Venus'   : 0,
              'Earth'   : 0,
              'Moon'    : 0,
              'Mars'    : 0,
              'Jupiter' : 0,
              'Saturn'  : 0,
              'Uranus'  : 0,
              'Neptune' : 0,
              'Pluto'   : 0
        },
        'dY': {
              'Sun'     : 0,
              'Mercury' : 0,
              'Venus'   : 0,
              'Earth'   : 0,
              'Moon'    : 0,
              'Mars'    : 0,
              'Jupiter' : 0,
              'Saturn'  : 0,
              'Uranus'  : 0,
              'Neptune' : 0,
              'Pluto'   : 0
        },
        'dZ': {
              'Sun'     : 0,
              'Mercury' : 0,
              'Venus'   : 0,
              'Earth'   : 0,
              'Moon'    : 0,
              'Mars'    : 0,
              'Jupiter' : 0,
              'Saturn'  : 0,
              'Uranus'  : 0,
              'Neptune' : 0,
              'Pluto'   : 0
        }
    }
    
    def degrees(self, pref=True):
        try:
            bool(pref)
        except:
            raise ValueError('Invalid degrees bool')
        else:
            if pref == True:
                try:
                    self.ta = np.deg2rad(self.ta)
                except:
                    pass
                else:
                    pass
                self.i = np.deg2rad(self.i)
                self.om = np.deg2rad(self.om)
                self.w = np.deg2rad(self.w)
            else:
                pass
            
    def units(self, unit='km'):
        self.a = conv.unit(self.a, 'm', unit)
        self.h_a = conv.unit(self.h_a, 'm', unit)
        self.Q = conv.unit(self.Q, 'm', unit)
        self.h_Q = conv.unit(self.h_Q, 'm', unit)
        self.q = conv.unit(self.q, 'm', unit)
        self.h_q = conv.unit(self.h_q, 'm', unit)
        self.R = conv.unit(self.R, 'm', unit)
        if self.r != None:
            self.r = conv.unit(self.r, 'm', unit)
        self.unit = unit
        
    def tp_solve(self,E):
        #if E < 0:
        #    return 99999
        #else:
        return np.sqrt((self.a**3)/(G*self.M_b))*(E - self.e*np.sin(E)) - self.tp
    
    def M_solve(self,E,*data):
        M, e = data
        return E - e*np.sin(E) - M
    
    def ta2t(self,ta,a):
        n = np.sqrt((G*dict_M[self.body])/(a**3))
        T = 2*np.pi*(1./n)
        E = np.arccos((self.e + np.cos(ta))/(1 + self.e*np.cos(ta)))
        #E0 = np.arccos((self.e + 1)/(1 + self.e))
        M = (E - self.e*np.sin(E))
        #M0 = (E0 - self.e*np.sin(E0))
        #dt = (M - M0)/n
        dt = M/n
        rest = self.res
        sc = np.pi
        mlt = 1
        mlt_l = []
        for idx,val in enumerate(dt):
            if mlt % 2 != 0:
                dt[idx] += (mlt-((mlt+1)/2))*T
            elif mlt % 2 == 0:
                dt[idx] = (mlt/2)*T - dt[idx]
            if ta[idx] > sc:
                mlt += 1
                sc += np.pi
            mlt_l += [mlt/2.]
        return {'dt':dt,'mlt':mlt_l}
    
    def orbit(self, a=None, h_a=None, e=None, ta=None, tp=None, i=0, om=0, w=0, Q=None, h_Q=None, q=None, h_q=None, unit_i='km', unit_o='km', deg=True, body='Earth'):
        self.a = a
        self.h_a = h_a
        self.e = e
        self.ta = ta
        self.tp = tp
        self.i = i
        self.om = om
        self.w = w
        self.Q = Q
        self.h_Q = h_Q
        self.q = q
        self.h_q = h_q
        self.unit = unit_i
        self.body = body
        self.r = None
        self.M = None
        self.E = None
        self.n = None
        self.V = None
        self.V_esc = None
        self.curr = False
        
        # Input value integrity check
        # Celestial body
        try:
            self.M_b = dict_M[body]
            self.R = dict_R[body]
            self.mu = self.M_b * G
        except:
            raise ValueError('Invalid celestial body')
        # SMA
        if self.a == None and self.h_a != None:
            try:
                self.h_a = float(self.h_a)
            except:
                raise ValueError('Invalid SMA altitude')
            else:
                self.a = self.R + unit(self.h_a, self.unit, 'm')
        elif self.a != None and self.h_a != None:
            raise SyntaxError('Conflicting SMA radius and altitude')
        elif self.a != None:
            try:
                self.a = float(self.a)
            except:
                raise ValueError('Invalid SMA')
            else:
                self.a = conv.unit(self.a, self.unit, 'm')
                self.h_a = self.a - self.R
        # Apoapsis
        if self.Q == None and self.h_Q != None:
            try:
                self.h_Q = float(self.h_Q)
            except:
                raise ValueError('Invalid apoapsis altitude')
            else:
                self.Q = self.R + unit(self.h_Q, self.unit, 'm')
        elif self.Q != None and self.h_Q != None:
            raise SyntaxError('Conflicting apoapsis radius and altitude')
        elif self.Q != None:
            try:
                self.Q = float(self.Q)
            except:
                raise ValueError('Invalid apoapsis')
            else:
                self.Q = conv.unit(self.Q, self.unit, 'm')
                self.h_Q = self.Q - self.R
        # Periapsis
        if self.q == None and self.h_q != None:
            try:
                self.h_q = float(self.h_q)
            except:
                raise ValueError('Invalid periapsis altitude')
            else:
                self.q = self.R + unit(self.h_q, self.unit, 'm')
                self.h_q = self.q - self.R
        elif self.q != None and self.h_q != None:
            raise SyntaxError('Conflicting periapsis radius and altitude')
        elif self.q != None:
            try:
                self.q = float(self.q)
            except:
                raise ValueError('Invalid periapsis')
            else:
                self.q = conv.unit(self.q, self.unit, 'm')
        # Eccentricity
        if self.e == None and self.Q >= self.q:
            try:
                self.e = float(self.e)
            except:
                try:
                    self.Q = float(self.Q)
                except:
                    raise ValueError('Invalid apoapsis radius')
                else:
                    try:
                        self.q = float(self.q)
                    except:
                        raise ValueError('Invalid periapsis radius')
                    else:
                        self.e = (self.Q - self.q)/(self.Q + self.q)
        elif self.e == None and self.Q < self.q:
            raise SyntaxError('Apoapsis must be smaller than periapsis')
        elif self.e != None and self.Q != None and self.q != None:
            raise SyntaxError('Conflicting eccentricity and apsides')
        else:
            pass
        # SMA (Semi-Major Axis) and Apsides conflict check  
        if self.a != None and (self.Q != None or self.q != None):
            raise SytaxError('Conflicting SMA and apsides')
        elif self.a == None and self.Q != None and self.q != None:
            self.a = (self.Q + self.q)/2
        # SMA and eccentricity to Apsides
        if (self.a != None and self.e != None) and (self.Q == None and self.q == None):
            self.Q = self.a*(1+self.e)
            self.h_a = self.a - self.R
            self.h_Q = self.Q - self.R
            self.q = self.a*(1-self.e)
            self.h_q = self.q - self.R
        else:
            self.h_a = self.a - self.R
            self.h_Q = self.Q - self.R
            self.h_q = self.q - self.R
        # Period calculation
        self.T = 2*np.pi*np.sqrt((self.a**3)/(self.mu))
        # Set True Anomly to 0 if not given
        if self.ta == None and self.tp == None:
            self.ta = 0
            self.curr = False
        else:
            self.curr = True
        # TP (Time from Periapsis) and True Anomaly conflict check
        if self.tp != None and self.ta != None:
            raise SytaxError('Conflicting TP and true anomaly')
        else:
            if self.ta != None:
                try:
                    self.ta = float(self.ta)
                except:
                    raise ValueError('Invalid true anomly')
                else:
                    self.degrees(deg)
                    self.E = 2*np.arctan(np.sqrt((1-self.e)/(1+self.e))*np.tan(self.ta/2))
            if tp != None:
                try:
                    self.tp = float(self.tp)
                except:
                    raise ValueError('Invalid TP')
                else:
                    self.degrees(deg)
                    E = optimize.fsolve(self.tp_solve, 1)
                    self.E = E[0]
                    self.ta = 2*np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(self.E/2))
                    if self.ta < 0:
                        self.ta += 2*np.pi
            if self.E != None:
                self.M = self.E - self.e*np.sin(self.E)
                self.n = np.sqrt(self.mu/(self.a**3))
                self.tfp = self.M/self.n
                if self.tfp < 0:
                    self.tfp += self.T
                elif self.tfp >= self.T:
                    self.tfp += -self.T
                self.ttp = self.T - self.tfp
                self.r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(self.ta)))
                self.V = np.sqrt(self.mu*((2./self.r)-(1./self.a)))/1000.
                self.V_esc = np.sqrt(2*self.mu/self.r)/1000.
        self.units(unit_o)
    
    def init(self, res=250, center='Sun'):
        self.res = res
        if center in l_CB:
            pass
        else:
            raise ValueError('Invalid plot center')
        self.center = center
        
    def view_init(self, res=250, center='Sun'):
        style.use('bmh')
        mpl.rcParams['legend.fontsize'] = 8
        self.fig = plt.figure()
        self.ax = plt.subplot(111, projection='3d')
        self.ax.set_aspect('equal')
        self.ax.set_adjustable('box')
        self.ax.set_anchor('C')
        self.ax.axis('equal')
        plt.rc('font', family='serif')
        self.res = res
        if center in l_CB:
            pass
        else:
            raise ValueError('Invalid plot center')
        self.center = center
    
    def kepler(self, rev=1, data=False, vis='v-line', marker=True, curr_background=True, curr_marker=True, curr_symbol='◉', curr_color='#D62829', curr_label='_nolegend_', curr_fontdict=None, linewidth=1.5, sgp=False, l1=None, l2=None, time='now', designation='NA', res=250):
        if res != None and type(res) == int:
            self.res = res
        try:
            designation = str(designation)
        except:
            raise ValueError('Invalid designation')
        try:
            marker = bool(marker)
        except:
            raise ValueError('Invalid marker bool')
        else:
            try:
                curr_marker = bool(curr_marker)
            except:
                raise ValueError('Invalid curr_marker bool')
            else:
                pass
        self.designation = designation
        self.marker = marker
        self.curr = curr_marker
        self.rev = rev
        self.vis = vis
        if 2*self.rev*self.res < 1:
            spacing = 1
            warnings.warn("Resolution too low; set to 1", UserWarning)
        else:
            spacing = np.round(2*self.rev*self.res)
        l_ta = np.linspace(0, 2*self.rev*np.pi, np.round(spacing))
        M = dict_M[self.body]
        ta1 = np.array([])
        r1 = np.array([])
        c1 = np.array([])
        x1 = np.array([])
        y1 = np.array([])
        z1 = np.array([])
        V1 = np.array([])
        a = conv.unit(self.a, self.unit, 'm')
        V_a = np.sqrt(G*M*((2./(a*(1+self.e)))-(1./a)))
        V_p = np.sqrt(G*M*((2./(a*(1-self.e)))-(1./a)))
        w = self.w + 0.5*np.pi
        om = self.om
        i = self.i
        for ta in l_ta:
            r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(ta)))
            r_temp = conv.unit(r, self.unit, 'm')
            V = np.sqrt(G*M*((2./r_temp)-(1./a)))
            if V_a == V_p:
                c_temp = 127
            else:
                c_temp = int(np.round(((V - V_a)/(V_p - V_a))*255.))
            x_temp = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
            x_temp += self.dCRD['dX'][self.body]
            y_temp = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
            y_temp += self.dCRD['dY'][self.body]
            z_temp = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
            z_temp += self.dCRD['dZ'][self.body]
            ta1 = np.append(ta1, ta)
            r1 = np.append(r1, r)
            c1 = np.append(c1, c_temp)
            x1 = np.append(x1, x_temp)
            y1 = np.append(y1, y_temp)
            z1 = np.append(z1, z_temp)
            V1 = np.append(V1, V/1000.)
        # Current location marker coordinates
        if self.ta != None:
            r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(self.ta)))
            ta = self.ta
            r_curr = r
            x_curr = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
            x_curr += self.dCRD['dX'][self.body]
            y_curr = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
            y_curr += self.dCRD['dY'][self.body]
            z_curr = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
            z_curr += self.dCRD['dZ'][self.body]
        # SGP4 TLE
        if sgp == True:
            if not sgp4_present:
                return warnings.warn("SPG4 not initialized: SGP4 propagation not functional", UserWarning)
            satellite = twoline2rv(l1, l2, wgs84)
            if time == 'now':
                currtime = pd.Timestamp.today(tz='Europe/Amsterdam').tz_convert('UTC')
            else:
                currtime = time
            currpos, currvel = satellite.propagate(currtime.year, currtime.month, currtime.day, currtime.hour, currtime.minute, currtime.second)
            x_curr = currpos[0]
            y_curr = currpos[1]
            z_curr = currpos[2]
        # Apoapsis coordinates
        ta = np.pi
        r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(ta)))
        w = self.w + 0.5*np.pi
        om = self.om
        i = self.i
        x_ap = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
        x_ap += self.dCRD['dX'][self.body]
        y_ap = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
        y_ap += self.dCRD['dY'][self.body]
        z_ap = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
        z_ap += self.dCRD['dZ'][self.body]
        # Periapsis coordinates
        ta = 0
        r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(ta)))
        x_pe = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
        x_pe += self.dCRD['dX'][self.body]
        y_pe = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
        y_pe += self.dCRD['dY'][self.body]
        z_pe = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
        z_pe += self.dCRD['dZ'][self.body]
        # Ascending node coodinates
        ta = -self.w
        r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(ta)))
        x_an = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
        x_an += self.dCRD['dX'][self.body]
        y_an = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
        y_an += self.dCRD['dY'][self.body]
        z_an = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
        z_an += self.dCRD['dZ'][self.body]
        # Descending node coodinates
        ta = -self.w + np.pi
        r = ((self.a*(1 - self.e**2)) / (1 + self.e * np.cos(ta)))
        x_dn = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.cos(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.cos(om) - r*np.cos(ta)*np.sin(w)*np.sin(om) - r*np.sin(ta)*np.cos(w)*np.sin(om)
        x_dn += self.dCRD['dX'][self.body]
        y_dn = r*np.cos(ta)*np.cos(w)*np.cos(i)*np.sin(om) - r*np.sin(ta)*np.sin(w)*np.cos(i)*np.sin(om) + r*np.cos(ta)*np.sin(w)*np.cos(om) + r*np.sin(ta)*np.cos(w)*np.cos(om)
        y_dn += self.dCRD['dY'][self.body]
        z_dn = -r*np.cos(ta)*np.cos(w)*np.sin(i) + r*np.sin(ta)*np.sin(w)*np.sin(i)
        z_dn += self.dCRD['dZ'][self.body]
        
        if vis != None:
            c1 = c1.astype(int)
            # Label integrity check
            try:
                curr_label = str(curr_label)
            except:
                raise ValueError('Invalid label')
            else:
                pass
            # Visualization
            try:
                linewidth = float(linewidth)
            except:
                raise ValueError('Invalid linewidth')
            if type(vis) == tuple:
                if self.curr==True:
                    self.ax.plot(x1, y1, z1, c=vis, label='_nolegend_', linewidth=linewidth, zorder=0)
                else:
                    self.ax.plot(x1, y1, z1, c=vis, label=curr_label, linewidth=linewidth, zorder=0)
            elif vis=='v-line':
                m = cm.ScalarMappable(cmap=cm.jet)
                m.set_array(V1)
                if self.curr==True:
                    self.ax.scatter(x1, y1, z1, c=plt.cm.jet(c1), label='_nolegend_', s=linewidth**2, zorder=0)
                else:
                    self.ax.scatter(x1, y1, z1, c=plt.cm.jet(c1), label=curr_label, s=linewidth**2, zorder=0)
                #cbaxes = self.fig.add_axes([0.9, 0.15, 0.0245, 0.7]) 
                cbaxes = self.fig.add_axes([0.884, 0.3, 0.0153, 0.4]) 
                #plt.colorbar(m, cax=cbaxes, label='Velocity [km/s]',use_gridspec=True)
                self.fig.colorbar(m, cax=cbaxes,label='Velocity [km/s]',use_gridspec=True)
            elif vis=='gradient':
                c_lin = np.around(np.linspace(0,255,(len(x1)))).astype(int)
                if self.curr==True:
                    self.ax.scatter(x1, y1, z1, c=plt.cm.Oranges(c_lin), label='_nolegend_', s=linewidth**2, zorder=0)
                else:
                    self.ax.scatter(x1, y1, z1, c=plt.cm.Oranges(c_lin), label=curr_label, s=linewidth**2, zorder=0)
            elif vis=='seq' or vis=='sequential':
                if self.curr==True:
                    self.ax.plot(x1, y1, z1, label='_nolegend_', linewidth=linewidth, zorder=0)
                else:
                    self.ax.plot(x1, y1, z1, label=curr_label, linewidth=linewidth, zorder=0)
            else:
                if vis in self.l_colors or '#' in vis:
                    if self.curr==True:
                        self.ax.plot(x1, y1, z1, c=vis, label='_nolegend_', linewidth=linewidth, zorder=0)
                    else:
                        self.ax.plot(x1, y1, z1, c=vis, label=curr_label, linewidth=linewidth, zorder=0)
                else:
                    raise ValueError('Invalid visual mode')
            # Marker bool
            try:
                self.marker = bool(self.marker)
            except:
                raise ValueError('Invalid marker bool')
            else:
                pass
            if self.marker == True:
                # Apoapsis marker
                self.ax.text(x_ap, y_ap, z_ap, 'Ap', horizontalalignment='center', verticalalignment='center', backgroundcolor='white', zorder=249,
                             fontdict={'family': 'serif', 'color':  'black', 'weight': 'bold', 'size': 8}).set_bbox(dict(facecolor='#B2E1FF', alpha=0.8, edgecolor='#B2E1FF', boxstyle='round', pad=0.2))
                # Periapsis marker
                self.ax.text(x_pe, y_pe, z_pe, 'Pe', horizontalalignment='center', verticalalignment='center', backgroundcolor='white', zorder=249,
                             fontdict={'family': 'serif', 'color':  'black', 'weight': 'bold', 'size': 8}).set_bbox(dict(facecolor='#B2E1FF', alpha=0.8, edgecolor='#B2E1FF', boxstyle='round', pad=0.2))
                # Ascending node marker
                self.ax.text(x_an, y_an, z_an, u"\u260A", horizontalalignment='center', verticalalignment='center', backgroundcolor='white', zorder=249,
                             fontdict={'family': "sans-serif", 'color':  'black', 'weight': 'black', 'size': 10}).set_bbox(dict(facecolor='#B2E1FF', alpha=0.8, edgecolor='#B2E1FF', boxstyle='round', pad=0.1))
                # Descending node marker
                self.ax.text(x_dn, y_dn, z_dn, u"\u260B", horizontalalignment='center', verticalalignment='center', backgroundcolor='white', zorder=249,
                             fontdict={'family': "sans-serif", 'color':  'black', 'weight': 'black', 'size': 10}).set_bbox(dict(facecolor='#B2E1FF', alpha=0.8, edgecolor='#B2E1FF', boxstyle='round', pad=0.1))
            # Current location marker◈◉◆▼⊛⊗
            if self.curr == True:
                self.ax.scatter(x_curr,y_curr,z_curr,color=curr_color,label=curr_label)
                if curr_background == False:
                    curr_fontdict = {'family': "sans-serif", 'color':  curr_color, 'weight': 'bold', 'size': 10}
                    self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                                 fontdict=curr_fontdict)
                else:
                    if curr_fontdict != None and type(curr_fontdict)==dict:
                        pass
                    else:
                        curr_fontdict = {'family': "sans-serif", 'color':  'white', 'weight': 'bold', 'size': 11}
                    if type(curr_color) == tuple:
                        self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                                     fontdict=curr_fontdict).set_bbox(dict(facecolor=curr_color, alpha=0.8, edgecolor=curr_color, boxstyle='round', pad=0.2))
                    elif curr_color in self.l_colors or '#' in curr_color:
                        self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                                     fontdict=curr_fontdict).set_bbox(dict(facecolor=curr_color, alpha=0.8, edgecolor=curr_color, boxstyle='round', pad=0.2))
                    else:
                        raise ValueError('Invalid current marker color')
                    #self.ax.scatter(x_curr, y_curr, z_curr, c='black', s=100, marker='X', zorder=1)
                
        if data == True:
            # General orbit details
            d1 = {
                'Designation'             : [self.designation],
                'Celestial Object'        : [self.body],
                'Body Radius'             : [dict_R[self.body]/1000.],
                'Body Mass'               : [self.M_b],
                'Gravitational Parameter' : [self.mu],
                'Unit'                    : ['km']
            }
            # Orbit
            d2 = {
                'Eccentricity'             : [self.e],
                'Semi-major Axis'          : [self.a],
                'Semi-minor Axis'          : [self.a*np.sqrt(1 - self.e**2)],
                'Semi-latus Rectum'        : [self.a*(1 - self.e**2)],
                'Period'                   : [self.T],
                'Apoapsis Radius'          : [self.Q],
                'Apoapsis Altitude'        : [self.h_Q],
                'Apoapsis Velocity'        : [V_a/1000.],
                'Periapsis Radius'         : [self.q],
                'Periapsis Altitude'       : [self.h_q],
                'Periapsis Velocity'       : [V_p/1000.],
                'Inclination'              : [np.rad2deg(self.i)],
                'Argument of periapsis'    : [np.rad2deg(self.w)],
                'Ascending node longitude' : [np.rad2deg(self.om)],
                'Unit'                     : [self.unit]
            }
            t = self.ta2t(ta1,a)
            # Coordinates
            d3 = {
                'Color'            : c1,
                'X'                : x1,
                'Y'                : y1,
                'Z'                : z1,
                'Radius'           : r1,
                'True Anomaly'     : np.rad2deg(ta1),
                'Time'             : t['dt'],
                'Revolution'       : np.ceil(t['mlt']),
                'Angular Velocity' : V1/r1, 
                'Velocity'         : V1
            }
            d = {
                'GNRL': pd.DataFrame(d1, index=['GNRL']),
                'ORBT': pd.DataFrame(d2, index=['ORBT']),
                'CRDS': pd.DataFrame(d3)
            }
            if self.ta != None:
                vx_temp = -(np.sqrt(self.mu*a)/(self.r*1000)) * np.sin(self.E)
                vy_temp = (np.sqrt(self.mu*a)/(self.r*1000)) * np.sqrt(1-self.e**2)*np.cos(self.E)
                vz_temp = 0

                vx_curr = vx_temp*np.cos(w)*np.cos(i)*np.cos(om) - vy_temp*np.sin(w)*np.cos(i)*np.cos(om) - vx_temp*np.sin(w)*np.sin(om) - vy_temp*np.cos(w)*np.sin(om)
                vy_curr = vx_temp*np.cos(w)*np.cos(i)*np.sin(om) - vy_temp*np.sin(w)*np.cos(i)*np.sin(om) + vx_temp*np.sin(w)*np.cos(om) + vy_temp*np.cos(w)*np.cos(om)
                vz_curr = -vx_temp*np.cos(w)*np.sin(i) + vy_temp*np.sin(w)*np.sin(i)
                d4 = {
                    'True Anomaly'        : np.rad2deg(self.ta),
                    'Radius'              : self.r,
                    'Eccentric Anomaly'   : np.rad2deg(self.E),
                    'Mean Anomaly'        : np.rad2deg(self.M),
                    'Mean Motion'         : self.n,
                    #'Time from Epoch'     : data_k['t_e'],
                    'Time from Periapsis' : self.tfp,
                    'Time to Periapsis'   : self.ttp,
                    'Velocity'            : self.V,
                    'Escape Velocity'     : self.V_esc,
                    'X'                   : x_curr,
                    'Y'                   : y_curr,
                    'Z'                   : z_curr,
                    'VX'                  : vx_curr/1000.,
                    'VY'                  : vy_curr/1000.,
                    'VZ'                  : vz_curr/1000.
                }
                d['SPEC'] = pd.DataFrame(d4, index=['SPEC'])
            return d
        
    def solarsystem(self, yr=2000, mo=1, dy=1, hr=0, mn=0, sec=0, now=False, time=None, sphere=False, unit='km', marker=True, curr_marker=False, view=True, res=None, linewidth=1.5):
        self.dCRDdict = {}
        try:
            marker = bool(marker)
        except:
            raise ValueError('Invalid marker bool')
        else:
            try:
                curr_marker = bool(curr_marker)
            except:
                raise ValueError('Invalid curr_marker bool')
            else:
                pass
        self.unit = unit
        if now:
            t = datetime.datetime.now().replace(tzinfo=tz.tzlocal()).astimezone(tz.gettz('UTC'))
        elif type(time)==datetime.datetime or type(time)==datetime.date:
            t = time.replace(tzinfo=tz.gettz('UTC'))
        else:
            t = datetime.datetime(yr,mo,dy,hr,mn,sec).replace(tzinfo=tz.gettz('UTC'))
        t_julian = t
        # SUN
        t = (t - datetime.datetime(2000, 1, 1, tzinfo=tz.gettz('UTC'))).total_seconds()
        R = conv.unit(dict_R['Sun'], 'm', self.unit)
        c = dict_c['Sun']
        if sphere==True and view==True:
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x_S = R*np.cos(u)*np.sin(v)
            y_S = R*np.sin(u)*np.sin(v)
            z_S = R*np.cos(v)
            self.ax.plot_wireframe(x_S, y_S, z_S, color=c, alpha=.5, zorder=0, label='_nolegend_')
        BODY = 'Sun'
        if view == True:
            self.ax.scatter(0, 0, 0, color=c, s=50, zorder=0, label=dict_sym[BODY] + u"\u202F" + BODY)
        # CENTURY TIME
        cent = 3155760000.0
        # ORBITAL DATA
        self.ssdata = []
        # PLANETS
        for PLNT in list_SS:
            T = t/cent
            a = conv.unit(dict_PLNTPOS[PLNT]['a'] + T*dict_PLNTPOS[PLNT]['da'], 'AU', self.unit)
            e = dict_PLNTPOS[PLNT]['e'] + T*dict_PLNTPOS[PLNT]['de']
            i = dict_PLNTPOS[PLNT]['I'] + T*dict_PLNTPOS[PLNT]['dI']
            om_bar = dict_PLNTPOS[PLNT]['o'] + T*dict_PLNTPOS[PLNT]['do']
            om = dict_PLNTPOS[PLNT]['O'] + T*dict_PLNTPOS[PLNT]['dO']
            w = om_bar - om
            L = dict_PLNTPOS[PLNT]['L'] + T*dict_PLNTPOS[PLNT]['dL']
            b = dict_PLNTPOS[PLNT]['b']
            c = dict_PLNTPOS[PLNT]['c']
            s = dict_PLNTPOS[PLNT]['s']
            f = dict_PLNTPOS[PLNT]['f']
            M = L - om_bar + b*T**2 + c*np.cos(np.deg2rad(f*T)) + s*np.sin(np.deg2rad(f*T))
            M_raw = M
            while M < -180:
                    M += 360
            while M > 180:
                    M += -360
            M = np.deg2rad(M)
            data = (M,e)
            E = optimize.fsolve(self.M_solve, 1, args=data)
            E = E[0]
            ta = np.rad2deg(2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)))
            self.orbit(a=a,
                       e=e,
                       i=i,
                       w=w,
                       om=om,
                       body='Sun',
                       ta=ta,
                       deg=True,
                       unit_i = self.unit,
                       unit_o = self.unit
            )
            if view == True:
                c = dict_c[PLNT]
            else:
                c = None
            dat = self.kepler(vis=c, data=True, marker=marker, curr_marker=curr_marker, curr_symbol=dict_sym[PLNT], curr_color=c, designation=PLNT, res=res, linewidth=linewidth)
            self.ssdata.append(dat)
            R = conv.unit(dict_R[PLNT], 'm', self.unit)
            self.dCRD['dX'][PLNT] = dat['SPEC']['X'][0]
            self.dCRD['dY'][PLNT] = dat['SPEC']['Y'][0]
            self.dCRD['dZ'][PLNT] = dat['SPEC']['Z'][0]
            if sphere==True and view==True:
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x_S = R*np.cos(u)*np.sin(v)
                y_S = R*np.sin(u)*np.sin(v)
                z_S = R*np.cos(v)
                for idx,val in enumerate(x_S):
                    for idx2,val2 in enumerate(x_S[idx]):
                        x_S[idx][idx2] += self.dCRD['dX'][PLNT]
                    for idx2,val2 in enumerate(y_S[idx]):
                        y_S[idx][idx2] += self.dCRD['dY'][PLNT]
                    for idx2,val2 in enumerate(z_S[idx]):
                        z_S[idx][idx2] += self.dCRD['dZ'][PLNT]
                self.ax.plot_wireframe(x_S, y_S, z_S, color=c, alpha=.5, zorder=0, label='_nolegend_')
            if view == True:
                self.ax.scatter(self.dCRD['dX'][PLNT], self.dCRD['dY'][PLNT], self.dCRD['dZ'][PLNT], color=c, s=50, zorder=0, label=dict_sym[PLNT] + ' ' + PLNT)
            #self.ax.plot(dat['CRDS']['X'], dat['CRDS']['Y'], dat['CRDS']['Z'], color=c, label='_nolegend_')
            self.ss = True
            self.epoch = t_julian
    
    def sgp4(self, l1, l2, time='now', orbit=False, marker=True, curr_background=True, curr_symbol='◉', curr_color='#92dd4f', curr_label='_nolegend_', curr_fontdict=None):
        if not sgp4_present:
            return warnings.warn("SPG4 not initialized: SGP4 propagation not functional", UserWarning)
        self.body = 'Earth'
        self.unit = 'km'
        if orbit==False:
            satellite = twoline2rv(l1, l2, wgs84)
            if time == 'now':
                currtime = pd.Timestamp.today(tz='Europe/Amsterdam').tz_convert('UTC')
            else:
                currtime = time
            currpos, currvel = satellite.propagate(currtime.year, currtime.month, currtime.day, currtime.hour, currtime.minute, currtime.second)
            x_curr = currpos[0]
            y_curr = currpos[1]
            z_curr = currpos[2]
            
        if marker == True:
            self.ax.scatter(x_curr,y_curr,z_curr,color=curr_color,label=curr_label)
            if curr_background == False:
                curr_fontdict = {'family': "sans-serif", 'color':  curr_color, 'weight': 'bold', 'size': 11}
                self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                             fontdict=curr_fontdict)
            else:
                if curr_fontdict != None and type(curr_fontdict)==dict:
                    pass
                else:
                    curr_fontdict = {'family': "sans-serif", 'color':  'white', 'weight': 'bold', 'size': 11}
                if type(curr_color) == tuple:
                    self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                                 fontdict=curr_fontdict).set_bbox(dict(facecolor=curr_color, alpha=0.8, edgecolor=curr_color, boxstyle='round', pad=0.2))
                elif curr_color in self.l_colors or '#' in curr_color:
                    self.ax.text(x_curr, y_curr, z_curr, curr_symbol, horizontalalignment='center', verticalalignment='center', zorder=250,
                                 fontdict=curr_fontdict).set_bbox(dict(facecolor=curr_color, alpha=0.8, edgecolor=curr_color, boxstyle='round', pad=0.2))
                else:
                    raise ValueError('Invalid current marker color')
        
    def getTLE(self):
        if not pycurl_present:
            return warnings.warning('Pycurl not initialized, getTLE function cannot be executed',UserWarning)
        USR = str(os.environ['USR'])
        PWD = str(os.environ['PWD'])

        # Current TLE
        QRY = "https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/>now-30/format/3le"
        c = pycurl.Curl()
        c.setopt(c.URL, 'https://www.space-track.org/ajaxauth/login')
        c.setopt(c.POSTFIELDS, 'identity='+USR+'&password='+PWD+'&query='+QRY)
        with BytesIO() as e:
            c.setopt(c.WRITEFUNCTION, e.write)
            c.setopt(pycurl.SSL_VERIFYPEER, 0)
            c.perform()
            c.close()
            filename = 'assets/current.3le'
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            f = open(filename, 'w', newline='')
            f.write(e.getvalue().decode('UTF-8'))

        # Recent SATS
        QRY = "https://www.space-track.org/basicspacedata/query/class/satcat/LAUNCH/>now-7/CURRENT/Y/orderby/LAUNCH DESC/format/csv"
        c = pycurl.Curl()
        c.setopt(c.URL, 'https://www.space-track.org/ajaxauth/login')
        c.setopt(c.POSTFIELDS, 'identity='+USR+'&password='+PWD+'&query='+QRY)
        with BytesIO() as e:
            c.setopt(c.WRITEFUNCTION, e.write)
            c.setopt(pycurl.SSL_VERIFYPEER, 0)
            c.perform()
            c.close()
            filename = 'assets/satcatrecent.csv'
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            f = open(filename, 'w', newline='')
            f.write(e.getvalue().decode('UTF-8'))
        
        # Current SATS
        QRY = "https://www.space-track.org/basicspacedata/query/class/satcat/CURRENT/Y/orderby/LAUNCH DESC/format/csv"
        c = pycurl.Curl()
        c.setopt(c.URL, 'https://www.space-track.org/ajaxauth/login')
        c.setopt(c.POSTFIELDS, 'identity='+USR+'&password='+PWD+'&query='+QRY)
        with BytesIO() as e:
            c.setopt(c.WRITEFUNCTION, e.write)
            c.setopt(pycurl.SSL_VERIFYPEER, 0)
            c.perform()
            c.close()
            filename='assets/satcatcurrent.csv'
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            f = open(filename, 'w', newline='')
            f.write(e.getvalue().decode('UTF-8'))

        # Launchsites
        QRY = "https://www.space-track.org/basicspacedata/query/class/launch_site/format/csv"
        c = pycurl.Curl()
        c.setopt(c.URL, 'https://www.space-track.org/ajaxauth/login')
        c.setopt(c.POSTFIELDS, 'identity='+USR+'&password='+PWD+'&query='+QRY)
        with BytesIO() as e:
            c.setopt(c.WRITEFUNCTION, e.write)
            c.setopt(pycurl.SSL_VERIFYPEER, 0)
            c.perform()
            c.close()
            filename = 'assets/launchsite.csv'
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            f = open(filename, 'w', newline='')
            f.write(e.getvalue().decode('UTF-8'))

    def parseTLE(self):
        self.SATCURR = pd.read_csv('assets/satcatcurrent.csv')
        self.LS = pd.read_csv('assets/launchsite.csv') 
        self.SATCURR.rename(columns={'NORAD_CAT_ID':'SCN'}, inplace=True)
        self.LS['LAUNCH_SITE'] = [x.upper() for x in self.LS['LAUNCH_SITE']]

        # Parse TLEs
        with open('assets/current.3le') as f:
            TLE = f.readlines()
        TLE = [x.strip() for x in TLE]
        TLECAT = {
            'SATNAME' : [],
            'SCN'     : [],
            'ELSET'   : [],
            'INTLDES' : [],
            'LINE1'   : [],
            'LINE2'   : [],
        }
        CAT = {
            'SATNAME'   : [],
            'SCN'    : [],
            'SCN2'   : [],
            'ELSET'  : [],
            'INTLDES' : [],
            'EPOCH'  : [],
            'YR'     : [],
            'DY'     : [],
            'DYFR'   : [],
            'N'      : [],
            'NS'     : [],
            'N1'     : [],
            'N1S'    : [],
            'N2'     : [],
            'N2S'    : [],
            'BSTAR'  : [],
            'ELTYPE' : [],
            'ELNUM'  : [],
            'CHKSM'  : [],
            'CHKSM2' : [],
            'IN'     : [],
            'OM'     : [],
            'EC'     : [],
            'W'      : [],
            'MA'     : [],
            'REV'    : []
        }
        RAWCAT = {
            'SATNAME'   : [],
            'SCN'    : [],
            'SCN2'   : [],
            'ELSET'  : [],
            'INTLDES' : [],
            'EPOCH'  : [],
            'YR'     : [],
            'DY'     : [],
            'DYFR'   : [],
            'N'      : [],
            'N1'     : [],
            'N2'     : [],
            'BSTAR'  : [],
            'ELTYPE' : [],
            'ELNUM'  : [],
            'CHKSM'  : [],
            'CHKSM2' : [],
            'IN'     : [],
            'OM'     : [],
            'EC'     : [],
            'W'      : [],
            'MA'     : [],
            'REV'    : []
        }
        SATCAT = {'SATNAME': [], 'SCN': [], 'ELSET': [], 'SCN2': [], 'INTLDES': []}
        LN = 0
        while LN < len(TLE):
            if int(TLE[LN][0]) == 0:
                SATNAME0 = TLE[LN][2:]
                RAWCAT['SATNAME'].append(SATNAME0)
                SATNAME = SATNAME0.strip()
                CAT['SATNAME'].append(SATNAME)
                SATCAT['SATNAME'].append(SATNAME)
                TLECAT['SATNAME'].append(SATNAME)
            elif int(TLE[LN][0]) == 1:
                LINE1 = TLE[LN]
                SCN0 = TLE[LN][2:7]
                SCN = int(SCN0)
                ELSET = TLE[LN][7]
                INTLDES0 = TLE[LN][9:17]
                INTLDES = INTLDES0.strip()
                EPOCH0 = TLE[LN][18:32]
                EPOCH = float(EPOCH0)
                YR0 = TLE[LN][18:20]
                YR = int(YR0)
                DY0 = TLE[LN][20:23]
                DY = int(DY0)
                DYFR0 = TLE[LN][23:32]
                DYFR = float(DYFR0)
                N01 = TLE[LN][33:43]
                N1 = float(N01)*2
                N1S = N1/86400.0
                N02 = TLE[LN][44:52]
                if len(N02[1:].split('-')) > 1:
                    N2 = float(N02[0] + '0.' + N02[1:].split('-')[0] + 'e-' + N02[1:].split('-')[1])*6
                elif len(N02[1:].split('+')) > 1:
                    N2 = float(N02[0] + '0.' + N02[1:].split('+')[0] + 'e' + N02[1:].split('+')[1])*6
                N2S = N2/86400.0
                BSTAR0 = TLE[LN][53:61]
                if len(BSTAR0[1:].split('-')) > 1:
                    BSTAR = float(BSTAR0[0] + '0.' + BSTAR0[1:].split('-')[0] + 'e-' + BSTAR0[1:].split('-')[1])
                elif len(BSTAR0[1:].split('+')) > 1:
                    BSTAR = float(BSTAR0[0] + '0.' + BSTAR0[1:].split('+')[0] + 'e' + BSTAR0[1:].split('+')[1])
                ELTYPE0 = TLE[LN][62]
                ELTYPE = int(ELTYPE0)
                ELNUM0 = TLE[LN][64:68]
                ELNUM = int(ELNUM0)
                CHKSM0 = TLE[LN][68]
                CHKSM = int(CHKSM0)
                # ADD TO CAT, RAW/TLECAT
                TLECAT['LINE1'].append(LINE1)
                RAWCAT['SCN'].append(SCN0)
                CAT['SCN'].append(SCN)
                TLECAT['SCN'].append(SCN)
                RAWCAT['ELSET'].append(ELSET)
                CAT['ELSET'].append(ELSET)
                TLECAT['ELSET'].append(ELSET)
                RAWCAT['INTLDES'].append(INTLDES0)
                CAT['INTLDES'].append(INTLDES)
                TLECAT['INTLDES'].append(INTLDES)
                RAWCAT['EPOCH'].append(EPOCH0)
                CAT['EPOCH'].append(EPOCH)
                RAWCAT['YR'].append(YR0)
                CAT['YR'].append(YR)
                RAWCAT['DY'].append(DY0)
                CAT['DY'].append(DY)
                RAWCAT['DYFR'].append(DYFR0)
                CAT['DYFR'].append(DYFR)
                RAWCAT['N1'].append(N01)
                CAT['N1'].append(N1)
                CAT['N1S'].append(N1S)
                RAWCAT['N2'].append(N02)
                CAT['N2'].append(N2)
                CAT['N2S'].append(N2S)
                RAWCAT['BSTAR'].append(BSTAR0)
                CAT['BSTAR'].append(BSTAR)
                RAWCAT['ELTYPE'].append(ELTYPE0)
                CAT['ELTYPE'].append(ELTYPE)
                RAWCAT['ELNUM'].append(ELNUM0)
                CAT['ELNUM'].append(ELNUM)
                RAWCAT['CHKSM'].append(CHKSM0)
                CAT['CHKSM'].append(CHKSM)
                # ADD TO SATCAT
                SATCAT['SCN'].append(SCN)
                SATCAT['ELSET'].append(ELSET)
                SATCAT['INTLDES'].append(INTLDES)
            elif int(TLE[LN][0]) == 2:
                LINE2 = TLE[LN]
                SCN02 = TLE[LN][2:7]
                SCN2 = int(SCN02)
                IN0 = TLE[LN][8:16]
                IN = float(IN0)
                OM0 = TLE[LN][17:25]
                OM = float(OM0)
                EC0 = TLE[LN][26:33]
                EC = float('0.' + EC0)
                W0 = TLE[LN][34:42]
                W = float(W0)
                MA0 = TLE[LN][43:51]
                MA = float(MA0)
                N0 = TLE[LN][52:63]
                N = float(N0)
                NS = N/86400.0
                REV0 = TLE[LN][63:68]
                REV = float(REV0)
                CHKSM02 = TLE[LN][68]
                CHKSM2 = int(CHKSM02)
                # APPEND TO CAT & RAWCAT
                TLECAT['LINE2'].append(LINE2)
                RAWCAT['SCN2'].append(SCN02)
                CAT['SCN2'].append(SCN2)
                RAWCAT['IN'].append(IN0)
                CAT['IN'].append(IN)
                RAWCAT['OM'].append(OM0)
                CAT['OM'].append(OM)
                RAWCAT['EC'].append(EC0)
                CAT['EC'].append(EC)
                RAWCAT['W'].append(W0)
                CAT['W'].append(W)
                RAWCAT['MA'].append(MA0)
                CAT['MA'].append(MA)
                RAWCAT['N'].append(N0)
                CAT['N'].append(N)
                CAT['NS'].append(NS)
                RAWCAT['REV'].append(REV0)
                CAT['REV'].append(REV)
                RAWCAT['CHKSM2'].append(CHKSM02)
                CAT['CHKSM2'].append(CHKSM2)
                # APPEND TO SATCAT
                SATCAT['SCN2'].append(SCN2)
            LN += 1

            
        self.CAT = pd.DataFrame(CAT)
        self.RAWCAT = pd.DataFrame(RAWCAT)
        self.SATCAT = pd.DataFrame(SATCAT)
        self.TLECAT = pd.DataFrame(TLECAT)

    def satorbit(self,SATNAME=None, SCN=None, INTLDES=None, cont=False, nocont=None, time='now', marker=False, curr_marker=True, curr_symbol='◉', vline=False, nums=True, abbrev=False, info=True, legend=True, label='SATNAME', res=25, rev=1, orbit=False, bg=True, orbitcolor=None, markercolor=None, SGP=False, linewidth=1.5):
        if SATNAME != None:
            QRYTYP = 'SATNAME'
            QRY = SATNAME
            try:
                cont = bool(cont)
            except:
                raise ValueError('Invalid cont bool')
            else:
                if cont==True:
                    if type(SATNAME)==list:
                        SATNAME = [x.upper() for x in SATNAME] #map(lambda x: x.upper(), SATNAME)
                        QRY = SATNAME
                        SATS=self.CAT
                        SATcontS = []
                        for X in SATNAME:
                            SATcontTMP = SATS[SATS['SATNAME'].str.contains(X)]
                            SATcontS.append(SATcontTMP)
                        SATS = pd.concat(SATcontS)
                    else:
                        SATNAME=SATNAME.upper()
                        QRY = SATNAME
                        SATS = self.CAT[self.CAT['SATNAME'].str.contains(SATNAME)]
                    if nocont != None:
                        if type(nocont)==list:
                            nocontS = SATS
                            nocont = [x.upper() for x in nocont] #map(lambda x: x.upper(), nocont)
                            for X in nocont:
                                nocontS = nocontS[~nocontS['SATNAME'].str.contains(X)]
                            SATS = nocontS
                        else:    
                            nocont=nocont.upper()
                            nocontS = self.CAT[~self.CAT['SATNAME'].str.contains(nocont)]
                    else:
                        pass
                elif cont==False:
                    if type(SATNAME)==list:
                        SATNAME = [x.upper() for x in SATNAME] #map(lambda x: x.upper(), SATNAME)
                        QRY = SATNAME
                        SATS=self.CAT
                        SATcontS = []
                        for X in SATNAME:
                            SATcontTMP = SATS[SATS['SATNAME']==X]
                            SATcontS.append(SATcontTMP)
                        SATS = pd.concat(SATcontS)
                    else:
                        SATNAME=SATNAME.upper()
                        QRY = SATNAME
                        SATS = self.CAT[self.CAT['SATNAME']==SATNAME]
        elif SCN != None:
            QRYTYP = 'SCN'
            QRY = SCN
            if type(SCN)==list:
                SATS=self.CAT
                SATcontS = []
                for X in SCN:
                    try:
                        X = int(X)
                    except:
                        raise ValueError('Invalid SCN')
                    else:
                        SATcontTMP = SATS[SATS['SCN']==X]
                        SATcontS.append(SATcontTMP)
                QRY = [x.upper() for x in QRY] #map(lambda x: str(x), QRY)
                SATS = pd.concat(SATcontS)
            else:
                try:
                    SCN = int(SCN)
                except:
                    raise ValueError('Invalid SCN')
                else:
                    SATS = self.CAT[self.CAT['SCN']==SCN]
        elif INTLDES != None:
            QRYTYP = 'INTLDES'
            QRY = INTLDES
            if type(INTLDES)==list:
                SATS=self.CAT
                SATcontS = []
                for X in INTLDES:
                    SATcontTMP = SATS[SATS['INTLDES']==X]
                    SATcontS.append(SATcontTMP)
                SATS = pd.concat(SATcontS)
            else:
                SATS = self.CAT[self.CAT['INTLDES']==INTLDES]
            
        if len(SATS) == 0:
            raise ValueError('No results for query')

        if time != None:
            if time == 'now':
                time = datetime.datetime.now().replace(tzinfo=tz.tzlocal())
            elif type(time) == tuple:
                try:
                    YR,MO,DY,HR,MIN,SEC=time
                    time = datetime.datetime(YR, MO, DY, HR, MIN, SEC, tzinfo=UTC)
                except:
                    raise ValueError('Invalid time')
                else:
                    pass
            else:
                raise ValueError('Invalid time')
        else:
            raise ValueError('time not specified')
        try:
            str(label)
        except:
            raise ValueError('Invalid label')
        else:
            if label == 'SATNAME':
                lbl = 0
            elif label == 'SCN':
                lbl = 1
            elif label == 'INTLDES':
                lbl = 2
        NO = 0
        SATDAT = {}
        self.res=res
        colors = plt.get_cmap("hsv")
        COLN = 1
        if len(SATS) > 10:
            colors = plt.get_cmap("hsv")
            clrdiv = len(SATS)
            if len(SATS) > 40:
                COLN = int(np.ceil(len(SATS)/40.))
        else:
            colors = plt.get_cmap("tab10")
            clrdiv = 10
            COLN = 1
        #plt.get_cmap("tab20")
        while NO < len(SATS):
            SAT = SATS.iloc[NO]
            SATNAME = SAT['SATNAME']
            SCN = SAT['SCN']
            INTLDES = SAT['INTLDES']
            YR = SAT['YR']
            DY = SAT['DY']
            SEC = SAT['DYFR']*86400.0
            UTC = tz.gettz('UTC')
            if YR < 59:
                FULLYEAR = 2000 + YR
            else:
                FULLYEAR = 1900 + YR
            EPOCH = datetime.datetime(FULLYEAR, 1, 1, tzinfo=UTC) + datetime.timedelta(int(DY)-1,SEC)
            if orbit == True:
                # Argument of Perigee
                W = SAT['W']
                # Right ascension
                OM = SAT['OM']
                # Mean Anomaly
                MA = SAT['MA']
                # Inclination
                IN = SAT['IN']
                # Eccentricity
                EC = SAT['EC']
                # Mean Motion
                NS = SAT['NS']
                N = SAT['N']
                # 1st time Der. of Mean Motion
                N1S = SAT['N1S']
                N1 = SAT['N1']
                # 2nd time Der. of Mean Motion
                N2S = SAT['N2S']
                # Pseudo-drag term B*
                BSTAR = SAT['BSTAR']

                time = datetime.datetime.now().replace(tzinfo=tz.tzlocal())
                DT = (time-EPOCH).total_seconds()/86400.0
                NCURR = N + N1*DT
                # GRAVITATIONAL PARAM KM3/DY2
                MU = 2.97554e15
                def MACURR(DT,MA,NCURR):
                    return MA + 360*(NCURR*DT-int(NCURR*DT)-int((MA+360*(NCURR*DT-int(NCURR*DT)))/360.))
                MACURR = np.deg2rad(MACURR(DT,MA,NCURR))
                def M_solve(E,*data):
                    M, e = data
                    return E - e*np.sin(E) - M
                data = (MACURR,EC)
                E = optimize.fsolve(M_solve, 1, args=data)
                E = E[0]
                TA = np.rad2deg(2*np.arctan(np.sqrt((1+EC)/(1-EC))*np.tan(E/2)))
                while TA < 0:
                    TA += 360
                while TA>360:
                    TA += -360
                A = (MU/((2*np.pi*NCURR)**2))**(1./3.)

                labels = [SATNAME,SCN,INTLDES]
                lblTYP = ['SATNAME','SCN','INTLDES'][lbl]
                curr_label = str(labels[lbl])
                if nums == True:
                    curr_symbol=str(NO)
                    curr_label = str(NO) + u' \u2014 ' + str(labels[lbl])
                    if len(SATS)>10:
                        if 0.1 <= (NO/float(clrdiv)) <= 0.55:
                            curr_fontdict = {'family': "monospace", 'color':  'black', 'weight': 'bold', 'size': 9}
                        else:
                            curr_fontdict = {'family': "monospace", 'color':  'white', 'weight': 'bold', 'size': 9}
                    else:
                        curr_fontdict = {'family': "monospace", 'color':  'white', 'weight': 'bold', 'size': 9}
                elif nums == False:
                    curr_fontdict=None
                else:
                    raise ValueError('Invalid nums bool')
                self.orbit(a=A,
                        e=EC,
                        ta=TA,
                        i=IN,
                        om=OM,
                        w=W,
                        unit_i='km',
                        unit_o='km',
                        body='Earth')
                try:
                    vline = bool(vline)
                except:
                    raise ValueError('Invalid vline bool')
                
                if SGP:
                    L1 = self.TLECAT[self.CAT['SCN']==SCN]['LINE1'].iloc[0]
                    L2 = self.TLECAT[self.CAT['SCN']==SCN]['LINE2'].iloc[0]
                else:
                    L1 = None
                    L2 = None
                
                if vline==True and len(SATS)==1:
                    curr_label = str(curr_symbol) + u' \u2014 ' + str(labels[lbl])
                    CURRDAT = self.kepler(vis='v-line', marker=marker, curr_background=bg, curr_marker=curr_marker, curr_symbol=curr_symbol, curr_fontdict=curr_fontdict, curr_label=curr_label, data=True, rev=rev, sgp=SGP, l1=L1, l2=L2, time=time, linewidth=linewidth)
                else:
                    orbitclr = colors(NO/float(clrdiv))
                    markerclr = colors(NO/float(clrdiv))
                    if orbitcolor != None:
                        orbitclr = orbitcolor
                    if markercolor != None:
                        markerclr = markercolor                   
                    CURRDAT = self.kepler(vis=orbitclr, marker=marker, curr_background=bg, curr_marker=curr_marker, curr_symbol=curr_symbol, curr_fontdict=curr_fontdict, curr_color=markerclr, curr_label=curr_label, data=True, rev=rev, sgp=SGP, l1=L1, l2=L2, time=time, linewidth=linewidth)
                SATDAT[NO] = CURRDAT
            else:
                labels = [SATNAME,SCN,INTLDES]
                lblTYP = ['SATNAME','SCN','INTLDES'][lbl]
                curr_label = str(labels[lbl])
                if nums == True:
                    curr_symbol=str(NO)
                    curr_label = str(NO) + u' \u2014 ' + str(labels[lbl])
                    if len(SATS)>10:
                        if 0.1 <= (NO/float(clrdiv)) <= 0.55:
                            curr_fontdict = {'family': "monospace", 'color':  'black', 'weight': 'bold', 'size': 9}
                        else:
                            curr_fontdict = {'family': "monospace", 'color':  'white', 'weight': 'bold', 'size': 9}
                    else:
                        curr_fontdict = {'family': "monospace", 'color':  'white', 'weight': 'bold', 'size': 9}
                elif nums == False:
                    curr_symbol='◉'
                    curr_fontdict=None
                else:
                    raise ValueError('Invalid nums bool')
                L1 = self.TLECAT[self.CAT['SCN']==SCN]['LINE1'].iloc[0]
                L2 = self.TLECAT[self.CAT['SCN']==SCN]['LINE2'].iloc[0]
                markerclr = colors(NO/float(clrdiv))
                if markercolor != None:
                    markerclr = markercolor  
                self.sgp4(L1, L2, time=time, marker=True, curr_background=bg, curr_symbol=curr_symbol, curr_fontdict=curr_fontdict, curr_color=markerclr, curr_label=curr_label)
                CURRDAT = None
                SATDAT[NO] = {'SCN':SCN,'INTLDES':INTLDES,'SATNAME':SATNAME}
                
            NO += 1
            
        SATDAT = pd.DataFrame(SATDAT).transpose()
        if legend == True:
            self.ax.legend(loc=1,prop={"family":"monospace"},ncol=COLN).set_zorder(500)
        elif legend == False:
            pass
        else:
            raise ValueError('Invalid legend bool')
        if len(SATS) > 1:
            RANGE = str(SATS.iloc[0][lblTYP])+u' \u2014 '+str(SATS.tail(1)[lblTYP].iloc[0])
        else:
            RANGE = 'NAN'
        if type(QRY)==list:
            QRY = '['+','.join(QRY)+']'
        if type(nocont)==list:
            nocont = '['+','.join(nocont)+']'
        if info == True:
            if len(SATS)==1:
                CURRSATINFO = self.SATCURR[self.SATCURR['SCN'] == SCN]
                OBJECT_TYPE = CURRSATINFO['OBJECT_TYPE'].iloc[0]
                COUNTRY = CURRSATINFO['COUNTRY'].iloc[0]
                LAUNCH = CURRSATINFO['LAUNCH'].iloc[0]
                SITE_CODE = CURRSATINFO['SITE'].iloc[0]
                LAUNCH_SITE = self.LS[self.LS['SITE_CODE']==SITE_CODE]['LAUNCH_SITE'].iloc[0]
                LAUNCH_PIECE = CURRSATINFO['LAUNCH_PIECE'].iloc[0]
                if abbrev==True:
                    LS_STR = SITE_CODE
                elif abbrev==False:
                    LS_STR = LAUNCH_SITE+' ('+SITE_CODE+')'
                else:
                    raise ValueError('Invalid abbrev bool')
                SATINFO = 'SATNAME : {}\nSCN     : {}\nINTLDES : {}\nOBJ TYPE: {}\nPIECE   : {}\nCOUNTRY : {}\nLAUNCH  : {}\nSITE    : {}'.format(SATNAME,SCN,INTLDES,OBJECT_TYPE,LAUNCH_PIECE,COUNTRY,LAUNCH,LS_STR)
                self.ax.text2D(0.025, 0.735, 'SAT INFO', transform=self.ax.transAxes, va='top',zorder=501,
                            fontdict={'family': 'monospace', 'color':  'black', 'weight': 'bold', 'size': 12}).set_bbox(dict(facecolor='#EEEEEE', alpha=0.75, edgecolor='#DBDBDB', boxstyle='round', pad=0.2))
                self.ax.text2D(0.025, 0.705, SATINFO, transform=self.ax.transAxes, va='top',zorder=501,
                                fontdict={'family': 'monospace', 'color':  'black', 'weight': 'normal', 'size': 10}).set_bbox(dict(facecolor='#EEEEEE', alpha=0.75, edgecolor='#DBDBDB', boxstyle='round', pad=0.2))
            QRYINFO = 'DATE   : {}\nTIME   : {}\nEPOCH  : {}\nQUERY  : {}\nTYPE   : {}\nCONT.  : {}\nNOCONT.: {}\nRANGE  : {}\nNUMBER : {}'.format(time.strftime('%Y-%m-%d'), time.strftime('%H:%M:%S:%f %Z'), EPOCH.strftime('%Y-%m-%d %H:%M:%S'), QRY, QRYTYP, str(cont).upper(), str(nocont).upper(), RANGE, NO)
            self.ax.text2D(0.025, 0.975, 'QUERY INFO', transform=self.ax.transAxes, va='top',zorder=501,
                            fontdict={'family': 'monospace', 'color':  'black', 'weight': 'bold', 'size': 12}).set_bbox(dict(facecolor='#EEEEEE', alpha=0.75, edgecolor='#DBDBDB', boxstyle='round', pad=0.2))
            self.ax.text2D(0.025, 0.94, QRYINFO, transform=self.ax.transAxes, va='top',zorder=501,
                            fontdict={'family': 'monospace', 'color':  'black', 'weight': 'normal', 'size': 10}).set_bbox(dict(facecolor='#EEEEEE', alpha=0.75, edgecolor='#DBDBDB', boxstyle='round', pad=0.2))
        return {'SATS':SATS,'SATDAT':SATDAT,'ORBTDAT':CURRDAT}

    def view(self, x=9, y=9, el=30, az=45, cust_txt=False, bodysurface=False, bodycolor=None, bodyalpha=0.5, show=True, save=False, name='EOS Output', extension='pdf',legend=False):
        try:
            bodyalpha = float(bodyalpha)
        except:
            raise ValueError('Invalid bodyalpha')
        else:
            if 0 <= bodyalpha <=1:
                pass
            else:
                raise ValueError('Invalid bodyalpha')
        try:
            bodysurface = bool(bodysurface)
        except:
            raise ValueError('Invalid bodysurface bool')
        self.cust_txt = cust_txt
        #plt.gca().invert_xaxis()
        #plt.gca().invert_yaxis()
        # Draw body
        if self.ss == False:
            R = conv.unit(dict_R[self.body], 'm', self.unit)
            if bodycolor != None:
                if (type(bodycolor) == tuple and len(bodycolor)==3) or bodycolor in self.l_colors or '#' in bodycolor:
                    c = bodycolor
            else:
                c = dict_c[self.body]
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x_S = R*np.cos(u)*np.sin(v)
            y_S = R*np.sin(u)*np.sin(v)
            z_S = R*np.cos(v)
            if bodysurface==False:
                self.ax.plot_wireframe(x_S, y_S, z_S, color=c, alpha=bodyalpha, label='_nolegend_', zorder=0)
            elif bodysurface==True:
                self.ax.plot_surface(x_S, y_S, z_S, color=c, alpha=bodyalpha, label='_nolegend_', zorder=0)
            x_limits = self.ax.get_xlim3d()
            y_limits = self.ax.get_ylim3d()
            z_limits = self.ax.get_zlim3d()
            x_range = abs(x_limits[1] - x_limits[0])
            x_middle = np.mean(x_limits)
            y_range = abs(y_limits[1] - y_limits[0])
            y_middle = np.mean(y_limits)
            z_range = abs(z_limits[1] - z_limits[0])
            z_middle = np.mean(z_limits)
            plot_radius = 0.5*max([x_range, y_range, z_range])
            self.ax.set_xlim3d([-plot_radius, plot_radius])
            self.ax.set_ylim3d([-plot_radius, plot_radius])
            self.ax.set_zlim3d([-plot_radius, plot_radius])
            if cust_txt == False:
                titletxt = 'Orbit'
                self.ax.set_title(r"$\mathrm{" + self.body + " \ " + titletxt + " \ Visualization}$", va='bottom')      
        elif self.ss == True:
            # (dX,dY,dZ) as center of plot
            x_limits = self.ax.get_xlim3d()
            y_limits = self.ax.get_ylim3d()
            z_limits = self.ax.get_zlim3d()
            x_range = abs(x_limits[1] - x_limits[0])
            x_middle = self.dCRD['dX'][self.center]
            y_range = abs(y_limits[1] - y_limits[0])
            y_middle = self.dCRD['dY'][self.center]
            z_range = abs(z_limits[1] - z_limits[0])
            z_middle = self.dCRD['dZ'][self.center]
            plot_radius = 0.5*max([x_range, y_range, z_range])
            self.ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
            self.ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
            self.ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
            if cust_txt == False:
                self.fig.suptitle(r'Solar System', fontsize=14, y=0.87, x=0.515)
                self.ax.set_title(r'{}'.format(self.epoch.strftime("%Y-%m-%d %H:%M:%S %Z")), fontsize=12, x=0.5)
        if type(cust_txt) == str:
            self.fig.suptitle(cust_txt, fontsize=12, y=0.965, x=0.515)
        if legend:
            self.ax.legend(loc=1,prop={"family":"sans\\-serif"}).set_zorder(500)
        self.ax.set_aspect(1)
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0,0))
        self.ax.xaxis.set_major_formatter(formatter) 
        self.ax.yaxis.set_major_formatter(formatter) 
        self.ax.zaxis.set_major_formatter(formatter)
        self.ax.set_xlabel(r'$X \ \mathrm{['+ self.unit +']}$')
        self.ax.set_ylabel(r'$Y \ \mathrm{['+ self.unit +']}$')
        self.ax.set_zlabel(r'$Z \ \mathrm{['+ self.unit +']}$')
        self.ax.view_init(el,az)
        self.fig.set_size_inches(x, y)
        #plt.tight_layout()
        if save == True:
            plt.savefig('{}.{}'.format(name,extension))
        if show == True:
            plt.show()