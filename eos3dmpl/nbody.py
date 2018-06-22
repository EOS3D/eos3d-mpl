import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from datetime import datetime, timedelta

from eos3dmpl.eos_aedata import *
from eos3dmpl.eos_core import *

class Body():
        def __init__(self,pos,vel,radius,mass,name,color='#e84d0b'):
            self.t = 0
            self.pos0 = pos
            self.posx = [pos[0]]
            self.posy = [pos[1]]
            self.posz = [pos[2]]
            self.pos = pos
            self.vpos = np.array(pos)
            self.vpos0 = np.array(pos)
            self.vel = vel
            self.vel0 = vel
            self.radius = radius
            self.mass = mass
            self.objs = []
            self.impulse = []
            #self.size = 50*radius/dict_R['Sun']
            self.size = 15
            self.name = name
            self.color = color
        
        def reset(self):
            self.posx = [self.pos0[0]]
            self.posy = [self.pos0[1]]
            self.posz = [self.pos0[2]]
            self.vel = self.vel0
            self.objs = []
            self.impulse = []

        def addobj(self,obj):
            if type(obj) == list:
                self.objs.extend(obj)
            else:
                self.objs.append(obj)

        def addimpulse(self,dv,time):
            self.impulse.append([dv,time,False])

class Nbody():
    G = 6.67408e-11 # m^3kg^-1s^-2
    AU = 149597870700 # m
    objs = []
    dots = []
    traces = []
    drawings = []

    def __init__(self,frames=1000,timestep=3600,frame0=0,interval=1,lims=[[-1.5*149597870700,1.5*149597870700],[-1.5*149597870700,1.5*149597870700]], yr=1, mo=1, dy=1, hr=0, mn=0, sec=0, time=None):
        self.frames = frames
        self.timestep = timestep
        self.frame0 = frame0
        self.interval = interval
        self.lims = lims
        self.yr = yr
        self.mo = mo
        self.dy = dy
        self.hr = hr
        self.mn = mn
        self.sec = sec
        if (type(time)==datetime.datetime or type(time)==datetime.date):
            self.time = time
            self.manual = True
        elif time=='manual':
            self.time = datetime.datetime(self.yr, self.mo, self.dy, self.hr, self.mn, self.sec)
            self.manual = True
        else:
            self.time = datetime.datetime(1,1,1,0,0,0)
            self.manual = False
            
    def reset(self):
        for obj in self.objs:
            obj.reset()
        self.objs = []
        self.dots = []
        self.traces = []
        self.drawings = []

    def gettime(self,t):
        d = self.time + timedelta(seconds=t)
        if self.manual:
            return '{: <2}yr:{: <2}mo:{: <2}dy:{: <2}hr:{: <2}mn:{: <2}s'.format(d.year,d.month,d.day,d.hour,d.minute,d.second)
        else:
            return '{: <2}yr:{: <2}mo:{: <2}dy:{: <2}hr:{: <2}mn:{: <2}s'.format(d.year-1,d.month-1,d.day-1,d.hour,d.minute,d.second)

    def propagate(self):
        G = self.G
        frames = self.frames + self.frame0
        dt = self.timestep
        for frame in range(frames):
            for obj in self.objs:
                obj.acc = [0.,0.,0.]
                obj.t = frame*dt
                for obj2 in self.objs:
                    relpos = obj2.vpos - obj.vpos
                    rdist = np.linalg.norm(relpos)
                    if rdist == 0:
                        pass
                    else:
                        vrdist = (relpos/rdist)
                        acc = vrdist*(G*obj2.mass/(rdist**2))
                        obj.acc[0] += acc[0]
                        obj.acc[1] += acc[1]
                        obj.acc[2] += acc[2]

            for obj in self.objs:
                obj.vel[0] += obj.acc[0]*dt
                obj.vel[1] += obj.acc[1]*dt
                obj.vel[2] += obj.acc[2]*dt

                if len(obj.impulse)>0:
                    for idx,impulse in enumerate(obj.impulse):
                        if impulse[1] <= obj.t and impulse[2] == False:
                            obj.vel[0] += impulse[0][0]
                            obj.vel[1] += impulse[0][1]
                            obj.vel[2] += impulse[0][2]
                            obj.impulse[idx][2] = True
                            

                obj.posx.append(obj.pos[0])
                obj.posy.append(obj.pos[1])
                obj.posz.append(obj.pos[2])

                obj.pos[0] += obj.vel[0]*dt
                obj.pos[1] += obj.vel[1]*dt
                obj.pos[2] += obj.vel[2]*dt

                obj.vpos = np.array(obj.pos)
    
    def addobj(self,obj):
        if type(obj) == list:
            self.objs.extend(obj)
        else:
            self.objs.append(obj)

    def connectobjs(self):
        objs = self.objs
        objsnew = []
        for obj in objs:
            objsall = list(self.objs)
            obj.addobj(objsall.remove(obj))
            objsnew.append(obj)
        self.objs = objsnew

    def ani_init(self):
        objs = self.objs
        lims = self.lims
        self.fig =  plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(True, linestyle = '-', color = '0.75')
        self.ax.set_xlim(lims[0])
        self.ax.set_ylim(lims[1])
        self.dots = []
        self.traces = []

        for idx,obj in enumerate(objs):
            self.dot = self.ax.scatter([],[],s=obj.size,c=obj.color,label=obj.name)
            self.trace, = self.ax.plot([],[],c=obj.color)
            self.dots.append(self.dot)
            self.traces.append(self.trace)
        
        time = '0yr:0mo:0dy:0hr:0mn:0s'

        self.timetext = self.ax.text(0.05, 0.05, time, transform=self.ax.transAxes, fontsize=10, zorder=255,
            verticalalignment='top', fontdict={'family': "monospace", 'color':  'black', 'weight': 'normal', 'size': 10})
        self.timetext.set_bbox(dict(facecolor='#EEEEEE', alpha=0.75, edgecolor='#DBDBDB', boxstyle='round', pad=0.2))

        self.propagate()


    def animatestep(self,idx):
        idx += self.frame0
        for n,obj in enumerate(self.objs):
            self.traces[n].set_xdata(obj.posx[:idx+1])
            self.traces[n].set_ydata(obj.posy[:idx+1])
            self.dots[n].set_offsets([[obj.posx[idx],obj.posy[idx]]])
        text = self.gettime(idx*self.timestep)
        self.timetext.set_text(text)
        self.drawings = [self.timetext]
        for drawing in self.dots+self.traces:
            self.drawings.append(drawing) 
        return self.drawings

    def animate(self,save=False,extension='mp4',filename='nbody_vis'):
        self.ani_init()
        anim = animation.FuncAnimation(self.fig, self.animatestep, interval=self.interval,
                                frames = self.frames-1, blit=True)
        self.ax.set_aspect('equal')
        if save==True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=5000, codec='h264')
            anim.save('{}.{}'.format(filename,extension), writer=writer)
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()