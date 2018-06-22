from eos3dmpl.nbody import *
from eos3dmpl.eos_core import *
import datetime
import os

AU = 149597870700

os.system('cls') 

def optionselect(option):
    if option < 1 or option > 6:
        print('Invalid option')
    elif option == 1 or option == 2:
        frame0 = 0
        interval = 1
        time = datetime.datetime.now()
        if option == 1:
            frames = 1000
            timestep = 1*24*3600
            lims = [[-2*AU,2*AU],[-2*AU,2*AU]]
        else:
            frames = 5000
            timestep = 10*24*3600
            lims = [[-35*AU,35*AU],[-35*AU,35*AU]]

        eos1 = Eos()

        eos1.solarsystem(time='now',view=False,res=1)

        Mercury_dat = eos1.ssdata[0]['SPEC']
        Venus_dat = eos1.ssdata[1]['SPEC']
        Earth_dat = eos1.ssdata[2]['SPEC']
        Mars_dat = eos1.ssdata[3]['SPEC']
        Jupiter_dat = eos1.ssdata[4]['SPEC']
        Saturn_dat = eos1.ssdata[5]['SPEC']
        Uranus_dat = eos1.ssdata[6]['SPEC']
        Neptune_dat = eos1.ssdata[7]['SPEC']
        Pluto_dat = eos1.ssdata[8]['SPEC']

        Sun = Body([0.,0.,0.],[0.,0.,0.],dict_R['Sun'],dict_M['Sun'],'Sun',dict_c['Sun'])
        Mercury = Body([Mercury_dat['X'][0]*1e3,Mercury_dat['Y'][0]*1e3,Mercury_dat['Z'][0]*1e3],[Mercury_dat['VX'][0]*1e3,Mercury_dat['VY'][0]*1e3,Mercury_dat['VZ'][0]*1e3],dict_R['Mercury'],dict_M['Mercury'],'Mercury',dict_c['Mercury'])
        Venus = Body([Venus_dat['X'][0]*1e3,Venus_dat['Y'][0]*1e3,Venus_dat['Z'][0]*1e3],[Venus_dat['VX'][0]*1e3,Venus_dat['VY'][0]*1e3,Venus_dat['VZ'][0]*1e3],dict_R['Venus'],dict_M['Venus'],'Venus',dict_c['Venus'])
        Earth = Body([Earth_dat['X'][0]*1e3,Earth_dat['Y'][0]*1e3,Earth_dat['Z'][0]*1e3],[Earth_dat['VX'][0]*1e3,Earth_dat['VY'][0]*1e3,Earth_dat['VZ'][0]*1e3],dict_R['Earth'],dict_M['Earth'],'Earth',dict_c['Earth'])
        Mars = Body([Mars_dat['X'][0]*1e3,Mars_dat['Y'][0]*1e3,Mars_dat['Z'][0]*1e3],[Mars_dat['VX'][0]*1e3,Mars_dat['VY'][0]*1e3,Mars_dat['VZ'][0]*1e3],dict_R['Mars'],dict_M['Mars'],'Mars',dict_c['Mars'])
        Jupiter = Body([Jupiter_dat['X'][0]*1e3,Jupiter_dat['Y'][0]*1e3,Jupiter_dat['Z'][0]*1e3],[Jupiter_dat['VX'][0]*1e3,Jupiter_dat['VY'][0]*1e3,Jupiter_dat['VZ'][0]*1e3],dict_R['Jupiter'],dict_M['Jupiter'],'Jupiter',dict_c['Jupiter'])
        Saturn = Body([Saturn_dat['X'][0]*1e3,Saturn_dat['Y'][0]*1e3,Saturn_dat['Z'][0]*1e3],[Saturn_dat['VX'][0]*1e3,Saturn_dat['VY'][0]*1e3,Saturn_dat['VZ'][0]*1e3],dict_R['Saturn'],dict_M['Saturn'],'Saturn',dict_c['Saturn'])
        Uranus = Body([Uranus_dat['X'][0]*1e3,Uranus_dat['Y'][0]*1e3,Uranus_dat['Z'][0]*1e3],[Uranus_dat['VX'][0]*1e3,Uranus_dat['VY'][0]*1e3,Uranus_dat['VZ'][0]*1e3],dict_R['Uranus'],dict_M['Uranus'],'Uranus',dict_c['Uranus'])
        Neptune = Body([Neptune_dat['X'][0]*1e3,Neptune_dat['Y'][0]*1e3,Neptune_dat['Z'][0]*1e3],[Neptune_dat['VX'][0]*1e3,Neptune_dat['VY'][0]*1e3,Neptune_dat['VZ'][0]*1e3],dict_R['Neptune'],dict_M['Neptune'],'Neptune',dict_c['Neptune'])
        Pluto = Body([Pluto_dat['X'][0]*1e3,Pluto_dat['Y'][0]*1e3,Pluto_dat['Z'][0]*1e3],[Pluto_dat['VX'][0]*1e3,Pluto_dat['VY'][0]*1e3,Pluto_dat['VZ'][0]*1e3],dict_R['Pluto'],dict_M['Pluto'],'Pluto',dict_c['Pluto'])

        objs = [Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto]

        nbdy = Nbody(timestep=timestep,frames=frames,frame0=frame0,interval=interval,lims=lims,time=time)
        nbdy.addobj(objs)
        nbdy.connectobjs()
        nbdy.animate(save=False,filename='solarsystem')
    elif option == 3 or option == 4:
        interval = 1
        if option == 3:
            timestep = 60
            frames = 15000
            frame0 = 0
            lims = [[-1e8,5e8],[-1e8,5e8]]
        else:
            timestep = 60
            frames = 1000
            frame0 = 5*24*60
            #lims = [[3.55e8,4.05e8],[-.05e8,.45e8]]
            lims = [[3.55e8,4.05e8],[.05e8,.55e8]]
        
        Earth = Body([0,0.,0.],[0,0,0],dict_R['Earth'],dict_M['Earth'],'Earth',dict_c['Earth'])
        Moon = Body([182780.060273e3,-359114.564224e3,0.],[0.852896e3,0.462594e3,0.],dict_R['Moon'],dict_M['Moon'],'Moon',dict_c['Moon'])
        #SC = Body([-6563.3e3,0.,0.],[0.,-(10.92267e3-2.7),0.],10,4e3,'SC','#45f442')
        SC = Body([-6563.3e3,0.,0.],[0.,-(10.92267e3),0.],10,4e3,'SC','#45f442')
        # GRAVITY ASSIST
        SC.addimpulse([0.0e3,.7e3,0],5*24*3600-4*3600+0*60)
        # ORBIT INJECTION 1
        SC.addimpulse([+0.3e3,-0.5e3,0],5*24*3600+8*3600+0*60)

        objs = [Earth,Moon,SC]

        nbdy = Nbody(timestep=timestep,frames=frames,frame0=frame0,interval=interval,lims=lims)
        nbdy.addobj(objs)
        nbdy.connectobjs()
        nbdy.animate(save=False,filename='tli')
    elif option == 5:
        interval = 1
        timestep = 3*3600
        frames = 5000
        frame0 = 0
        lims = [[-5e8,5e8],[-5e8,10e8]]

        Earth = Body([0,0.,0.],[0,0,0],dict_R['Earth'],dict_M['Earth'],'Earth',dict_c['Earth'])
        Moon = Body([0.4055e9,0.,0.],[0,.970e3,0.],dict_R['Moon'],dict_M['Moon'],'Moon',dict_c['Moon'])
        
        objs = [Earth,Moon]

        nbdy = Nbody(timestep=timestep,frames=frames,frame0=frame0,interval=interval,lims=lims)
        nbdy.addobj(objs)
        nbdy.connectobjs()
        nbdy.animate(save=False,filename='em-system')

def mainloop():
    running = True
    print("    __________  __________ ____        __  _______  __ ")
    print("   / ____/ __ \/ ___/__  // __ \      /  |/  / __ \/ / ")
    print("  / __/ / / / /\__ \ /_ </ / / /____ / /|_/ / /_/ / /  ")
    print(" / /___/ /_/ /___/ /__/ / /_/ //___// /  / / ____/ /___")
    print("/_____/\____//____/____/_____/     /_/  /_/_/   /_____/")
    print()
    print('    Extendable Orbit System 3D - Matplotlib version')
    print()
    print()
    print('    ++++++++++++++++++++++')
    print('     Nbody class showcase')
    print('    ++++++++++++++++++++++')
    print()
    print()
    print('ABOUT:')
    print('The Nbody class allows for n-body simulations of all scales, including those of planets and satellites.')
    print('This demo serves to showcase some of the functions the Nbody class has, including applications to real-life.')
    print('Along with a Solar System simulation, a simulation of an Apollo 11-like lunar rendezvous trajectory is provided.')
    print()
    print('OPTIONS:')
    print('+ Solar System:')
    print('  (1) Solar System n-body simulation (inner planets)')
    print('  (2) Solar System n-body simulation (overall)')
    print('+ Lunar Rendezvous (least energy):')
    print('  (3) Translunar Injection & Lunar Orbit Insertion (overall)')
    print('  (4) Translunar Injection & Lunar Orbit Insertion (close-up)')
    print('+ Earth-Moon:')
    print('  (5) Earth-Moon system n-body simulation (overall)')
    print('  (6) Quit')
    print()
    option = int(input('SELECT OPTION [1-6]: '))
    if option == 6:
        running = False
    else:
        optionselect(option)
    

mainloop()