from eos3dmpl.nbody import *
from eos3dmpl.eos_core import *
import datetime
import os

# Surpress any warnings
import warnings
warnings.filterwarnings("ignore")

os.system('cls') 

def option1():
    eos1 = Eos()
    eos1.view_init()
    eos1.orbit(a=26554,e=0.72,tp=4000,i=63.4,om=0,w=-90,unit_i='km',unit_o='km',body='Earth')
    eos1.kepler(vis='v-line',marker=True,curr_marker=True)
    eos1.view(x=8,y=8,bodysurface=True,bodyalpha=0.2,cust_txt='Molniya Orbit Visualization')

def option2():
    print()
    print('ABOUT:')
    print('The positions of the two spacecraft are those as of Mar 13, 2018.')
    print('Orbital elements have been retrieved from JPL\'s HORIZONS system.')
    eos1 = Eos()
    eos1.view_init()
    eos1.orbit(a=2.438805150792154E-05,e=5.662947956850604E-03,i=1.179742492581005E+02,om=1.033013426844067E+02,w=2.888728222352130E+02,th=3.387081555968664E+02,unit_i='AU',unit_o='km',body='Mars')
    eos1.kepler(vis='green',curr_color='green',marker=True,curr_marker=True,curr_label='Mars Reconaissance Orbiter')
    eos1.orbit(a=5.891726051080387E-05,e=5.735242533925626E-01,i=1.179742492581005E+02,om=2.811918452951060E+02,w=1.158117224248728E+02,th=1.047987681678214E+02,unit_i='AU',unit_o='km',body='Mars')
    eos1.kepler(vis='red',curr_color='red',marker=True,curr_marker=True,curr_label='Mars Express')
    eos1.view(x=8,y=8,bodysurface=True,bodyalpha=0.2,cust_txt='MRO & MEX Orbits', legend=True)

def option3(satname,legend):
    eos1 = Eos()
    #eos1.getTLE()
    eos1.view_init()
    eos1.parseTLE()
    eos1.satorbit(SATNAME=satname,cont=True, res=50, nums=False, info=True, marker=False, abbrev=False, vline=False, label='SATNAME', legend=legend, bg=False, orbit=True, orbitcolor=None, markercolor=None, SGP=False, curr_marker=True)
    eos1.view(x=8, y=8, bodysurface=True, bodyalpha=0.2, bodycolor='white', cust_txt='{} Orbit(s)'.format(satname), save=True)

def option3handler():
    try:
        satname = input('GIVE SATNAME(S) [STR/LIST]: ')
    except:
        raise ValueError('Invalid SATNAME')
    else:
        try:
            legend = str(input('SHOW LEGEND     [y/N]: '))
        except:
            legend = False
        else:
            if legend.upper() == 'Y':
                legend = True
            elif legend.upper() == 'N':
                legend = False
            else:
                legend = False
        option3(satname,legend)

def option4():
    eos1 = Eos()
    #eos1.getTLE()
    eos1.view_init()
    eos1.parseTLE()
    eos1.satorbit(SATNAME=['NAVSTAR'],cont=True, nocont=['DEB','R/B'], res=25, nums=False, info=True, marker=False, abbrev=False, vline=False, label='SATNAME', legend=False, bg=False, orbit=True, orbitcolor=None, markercolor=None, SGP=False, curr_marker=True)
    eos1.view(x=8, y=8, bodysurface=True, bodyalpha=0.2, bodycolor='white', cust_txt='NAVSTAR Orbits')

def option5():
    eos1 = Eos()
    #eos1.getTLE()
    eos1.view_init()
    eos1.parseTLE()
    eos1.satorbit(SATNAME=['DELFI'],cont=True, nocont=['DEB','R/B'], res=25, nums=False, info=True, marker=False, abbrev=False, vline=False, label='SATNAME', legend=True, bg=True, orbit=True, orbitcolor=None, markercolor=None, SGP=False, curr_marker=True)
    eos1.view(x=8, y=8, bodysurface=True, bodyalpha=0.2, bodycolor='white', cust_txt='Delfi Satellite Orbits')

def option6():
    if not eos.eos_core.sgp4_present:
        raise ValueError('SGP4 not installed, option unavailable')
    eos1 = Eos()
    #eos1.getTLE()
    eos1.view_init()
    eos1.parseTLE()
    eos1.satorbit(SATNAME=['N3XT'],cont=True, nocont=['DEB','R/B'], res=25, nums=False, info=True, marker=False, abbrev=False, vline=False, label='SATNAME', legend=True, bg=False, orbit=False, orbitcolor=None, markercolor=None, SGP=True, curr_marker=True)
    eos1.view(x=8, y=8, bodysurface=True, bodyalpha=0.2, bodycolor='white', cust_txt='Delfi N3XT Current Position')


def option7():
    eos1 = Eos()
    #eos1.getTLE()
    eos1.view_init()
    eos1.parseTLE()
    eos1.satorbit(SATNAME=['FALCON 9 R/B'],cont=True, res=50, nums=False, info=True, marker=False, abbrev=False, vline=False, label='SATNAME', legend=True, bg=False, orbit=True, orbitcolor=None, markercolor=None, SGP=False, curr_marker=True)
    eos1.view(x=8, y=8, el=90, az=0, bodysurface=True, bodyalpha=0.2, bodycolor='black', cust_txt='Falcon 9 Upper Stage Orbits')

def option8():
    eos1 = Eos()
    eos1.view_init()
    eos1.solarsystem(now=True,view=True,res=50)
    eos1.view(x=8, y=8, el=90, az=0, legend=True)

def option9():
    eos1 = Eos()
    eos1.view_init(center='Jupiter')
    eos1.solarsystem(now=True,view=True,res=50,unit='AU',sphere=True)
    eos1.orbit(q=1.35,Q=5.4,i=79.11,unit_i='AU',unit_o='AU',body='Sun')
    eos1.kepler(vis='yellow',res=200,curr_marker=False,curr_label='Ulysses')
    eos1.orbit(q=75600,Q=8.1e6,i=90,unit_i='km',unit_o='AU',body='Jupiter')
    eos1.kepler(vis='green',res=200,curr_marker=False,curr_label='Juno')
    eos1.view(x=8, y=8,cust_txt='Ulysses and Juno Orbits',legend=True)

def option10():
    eos1 = Eos()
    eos1.view_init()

    eos1.orbit(Q=0.4055e6,
                q=0.3633e6,
                tp=-395734.454737,
                i=0,
                om=0,
                w=0,
                unit_i='km',
                unit_o='km',
                body='Earth')
    eos1.kepler(vis='blue', marker=False, curr_color='blue', curr_marker=True, curr_background=False, curr_label='Moon Orbit')
    eos1.orbit(q=dict_R['Earth']/1000.+185.2,
                Q=0.3633e6,
                tp=395734.454737,
                i=0,
                om=0,
                w=180,
                unit_i='km',
                unit_o='km',
                body='Earth')
    TLO = eos1.kepler(vis='red', marker=False, data=True, curr_color='red', curr_marker=True, curr_background=False, curr_label='Lunar Transfer Orbit (LTO)')
    eos1.orbit(a=dict_R['Earth']/1000.+185.2,
                e=0,
                ta=0,
                i=0,
                om=0,
                w=0,
                unit_i='km',
                unit_o='km',
                body='Earth')
    eos1.kepler(vis='green', marker=False, curr_color='green', curr_marker=True, curr_background=False, curr_label='Earth Parking Orbit')
    eos1.view(x=8,y=8,el=90,az=0,bodysurface=True,bodycolor='white',legend=True,cust_txt='Lunar Rendezvous Profile')

def optionselect(option):
    if option == 1:
        option1()
    elif option == 2:
        option2()
    elif option == 3:
        option3handler()
    elif option == 4:
        option4()
    elif option == 5:
        option5()
    elif option == 6:
        option6()
    elif option == 7:
        option7()
    elif option == 8:
        option8()
    elif option == 9:
        option9()
    elif option == 10:
        option10()

def mainmenu():
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
    print('      Eos class showcase')
    print('    ++++++++++++++++++++++')
    print()
    print()
    print('ABOUT:')
    print('The Eos class allows for orbit calculation and visualization. This is not limited to using keplerian elements.')
    print('SPACETRACK satellite data may be used to track satellites currently in orbit, with the option of applying SGP4.')
    print('This demo serves to showcase some of the functions the Eos class has, demonstrating SPACETRACK capabilities,')
    print('orbit visualization including velocity colormapping and Solar System visualization.')
    print()
    print('OPTIONS:')
    print('+ General Orbit Visualization (orbit,kepler,view):')
    print('  (1) Molniya orbit with velocity visualization')
    print('  (2) Mars Reconaissance Orbiter & Mars Express orbits')
    print('+ Satellite Tracking (satorbit):')
    print('  (3) Search SATNAME & display orbit')
    print('  (4) NAVSTAR satellite orbits')
    print('  (5) Delfi satellite orbits')
    print('  (6) Current Delfi N3XT postion (requires SGP4)')
    print('  (7) Falcon 9 upper stage orbits')
    print('+ Solar System:')
    print('  (8) Current solar system celestial body location')
    print('  (9) Ulysses and Juno SC orbits in solar system')
    print('+ Earth-Moon-SC Orbits')
    print('  (10) Orbits used in lunar trajectory calculations')
    print('  (11) Quit')
    print()
    option = int(input('SELECT OPTION [1-11]: '))
    if option == 11:
        running = False
    elif option < 1 or option > 11:
        print('Invalid option')
    else:
        optionselect(option)
    

mainmenu()