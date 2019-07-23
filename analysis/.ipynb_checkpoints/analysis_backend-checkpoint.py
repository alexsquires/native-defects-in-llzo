############## import modules and functions ###################
import seaborn
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy.stats import linregress
from vasppy.calculation import *
import numpy as np
from math import log

############## import constants ###################
from scipy.constants import Boltzmann as kb
from scipy.constants import N_A as avogadro
from scipy.constants import atm as atm
from scipy.constants import physical_constants

############## define additional constants ###################
kb_e = physical_constants['Boltzmann constant in eV/K'][0]
j_to_ev = 1/physical_constants['electron volt-joule relationship'][0]
S_0 = 205 * j_to_ev / avogadro
Cp = (7/2)*kb_e
E_vbm = 3.6468                    # valence band maximum from reference caluclation

############## import DFT data #############
defects =  import_calculations_from_file('data/defects-la-zr.yaml') # need a final update
elements = import_calculations_from_file('data/elements.yaml')
interest = import_calculations_from_file('data/interest.yaml')

pbe_defects =  import_calculations_from_file('data/pbe_defects.yaml') # need a final update
pbe_elements = import_calculations_from_file('data/pbe_elements.yaml')

############## define classes ###################
class ChargeState:

    def __init__( self, charge, formation_energy, degeneracy ):
        self.charge = charge
        self.formation_energy = formation_energy
        self.degeneracy = degeneracy

class ChargeState:

    def __init__( self, charge, formation_energy, degeneracy ):
        self.charge = charge
        self.formation_energy = formation_energy
        self.degeneracy = degeneracy

class PotentialsSet:

    def __init__(self, Li_pot, O_pot, La_pot, Zr_pot, o_vac, li_vac, li_int, la_vac, zr_vac, la_zr, zr_la, zr_li, o_int, zr_li_tet):
        self.Li_pot = Li_pot
        self.O_pot = O_pot
        self.La_pot = La_pot
        self.Zr_pot = Zr_pot
        self.o_vac = o_vac
        self.li_vac = li_vac
        self.li_int = li_int
        self.la_vac = la_vac
        self.zr_vac = zr_vac
        self.la_zr = la_zr
        self.zr_la = zr_la
        self.zr_li = zr_li
        self.o_int = o_int
        self.zr_li_tet = zr_li_tet


class PotentialsSetResult:

    def __init__(self, Li_pot, O_pot, o_vac_0, o_vac_1, o_vac_2, li_vac_0, li_vac_1, li_int_0, li_int_1, la_vac_3, zr_vac_4, la_zr_0, la_zr_1, zr_la_0, zr_la_1, zr_li_0, zr_li_1, zr_li_2, zr_li_3, o_int_1, o_int_0, zr_li_tet_0, zr_li_tet_1, zr_li_tet_2, zr_li_tet_3):
        self.Li_pot = Li_pot
        self.O_pot = O_pot
        self.o_vac_0 = o_vac_0
        self.o_vac_1 = o_vac_1
        self.o_vac_2 = o_vac_2
        self.li_vac_0 = li_vac_0
        self.li_vac_1 = li_vac_1
        self.li_int_0 = li_int_0
        self.li_int_1 = li_int_1
        self.la_vac_3 = la_vac_3
        self.zr_vac_4 = zr_vac_4
        self.la_zr_0 = la_zr_0
        self.la_zr_1 = la_zr_1
        self.zr_la_0 = zr_la_0
        self.zr_la_1 = zr_la_1
        self.zr_li_0 = zr_li_0
        self.zr_li_1 = zr_li_1
        self.zr_li_2 = zr_li_2
        self.zr_li_3 = zr_li_3
        self.o_int_1 = o_int_1
        self.o_int_0 = o_int_0
        self.zr_li_tet_0 = zr_li_tet_0
        self.zr_li_tet_1 = zr_li_tet_1
        self.zr_li_tet_2 = zr_li_tet_2
        self.zr_li_tet_3 = zr_li_tet_3

class Defect:

    def __init__( self, label, charge_states, n_sites ):
        self.label = label
        self.charge_states = charge_states
        self.n_sites = n_sites

    @property
    def n_charge_states( self ):
        return len( self.charge_states )

class SCFermi:

    def __init__( self, defects, nelect, e_gap, temperature, spin_polarised=False ):
        self.defects = defects
        self.nelect = nelect
        self.e_gap = e_gap
        self.temperature = temperature
        self.spin_polarised = spin_polarised

    @property
    def n_defects( self ):
        return len( self.defects )

    def output( self ):

            with open('input-fermi.dat', 'w') as f:

                if self.spin_polarised:
                    f.write( '2' + '\n')
                else:
                    f.write( '1' + '\n' )
                f.write( str(self.nelect) + '\n' )
                f.write( str(self.e_gap) + '\n')
                f.write( str(self.temperature) + '\n')
                f.write( str(self.n_defects) + '\n' )
                for d in self.defects:
                    f.write( '{} {} {}'.format( d.label, d.n_charge_states, d.n_sites ) + '\n')
                    for c in d.charge_states:
                        f.write( '{} {} {}'.format( c.charge, c.formation_energy, c.degeneracy ) + '\n')


                f.close()

############## define functions ########################

def sc_fermi_vacancy_wrap( calc, lattice_site, stoich, limit, cor ):
    if re.search('\-\-\-\-', str(calc.title)) is not None:
        charge = -4
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (-4 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  -4, formation_energy, 4)

    elif re.search('\-\-\-', str(calc.title)) is not None:
        charge = -3
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (-3 * ( E_vbm )) + cor
        #print( calc.title,round(formation_energy,1), charge )
        return ChargeState(  -3, formation_energy, 4 )

    elif re.search('\+\+', str(calc.title)) is not None:
        charge = +2
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (2 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        #print(limit,cor)
        return ChargeState(  2, formation_energy, 4)

    elif re.search('\+', str(calc.title)) is not None:
        charge = 1
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  1, formation_energy, 2)

    elif re.search('\-', str(calc.title)) is not None:
        charge = -1
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (-1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  -1, formation_energy, 2)

    else:
        charge = 0
        formation_energy = ( calc.energy - stoich ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit))  + (0 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  0, formation_energy, 1)

def sc_fermi_interstitial_wrap( calc, lattice_site, stoich, limit, cor ):
    charge = 2
    if re.search('\+\+', str(calc.title)) is not None:
        formation_energy = ( calc.energy - stoich ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (2 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  2, formation_energy, 3)

    elif re.search('\+', str(calc.title)) is not None:
        charge = 1
        formation_energy = ( calc.energy - stoich ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  1, formation_energy, 2)

    elif re.search('\-', str(calc.title)) is not None:
        charge = -1
        formation_energy = ( calc.energy - stoich ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (-1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge)
        return ChargeState(  -1, formation_energy, 2)

    else:
        charge = 0
        formation_energy = ( calc.energy - stoich ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (0 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  0, formation_energy, 1)

def sc_fermi_sub_wrap( calc, non_native, native, stoich, non_native_limit, native_limit, cor ):
    if re.search('\+\+\+', str(calc.title)) is not None:
        charge = +3
        formation_energy = ( calc.energy - stoich ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (3 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), 3 )
        return ChargeState(  3, formation_energy, 2)

    elif re.search('\+\+', str(calc.title)) is not None:
        charge = +2
        formation_energy = ( calc.energy - stoich ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (2 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge)
        return ChargeState(  2, formation_energy, 2)

    elif re.search('\+', str(calc.title)) is not None:
        charge = 1
        formation_energy = ( calc.energy - stoich ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  1, formation_energy, 2)

    elif re.search('\-', str(calc.title)) is not None:
        charge = -1
        formation_energy = ( calc.energy - stoich ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (-1 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), -1 )
        return ChargeState(  -1, formation_energy, 2)

    else:
        charge = 0
        formation_energy = ( calc.energy - stoich ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (0 * ( E_vbm )) + cor
        #print( calc.title, round(formation_energy,1), charge )
        return ChargeState(  0, formation_energy, 1)


def dependance(P,T):
    chem_pot = 0.5 * ( (Cp * (T - 298))
                      - T * ( (S_0 + (Cp * np.log(T/298)) + (kb_e * np.log((1/P)) ) ) ))  #### This function gives dependance of mu_O(T,P)
    return chem_pot                                                                       #### BEWARE, CURRENTLY BACKWARDS

li_limit_1 = linregress((-4.7561,-0.0561),(-0.0242,-2.3742))  # li upper and lower bounds
li_limit_2 = linregress((-4.5561,-0.0561),(-0.3242,-2.5742))  # give me the upper and lower bounds of the li stabilty region

la_limit_1 = linregress((-4.7940,0),(0,-7.5181))  # la upper and lower bounds

zr_limit_1 = linregress((-4.7940,0),(0,-9.2133))  # zr upper and lower bounds
zr_limit_2 = linregress((-4.7940,0),(-0.2466,-10.0085))  # zr upper and lower bounds

def line(m,c,x): # straight line, duh
    y = m*x+c
    return y

def li_limits(o_limit):                  # from my mu_o value (typically as given by dependance fn), work out the corresponding Li limits
    m_1 = li_limit_1.slope
    c_1 = li_limit_1.intercept
    m_2 = li_limit_2.slope
    c_2 = li_limit_2.intercept
    li_limits_1 = line(m_1,c_1,o_limit)
    li_limits_2 = line(m_2,c_2,o_limit)
    return( li_limits_1, li_limits_2)

def la_limits(o_limit):                  # from my mu_o value (typically as given by dependance fn), work out the corresponding La limits
    m_1 = la_limit_1.slope
    c_1 = la_limit_1.intercept
    la_limits_1 = line(m_1,c_1,o_limit)
    return( la_limits_1)

def zr_limits(o_limit):                  # from my mu_o value (typically as given by dependance fn), work out the corresponding zr limits
    m_1 = zr_limit_1.slope
    c_1 = zr_limit_1.intercept
    m_2 = zr_limit_2.slope
    c_2 = zr_limit_2.intercept
    zr_limits_1 = line(m_1,c_1,o_limit)
    zr_limits_2 = line(m_2,c_2,o_limit)
    return( zr_limits_1, zr_limits_2)

def give_me_concs(P,T):                   #### given a pressure, and a temperature:
    o = dependance(P,T)                   ### give me my mu_O....
    li = li_limits(o)                     ### and the corresponding li limits....
    la = la_limits(o)
    zr = zr_limits(o)
    li_mid = (li[0] + li[1])/2            ### work out the mid-point of the limits, TAKE CARE, LI LIMTIS ARE INCORRECT AT EXTREME ENDS OF STAB REGION
    li_avg = make_calcs(o,li_mid,la,zr[0])         ###
    sc_fermi_run_avg = run_some_fermi(li_avg,T)
    d_avg = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_avg.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_avg.Li_pot,
         'vo2':sc_fermi_run_avg.o_vac_2,
         'vo1':sc_fermi_run_avg.o_vac_1,
         'vo0':sc_fermi_run_avg.o_vac_0,
         'vli0':sc_fermi_run_avg.li_vac_0,
         'vli1':sc_fermi_run_avg.li_vac_1,
         'ili1':sc_fermi_run_avg.li_int_1,
         'ili0':sc_fermi_run_avg.li_int_0,
         'zrli1':sc_fermi_run_avg.zr_li_1,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_avg = pd.DataFrame(data=d_avg, index=['Li_avg'])
    li_1 = make_calcs(o,li[0],la,zr[0])         ###
    sc_fermi_run_1 = run_some_fermi(li_1,T)
    d_1 = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_1.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_1.Li_pot,
         'vo2':sc_fermi_run_1.o_vac_2,
         'vo1':sc_fermi_run_1.o_vac_1,
         'vo0':sc_fermi_run_1.o_vac_0,
         'vli0':sc_fermi_run_1.li_vac_0,
         'vli1':sc_fermi_run_1.li_vac_1,
         'ili1':sc_fermi_run_1.li_int_1,
         'ili0':sc_fermi_run_1.li_int_0,
         'zrli1':sc_fermi_run_avg.zr_li_1,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_1 = pd.DataFrame(data=d_1, index=['Li_avg'])
    li_2 = make_calcs(o,li[1],la,zr[1])         ###
    sc_fermi_run_2 = run_some_fermi(li_2,T)
    d_2 = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_2.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_2.Li_pot,
         'vo2':sc_fermi_run_2.o_vac_2,
         'vo1':sc_fermi_run_2.o_vac_1,
         'vo0':sc_fermi_run_2.o_vac_0,
         'vli0':sc_fermi_run_2.li_vac_0,
         'vli1':sc_fermi_run_2.li_vac_1,
         'ili1':sc_fermi_run_2.li_int_1,
         'ili0':sc_fermi_run_2.li_int_0,
         'zrli1':sc_fermi_run_avg.zr_li_1,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_2 = pd.DataFrame(data=d_2, index=['Li_avg'])
    return df_avg, df_1, df_2
    #return df

from scipy.spatial import ConvexHull
def encircle(x,y, ax=None, **kw):
    if not ax: ax=plt.gca()
    p = np.c_[x,y]
    hull = ConvexHull(p)
    poly = plt.Polygon(p[hull.vertices,:], **kw)
    ax.add_patch(poly)

def give_me_concs_sometimes(P,T):                   #### given a pressure, and a temperature:
    o = dependance(P,T)                   ### give me my mu_O....
    li = li_limits(o)                     ### and the corresponding li limits....
    la = la_limits(o)
    zr = zr_limits(o)
    li_mid = (li[0] + li[1])/2
    zr_mid = (zr[0] + zr[1])/2            ### work out the mid-point of the limits, TAKE CARE, LI LIMTIS ARE INCORRECT AT EXTREME ENDS OF STAB REGION
    zr_avg = make_calcs(o,li_mid,la,zr_mid)         ###
    sc_fermi_run_avg = run_some_fermi(zr_avg,T)
    d_avg = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_avg.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_avg.Li_pot,
         'vo2':sc_fermi_run_avg.o_vac_2,
         'vo1':sc_fermi_run_avg.o_vac_1,
         'vo0':sc_fermi_run_avg.o_vac_0,
         'vli0':sc_fermi_run_avg.li_vac_0,
         'vli1':sc_fermi_run_avg.li_vac_1,
         'ili1':sc_fermi_run_avg.li_int_1,
         'ili0':sc_fermi_run_avg.li_int_0,
         'zrla0':sc_fermi_run_avg.zr_la_0,
         'zrla1':sc_fermi_run_avg.zr_la_1,
         'vla3':sc_fermi_run_avg.la_vac_3,
         'vzr4':sc_fermi_run_avg.zr_vac_4,
         'lazr0':sc_fermi_run_avg.la_zr_0,
         'lazr1':sc_fermi_run_avg.la_zr_1,
         'zrlioct0':sc_fermi_run_avg.zr_li_0,
         'zrlioct1':sc_fermi_run_avg.zr_li_1,
         'zrlioct2':sc_fermi_run_avg.zr_li_2,
         'zrlioct3':sc_fermi_run_avg.zr_li_3,
         'zrlitet0':sc_fermi_run_avg.zr_li_tet_0,
         'zrlitet1':sc_fermi_run_avg.zr_li_tet_1,
         'zrlitet2':sc_fermi_run_avg.zr_li_tet_2,
         'zrlitet3':sc_fermi_run_avg.zr_li_tet_3,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_avg = pd.DataFrame(data=d_avg, index=['Li_avg'])
    li_1 = make_calcs(o,li[1],la,zr[0])         ###
    sc_fermi_run_1 = run_some_fermi(li_1,T)
    d_1 = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_1.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_1.Li_pot,
         'vo2':sc_fermi_run_1.o_vac_2,
         'vo1':sc_fermi_run_1.o_vac_1,
         'vo0':sc_fermi_run_1.o_vac_0,
         'vli0':sc_fermi_run_1.li_vac_0,
         'vli1':sc_fermi_run_1.li_vac_1,
         'ili1':sc_fermi_run_1.li_int_1,
         'ili0':sc_fermi_run_1.li_int_0,
         'zrli1':sc_fermi_run_avg.zr_li_1,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_1 = pd.DataFrame(data=d_1, index=['Li_avg'])
    li_2 = make_calcs(o,li[0],la,zr[1])         ###
    sc_fermi_run_2 = run_some_fermi(li_2,T)
    d_2 = {'P / atm': P, 'T / K': T, 'mu_O':sc_fermi_run_2.O_pot, '$\mu_\mathrm{Li}$ / eV':sc_fermi_run_2.Li_pot,
         'vo2':sc_fermi_run_2.o_vac_2,
         'vo1':sc_fermi_run_2.o_vac_1,
         'vo0':sc_fermi_run_2.o_vac_0,
         'vli0':sc_fermi_run_2.li_vac_0,
         'vli1':sc_fermi_run_2.li_vac_1,
         'ili1':sc_fermi_run_2.li_int_1,
         'ili0':sc_fermi_run_2.li_int_0,
         'zrli1':sc_fermi_run_avg.zr_li_1,
         'io1':sc_fermi_run_avg.o_int_1,
         'io0':sc_fermi_run_avg.o_int_0}
    df_2 = pd.DataFrame(data=d_2, index=['Li_avg'])
    return df_avg, df_1, df_2
    #return df

    ###### The seemingly random numbers added to the end of the 'sc_fermi' calls are corrections from sxdefectalign

def make_calcs(o,li,la,zr):

    defects =  import_calculations_from_file('data/defects-la-zr.yaml') # need a final update
    elements = import_calculations_from_file('data/elements.yaml')
    interest = import_calculations_from_file('data/interest.yaml')

    O = sc_fermi_vacancy_wrap( defects['O_vac'],elements['O'],interest['LLZO'].energy,o,0)
    O1= sc_fermi_vacancy_wrap( defects['O_vac+'],elements['O'],interest['LLZO'].energy,o,0.0408347)
    O2 = sc_fermi_vacancy_wrap( defects['O_vac++'],elements['O'],interest['LLZO'].energy,o,0.163339)
    Li = sc_fermi_vacancy_wrap( defects['Li_vac'],elements['Li'],interest['LLZO'].energy,li,0)
    Li1 = sc_fermi_vacancy_wrap( defects['Li_vac-'],elements['Li'],interest['LLZO'].energy,li,0.0408347)
    La = sc_fermi_vacancy_wrap( defects['La_vac---'],elements['La'],interest['LLZO'].energy,la,0.367513)
    Zr = sc_fermi_vacancy_wrap( defects['Zr_vac----'],elements['Zr'],interest['LLZO'].energy,zr,0.653356)
    Lii = sc_fermi_interstitial_wrap( defects['Li_int'],elements['Li'],interest['LLZO'].energy,li,0)
    Lii1 = sc_fermi_interstitial_wrap( defects['Li_int+'],elements['Li'],interest['LLZO'].energy,li,0.0408347)
    La_Zr = sc_fermi_sub_wrap( defects['La_Zr'],elements['La'],elements['Zr'],interest['LLZO'].energy,la,zr,0)
    Zr_La = sc_fermi_sub_wrap( defects['Zr_La'],elements['Zr'],elements['La'],interest['LLZO'].energy,zr,la,0.0408347)
    La_Zr1 = sc_fermi_sub_wrap( defects['La_Zr-'],elements['La'],elements['Zr'],interest['LLZO'].energy,la,zr,0.0408347)
    Zr_La1 = sc_fermi_sub_wrap( defects['Zr_La+'],elements['Zr'],elements['La'],interest['LLZO'].energy,zr,la,0.0408347)
    Zr_Li3 = sc_fermi_sub_wrap( defects['Zr_Li+++'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li,0.367513)
    Zr_Li2 = sc_fermi_sub_wrap( defects['Zr_Li++'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0.163339)
    Zr_Li1 = sc_fermi_sub_wrap( defects['Zr_Li+'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0.0408347)
    Zr_Li = sc_fermi_sub_wrap( defects['Zr_Li'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0)
    Zr_Li_tet3 = sc_fermi_sub_wrap( defects['Zr_Li_tet+++'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li,0.367513)
    Zr_Li_tet2 = sc_fermi_sub_wrap( defects['Zr_Li_tet++'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0.163339)
    Zr_Li_tet1 = sc_fermi_sub_wrap( defects['Zr_Li_tet+'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0.0408347)
    Zr_Li_tet = sc_fermi_sub_wrap( defects['Zr_Li_tet'],elements['Zr'],elements['Li'],interest['LLZO'].energy,zr,li, 0)
    Oi = sc_fermi_interstitial_wrap( defects['O_i'],elements['O'],interest['LLZO'].energy,o,0)
    Oi1 = sc_fermi_interstitial_wrap( defects['O_int-'],elements['O'],interest['LLZO'].energy,o,0)
    o_vac = Defect( 'V_O', [O,O1,O2], n_sites=3 )
    li_vac = Defect( 'V_Li', [Li,Li1], n_sites=3 )
    li_int = Defect( 'Li_i', [Lii,Lii1], n_sites=1 )
    la_vac = Defect( 'V_La', [La], n_sites=2 )
    zr_vac = Defect( 'V_Zr', [Zr], n_sites=1 )
    la_zr = Defect('La_Zr', [La_Zr,La_Zr1], n_sites=1 )
    zr_la = Defect('Zr_La', [Zr_La,Zr_La1], n_sites=2 )
    zr_li = Defect('Zr_Li', [Zr_Li,Zr_Li1,Zr_Li2,Zr_Li3], n_sites=2 )
    o_int = Defect('O_i', [Oi,Oi1], n_sites=1 )
    zr_li_tet = Defect('Zr_Li_tet', [Zr_Li_tet,Zr_Li_tet1,Zr_Li_tet2,Zr_Li_tet3], n_sites=2 )
    calcs = PotentialsSet( li, o, la, zr, o_vac, li_vac, li_int, la_vac, zr_vac, la_zr, zr_la, zr_li, o_int, zr_li_tet)
    return calcs

def off_stoich(conc):
    per_cubic_angstrom = conc / 1e+24
    per_unit_cell = per_cubic_angstrom * (13.003 * 13.003 * 12.498)
    return per_unit_cell

def run_some_fermi(for_analysis,T):
    scf = SCFermi( [ for_analysis.o_vac, for_analysis.li_vac, for_analysis.li_int, for_analysis.la_vac, for_analysis.zr_vac, for_analysis.la_zr, for_analysis.zr_la, for_analysis.zr_li, for_analysis.o_int, for_analysis.zr_li_tet], nelect=544, e_gap=5.9034, temperature=T, spin_polarised=True )
    scf.output()
    with open("blah.txt", 'w') as f:
        sp = subprocess.run(["./sc-fermi"],stdout=f)
        text_file = open(("blah.txt") , "r")
        lines =  text_file.readlines()
        zeros = []
        ones = []
        twos = []
        neg_ones = []
        fours = []
        threes =[]
        neg_threes =[]
        neg_fours =[]
        for line in lines:
                if re.search(r"\s0.0\s", line):
                    zero = (line)
                    zeros.append(zero)
                if re.search(r"\s1.0\s", line):
                    one = (line)
                    ones.append(one)
                if re.search(r"\s2.0\s", line):
                    two = (line)
                    twos.append(two)
                if re.search(r"\s3.0\s", line):
                    three = (line)
                    threes.append(three)
                if re.search(r"\s-3.0\s", line):
                    neg_three = (line)
                    neg_threes.append(neg_three)
                if re.search(r"\s-4.0\s", line):
                    neg_four = (line)
                    neg_fours.append(neg_four)
                if re.search(r"\s-1.0\s", line):
                    neg_one = (line)
                    neg_ones.append(neg_one)
        m = re.search(r"[^\s]+E[^\s]+", zeros[0])
        o0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[1])
        liv0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[2])
        lii0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[3])
        la_zr0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[4])
        zr_la0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[5])
        zr_li0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[6])
        o_0_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", ones[0])
        o1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", ones[1])
        lii1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", ones[2])
        zr_la1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", twos[0])
        o2_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", neg_ones[0])
        liv1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", neg_ones[1])
        la_zr1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", neg_ones[2])
        o_1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", neg_threes[0])
        la_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", neg_fours[0])
        zr_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", threes[0])
        zr_li3_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", twos[1])
        zr_li2_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", ones[3])
        zr_li1_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", threes[1])
        zr_li3_tet_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", zeros[7])
        zr_li0_tet_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", ones[4])
        zr_li1_tet_concs = float(m.group(0))
        m = re.search(r"[^\s]+E[^\s]+", twos[2])
        zr_li2_tet_concs = float(m.group(0))
        result = PotentialsSetResult(for_analysis.Li_pot, for_analysis.O_pot, o0_concs, o1_concs, o2_concs, liv0_concs, liv1_concs, lii0_concs, lii1_concs, zr_concs, la_concs, la_zr0_concs, la_zr1_concs, zr_la0_concs, zr_la1_concs, zr_li0_concs, zr_li1_concs, zr_li2_concs, zr_li3_concs, o_1_concs, o_0_concs, zr_li0_tet_concs, zr_li1_tet_concs, zr_li2_tet_concs, zr_li3_tet_concs)
        return(result)
