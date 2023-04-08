#!/usr/bin/env python3
# coding=UTF-8
"""
Name: doscar_IO (ver 2.0)
Created on Fri Sep 18 10:57:09 2019
Developer: Ming-Wen Chang
E-mail: m.chang@tue.nl
"""
#import re
#import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('Agg')
#from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import simps
from scipy.interpolate import make_interp_spline
from scipy.signal import hilbert
from collections import OrderedDict

TDOS_channels =  OrderedDict([
                    (2, ['energy', 'tdos']),
                    (3, ['energy', 'tdos_up', 'tdos_down'])
                    ])      


PDOS_channels = OrderedDict([
                (4, ['energy','s', 'p', 'd']),  
                 
                (5, ['energy','s', 'p', 'd', 'f']),  
                
                (7, ['energy','s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']),
                
                (9, ['energy','s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down', 'f_up', 'f_down']),
                 
                (10, ['energy','s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z^2', 'd_xz', 'd_x^2-y^2']),  
                 
                (17, ['energy','s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z^2', 'd_xz', 'd_x^2-y^2', 
                     'f_y(3x^2-y^2)', 'f_xyz', 'f_yz^2', 'f_z^3', 'f_xz^2', 'f_z(x^2-y^2)', 'f_x(x^2-3y^2)']),
                
                (19, ['energy','s_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down', 
                    'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z^2_up', 'd_z^2_down', 'd_xz_up', 
                    'd_xz_down', 'd_x^2-y^2_up', 'd_x^2-y^2_down']), 
                
                (33, ['energy','s_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                    'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z^2_up', 'd_z^2_down', 'd_xz_up',
                    'd_xz_down', 'd_x^2-y^2_up', 'd_x^2-y^2_down', 
                    'f_y(3x^2-y^2)_up', 'f_y(3x^2-y^2)_down', 'f_xyz_up', 'f_xyz_down', 'f_yz^2_up', 'f_yz^2_down', 
                    'f_z^3_up', 'f_z^3_down', 'f_xz^2_up', 'f_xz^2_down', 'f_z(x^2-y^2)_up', 'f_z(x^2-y^2)_down', 
                    'f_x(x^2-3y^2)_up', 'f_x(x^2-3y^2)_down'])
                ])   

#Defalut style
lc  = '#212121' 
rc1 = {'axes.axisbelow': True,
       'axes.edgecolor': '#212121',
       'axes.facecolor': 'white',    
       'axes.labelcolor': '#212121',
       'axes.linewidth': 1.00,
       'figure.facecolor': 'white',
       'font.family': ['sans-serif'],
       'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif'],
       'axes.grid': True,
       'grid.color': '.8',
       'grid.linewidth':0.50,
       'grid.linestyle': '--',
       'grid.alpha':0.80,
       'axes.spines.left':True,
       'axes.spines.right':True,
       'axes.spines.top':True,
       'axes.spines.bottom':True,
       'legend.frameon': True,
       'legend.numpoints': 1,
       'legend.scatterpoints': 1,
       'lines.solid_capstyle': 'round',
       'text.color': '#212121',
       'xtick.color': '#212121',
       'xtick.direction': 'out',
       'xtick.major.size': 3.5,
       'xtick.minor.size': 3.5,
       'ytick.color': '#212121',
       'ytick.direction': 'out',
       'ytick.major.size': 3.5,
       'ytick.minor.size': 3.5}  


cmap =  ['blue', 'red', 'gold', 'salmon', 'lightcoral', 'lightskyblue', 
          'darkgreen', 'black', 'orange','powderblue','olivedrab', 'burlywood',
          'indianred', 'steelblue', 'lawngreen', 'y', 'hotpink', 'slategrey', 
          'yellowgreen','palegreen','sandybrown', 'tomato', 'darkviolet', 
          'lightgreen', 'tan','maroon']
fig1 = {'figsize':(6, 1.5),
        'figdpi': 150,
        'figstack':'row',
        'fontsize': 12,
        'linewidth':0.50,
        'linealpha': 1.00,
        'linestyle':'-',            
        'fill':True,
        'fillalpha': 0.50,
        'axvlinewidth':1.00,
        'axvlinealpha':0.90,
        'axvlinestyle':'-',
        'axvlinecolor':'#393E46',
        'axhlinewidth':1.00,
        'axhlinealpha':0.90,
        'axhlinestyle':'-',
        'axhlinecolor':'#393E46',
        'nxticks': 10,
        'nyticks':5,
        'xlabel':'$E - E_{f}\ (eV)$',
        'xlabelfontsize':12,          
        'ylabel':'pDOS',
        'ylabelfontsize':12,
        'tickwhich':'major',
        'ticklabelleft':True,
        'tickleft':True,
        'ticklabelbottom':True,
        'tickbottom':True,
        'tickdirection':'in',
        'tickfontsize':12,
        'showdbcline':True,
        'dbclinewidth':1.00,
        'dbclinealpha':0.50,
        'dbclinestyle':'--',
        'dbclinecolor':'#393E46',
        'showlegend':True,
        'legendboxposition':(1.00, 1.02),
        'legendboxframeon':True,
        'legendfontsize':10,
        'cmap':cmap}

class Doscar:
    
    def __init__(self, filename=None, efermi=None, undecomposed=True):  
        self.filename = filename
        self.undecomposed = undecomposed
        self.efermi = efermi
        self._efermi = 0.00
        self._energy = None
        self.data = None
        self.pdos = None
        self.analdbc = False
        
        if filename is not None:
            self.data = self.read_DOSCAR(self.filename)
            self._energy = self.data[0]['energy']
            
        if self.efermi is not None:
            self.set_efermi(self.efermi)
            
    def set_efermi(self, efermi):
        delta = efermi - self._efermi #shift value
        self._efermi  = efermi
        self._energy += delta 
        if self.data is not None:
            for df in self.data:
                df['energy'] = self._energy
                
        if self.pdos is not None:
            self.pdos['energy'] = self._energy
        
    def get_efermi(self):
        return self._efermi
            
    def read_DOSCAR(self, filename='DOSCAR'):   
        with open(filename, 'r' ) as txt: 
            natoms = int(txt.readline().split()[0]) #read line 1
            [txt.readline() for i in range(3)] #skip line2 to line4
            
            if 'LOBSTER' in txt.readline(): # read line5 
                lobster = True
            else:
                lobster = False
            
            data = []
            for i in range(0, natoms+1): #read line6 - end
                dos = []
                head = txt.readline().split()
                nedos = int(head[2])
                
                if not lobster:
                    efermi = float(head[3])
                else:
                    efermi = 0.00
            
                if len(head) == 5:#for vasp
                    factors = None
                else:#for lobster 
                    factors = head[7:]
    
                for j in range(nedos):
                    dos.append([float(value) for value in txt.readline().split()])
                dos = pd.DataFrame(data=dos)   
                    
                nn =  len(dos.columns)
                if i == 0:
                    if nn == 5: #spin-ploarized 
                        nn -= 2
                        self.ispin= True
                        dos = dos.iloc[:, 0:nn] #strip int_tdos_up and int_tdos_down 
                    else:#non spin ploarized 
                        nn -= 1
                        self.ispin= False
                        dos = dos.iloc[:, 0:nn] #strip int_tdos_up and int_tdos_down 
                        
                    dos.iloc[:, 0] -= efermi 
                    energy = dos.iloc[:, 0]
                    factors = TDOS_channels[nn]
                    
                else:
                    if factors is None: #vasp
                        factors = PDOS_channels[nn]
                    else: #lobster 
                        if nn - 1 > len(factors): #spin-ploarized 
                            ll = [ ]
                            for ft in factors:
                                ll.append(ft + '_up')
                                ll.append(ft + '_down')
                            factors = ll
                            factors.insert(0, 'energy')
                        else:
                            factors.insert(0, 'energy')
                    dos.iloc[:, 0] =energy
                dos.columns = factors  
                dos = self._invert_dos_values_of_spin_down(dos)
                data.append(dos)
                
            return data
                     
    def sum_atomic_pdos(self, dfobjs):
        energy = dfobjs[0]['energy']
        df = sum(dfobjs)
        if self.undecomposed :
            df = self._reduce_to_undecomposed(df)
        df['energy'] = energy
        return df
    
    def concat_atomic_pdos(self, dfobjs):
        energy = dfobjs[0]['energy']
        dfobjs = [df.iloc[:, 1:] for df in dfobjs]
        df = pd.concat(dfobjs, axis=1) 
        df.insert(loc=0, column='energy', value=energy)
        return df
                    
    def to_tdos(self, df=None, columns=None):
        if df is None:
            df = self.pdos
        if columns is None:
            if self.ispin:
                columns=['energy', 'tdos_up', 'tdos_down']
            else:
                columns=['energy', 'tdos']
         
        df2 = pd.DataFrame(data=0.00, index=range(len(df)), columns=columns)
        df2.iloc[:, 0] = df['energy']
        for idx in range(1, len(df.columns[1:]) + 1):
            if self.ispin:
                if idx%2 !=0:#spin up
                    df2.iloc[:, 1] += df.iloc[:,idx]
                else:#spin down
                    df2.iloc[:, 2] += df.iloc[:,idx]
            else:
                df2.iloc[:, 1] += df.iloc[:,idx]

        return df2
    
    def select_orbital(self, df, orb='d'):
        
        if len(orb) == 1: #i.e. s, p, d or f     
            #pattern = r'.*-?\d?%s_.*' %(orb)    
            pattern = r'.*-?\d?%s_.*' %(orb)    
            ids = [0]+ [idx for idx, ft in enumerate(df.columns) if re.match(pattern, ft)]
            df = df.iloc[:,ids]
        else:
            #pattern = re.sub('(_up|_down)','', orb) #strip "_up" and "_down" words
            pattern = '.*-?\d?' + re.sub('\^','\^', orb)
            ids = [0]+ [idx for idx, ft in enumerate(df.columns) if re.match(pattern, ft)]
            df = df.iloc[:,ids]
            pass
            
        return df  
    
    def get_pdos(self, kwargs, data=None):
        if data is None:
            data = self.data
            
        for idx, key in enumerate(kwargs.keys()):
            atom = key
            orb =  kwargs[key][0] 
            #sum the pdos of selected atoms 
            df = self.sum_atomic_pdos([data[idx] for idx in kwargs[key][1]]) 
            #Select target orbital and re-lable columns
            df = self.select_orbital(df, orb)
            df.columns = ['energy'] + [atom + '-' + ft for ft in df.columns[1:]]
            if idx == 0:
                df2 = df
            else:
                df2 = self.concat_atomic_pdos([df2, df])   
        self.pdos = df2
        return df2             
                
    def calculate_dbc(self, df=None, orbital='d', erange=None):
        if df is None:
            df = self.pdos
            
        df = self.select_orbital(df, orbital)
        epsilon = df['energy']  
        
        if erange is None:
            erange = (epsilon.iloc[0], self._efermi+4.00)
            #erange = (epsilon.iloc[0], epsilon.iloc[-1])
            
        if self.ispin: #spin polarize 
            rho1 = df.iloc[:,1].values
            rho2 = df.iloc[:,2].values
        else:#non-spin polarize 
            rho1 = df.iloc[:,1].values
            rho2 = df.iloc[:,1].values
          
        #Set intergration range for the calculation of nunber of electrons in d-band(ned)  
        mask = (epsilon >= erange[0]) & (epsilon <= self._efermi) 
        x = epsilon[mask] #integrating the DOS to the Fermil level
        y1 = rho1[mask] #density up
        y2 = rho2[mask] #density down
        ned_up = abs(simps(y1, x))  #d-band filling by up electrons 
        ned_down = abs(simps(y2, x)) #d-band filling by down electrons
        if self.ispin:
            self.ned = ned_up + ned_down
        else:
            self.ned = ned_up 
                                    
        #Set intergration range for the calculation of d-band center and width
        mask = (epsilon >= erange[0]) & (epsilon <= erange[1])
        x = epsilon[mask] #integrating the DOS to the Fermil level + n_eV
        y1 = rho1[mask] #density up
        y2 = rho2[mask] #density down
        
        dbc_up = simps(y1*x, x) / simps(y1, x)  #d-band center up 
        dbc_down = simps(y2*x, x) / simps(y2, x) #d-band center down 
        dbw_up =   np.sqrt(simps(y1*x**2, x) / simps(y1, x))  #d-band width up 
        dbw_down = np.sqrt(simps(y2*x**2, x) / simps(y2, x))  #d-band width down 
        self.dbc =  dbc_up, dbc_down
        self.dbw = dbw_up,  dbw_down
        
        #Do Hilbert transform to get d-band edge
        #cf. https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.115114
        
        h_up = hilbert(1*y1)
        h_down = hilbert(-1*y2)
        #envelope_up =  np.abs(h_up)
        #envelope_down =  np.abs(h_down)
        dbe_up = x[np.argmax(h_up)] #d-band edge up
        dbe_down = x[np.argmax(h_down)] #d-band edge down
        self.dbe =  dbe_up, dbe_down
 
        self.analdbc = True   
        return self.dbc     
    
    def split_doscar(self):
        for idx, df in enumerate(self.data):
            if self.undecomposed:
                df = self._reduce_to_undecomposed(df)
            filename='dos' + str(idx)
            self.save_df_to_txt(df, filename)
            
    def save_df_to_txt(self, df, filename='df.txt'):
        #df = self._invert_dos_values_of_spin_down(df)
        with open(filename, 'w') as txt:
            txt.write(df.to_string(index=False))
            
    def save_dos(self, df=None, filename='pdos'):
        if df is None:
            try:
                df = self.pdos
                self.save_df_to_txt(df, 'pdos')
            except:
                df = self.data[0]
                self.save_df_to_txt(df, 'tdos')
        else:
            self.save_df_to_txt(df, filename)
        
    def _reduce_to_undecomposed(self, df):
        factors = [ ]
        for ft1 in df.columns[1:]:
            fft1 = ft1.split('_')[0] + '_' + ft1.split('_')[-1] 
            if fft1 not in  factors:
                factors.append(fft1)
        factors.insert(0, 'energy')
                
        df2 = pd.DataFrame(data=0.00, index=range(len(df)), columns=factors)
        for ft2 in df2.columns[1:]:
            for ft1 in df.columns[1:]:
                name = ft1.split('_')[0] + '_' + ft1.split('_')[-1]
                if ft2 == name:
                    df2[ft2] += df[ft1]
        df2['energy'] = df['energy']
        return df2
    
    def _invert_dos_values_of_spin_down(self, df):
        for ft in df.columns:
            if 'down' in ft:
                df[ft] *= -1
        return df
    
    def spline_df(self, df):
        lb = df.columns[0]
        df2 = pd.DataFrame()
        for ft in df.columns[1:]:
             xnew, ynew  = smooth_line(df.iloc[:,0], df[ft])
             df2[lb] = xnew
             df2[ft] = ynew
        return df2
    
    
    def plot_pdos(self, df=None, filename='pdos.png', 
                  erange=None, dosrange=None,smooth=True,
                  figParams=None, rcParams=None):
        
        if df is None:
            df = self.pdos
                        
        if smooth:
            df = self.spline_df(df)
        
        if erange is None:
            emin, emax  = (self._efermi-6.00, self._efermi+4.00)
        else:
            emin, emax = erange[0], erange[1] 

 
        import matplotlib as mpl   
        if rcParams is None:
            rcParams = rc1 
        for key, value in rcParams.items(): 
            mpl.rcParams[key] = value
            
        if figParams is None:
            figParams = fig1
            
        #df = self._invert_dos_values_of_spin_down(df)
        mask = ( df.iloc[:,0].values  >= emin) & ( df.iloc[:,0].values  <= emax)  
        evalues =  df.iloc[:,0].values[mask] 
        dmax = 0.00

        coloridx = 0 
        fig, ax = plt.subplots(figsize=figParams['figsize'], dpi=figParams['figdpi'])
        for ft in df.columns[1:]:
 
            #Color setting
            if self.undecomposed:
                label = ft.split('_')[0]
            else:
                label = ft.rstrip('_up').rstrip('_down')

            
            #Values for plotting
            cc = figParams['cmap'][coloridx]
            dosvalues = df[ft].values[mask] 
            
            if figParams['figstack'] == 'column':
                x , y = dosvalues, evalues
                if 'down' not in ft:  
                    ax.plot(x, y, linewidth=figParams['linewidth'], alpha=figParams['linealpha'], label=label, color=cc)
                    coloridx -= 1
                else:
                    ax.plot(x, y, linewidth=figParams['linewidth'], alpha=figParams['linealpha'], label='', color=cc)
                    
                if figParams['fill']:
                    ax.fill_betweenx(y, x, -0.000,  where=y <= 0, interpolate=True, color=cc, alpha=figParams['fillalpha']) 
                    ax.fill_betweenx(y, x, +0.000,  where=y >= 0, interpolate=True, color=cc, alpha=figParams['fillalpha']*0.1)
                    
            else:
                x , y = evalues, dosvalues
                if 'down' not in ft:  
                    ax.plot(x, y, linewidth=figParams['linewidth'], alpha=figParams['linealpha'], label=label, color=cc)
                    coloridx -= 1
                else:
                    ax.plot(x, y, linewidth=figParams['linewidth'], alpha=figParams['linealpha'], label='', color=cc)
                    
                if figParams['fill']: 
                    ax.fill_between(x, y, -0.000,  where=x <= 0, interpolate=True, color=cc, alpha=figParams['fillalpha']) 
                    ax.fill_between(x, y, +0.000,  where=x >= 0, interpolate=True, color=cc, alpha=figParams['fillalpha']*0.1) 
           
            coloridx += 1 
            if np.max(np.abs(dosvalues)) > dmax:
                dmax = np.max(np.abs(dosvalues))

        if dosrange is None:
            if self.ispin:
                dmin, dmax = -1.10*dmax, 1.10*dmax
            else:
                dmin, dmax = -0.00, 1.10*dmax
        else:
            dmin, dmax = dosrange[0], dosrange[1]
            
        #Set range, lable, ticks    
        ax.set_ylabel(figParams['ylabel'], size=figParams['ylabelfontsize'])
        ax.set_xlabel(figParams['xlabel'], size=figParams['xlabelfontsize'])
        ax.xaxis.set_major_locator(plt.MaxNLocator(figParams['nxticks']))
        ax.yaxis.set_major_locator(plt.MaxNLocator(figParams['nyticks']))
        ax.tick_params(which=figParams['tickwhich'], 
                       labelleft=figParams['ticklabelleft'], 
                       left=figParams['tickleft'], 
                       labelbottom=figParams['ticklabelbottom'], 
                       bottom=figParams['tickbottom'], 
                       direction =figParams['tickdirection'], 
                       labelsize=figParams['tickfontsize']) 
            

        if figParams['figstack'] == 'column':
            ax.set_xlim([dmin, dmax])
            ax.set_ylim([emin, emax])
            ax.axvline(x=0, 
                       linewidth=figParams['axvlinewidth'], 
                       linestyle=figParams['axvlinestyle'], 
                       color=figParams['axvlinecolor'], 
                       alpha=figParams['axvlinealpha'])
            
            ax.axhline(y=self._efermi, 
                       linewidth=figParams['axhlinewidth'], 
                       linestyle=figParams['axhlinestyle'], 
                       color=figParams['axhlinecolor'], 
                       alpha=figParams['axhlinealpha'])
            
            if self.analdbc and figParams['showdbcline']:
                ax.axhline(y=max(self.dbc),
                           xmin=0.00,
                           xmax=1.00, 
                           linewidth=figParams['dbclinewidth'],
                           linestyle=figParams['dbclinestyle'],
                           color=figParams['dbclinecolor'],
                           alpha=figParams['dbclinealpha'])
        else:
            ax.set_xlim([emin, emax])
            ax.set_ylim([dmin, dmax])
            ax.axvline(x=self._efermi, 
                       linewidth=figParams['axvlinewidth'], 
                       linestyle=figParams['axvlinestyle'], 
                       color=figParams['axvlinecolor'], 
                       alpha=figParams['axvlinealpha'])
            
            ax.axhline(y=0, 
                       linewidth=figParams['axhlinewidth'], 
                       linestyle=figParams['axhlinestyle'], 
                       color=figParams['axhlinecolor'], 
                       alpha=figParams['axhlinealpha'])
            
            if self.analdbc and figParams['showdbcline']:
                ax.axvline(x=max(self.dbc), 
                           ymin=0.00, 
                           ymax=1.00, 
                           linewidth=figParams['dbclinewidth'],
                           linestyle=figParams['dbclinestyle'],
                           color=figParams['dbclinecolor'],
                           alpha=figParams['dbclinealpha'])
                
        if figParams['showlegend']: 
            legend =ax.legend(bbox_to_anchor=figParams['legendboxposition'],
                              #loc= 'upper left', 
                              prop={'size':figParams['legendfontsize']}, 
                              frameon=figParams['legendboxframeon']) 
            [i.set_linewidth(3) for i in legend.legendHandles]
        else:
            ax.legend().remove()

        fig.savefig(filename, bbox_inches="tight")    
    
    
                        
def smooth_line(x, y, ngrids=None):
    if ngrids is None:
        ngrids = len(x) * 100
    xnew = np.linspace(x.min(), x.max(), ngrids) 
    bsplobj = make_interp_spline(x, y, k=3) #BSpline object
    ynew = bsplobj(xnew) #smoothed data
    return xnew, ynew    


def spline_df(df, ngrids=None):
    if ngrids is None:
        ngrids = 10 * len(df)
    lb = df.columns[0]
    df2 = pd.DataFrame()
    for ft in df.columns[1:]:
         xnew, ynew  = smooth_line(df.iloc[:,0], df[ft], ngrids)
         df2[lb] = xnew
         df2[ft] = ynew
    return df2

debug = False
if debug:
    atoms = list(range(1, 46))
    #atoms = list([168])
    #kwargs = {'Pt':('d_x^2-y^2', atoms)}
    kwargs = {'Pd':('d', atoms)}
    a = Doscar('../DOSCAR', undecomposed=True)
    a.get_pdos(kwargs)
    a.pdos
    a.calculate_dbc(orbital='d')
    a.plot_pdos(erange=(-10,2), stack='row', size=(6, 1.5))
    a.save_dos()
    print(a.dbc)
    print(a.dbw)
    print(a.dbe)
    #print(a.nelectrons_in_pdos)

 
        
        
 

    

