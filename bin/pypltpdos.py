#!/usr/bin/env python3
# coding=UTF-8
"""
Name: pypltpdos
Created on Fri Sep 18 10:57:09 2019
Developer: Ming-Wen Chang
E-mail: m.chang@tue.nl
"""

import re
import sys
import vasp_modules.doscar_io as dio
from collections import OrderedDict

class State:
    #DOSCAR parameters
    dosfile = 'DOSCAR' 
    system = ''
    efermi = None
    undecomposed = True
    splitdos = False
    savepdos = True
    
    #For atom obj
    atomobj = False
    kwargs = OrderedDict() # {}
    obj = []
    
    #For mole obj
    moleobj = False
    kwargs_mole = OrderedDict()
    
    #DBC 
    dbc = False
    dbcorbital= 'd'
    dbcintegrange = None
    
    #Plot parameters 
    erange = None
    dosrange = None
    smooth = True
    plotpdos = True
 
    #
    atomobj = False
    moleobj = False
    nxobj = 0
    
    cmap =  ['blue', 'red', 'gold', 'salmon', 'lightcoral', 'lightskyblue', 
             'darkgreen', 'black','orange','powderblue','olivedrab', 'burlywood',
             'indianred', 'steelblue', 'lawngreen', 'yellow', 'hotpink', 'slategrey', 
              'yellowgreen','palegreen', 'sandybrown', 'tomato', 'darkviolet', 
              'lightgreen', 'tan','maroon']
        
    #Defalut style
    figParams = {
            'figsize':(6, 1.5),
            'figdpi': 150,
            'figstack':'row',
            'independent':True,
            
            
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
  
            'xlabel':'E - ${E_f}$ (eV)',
            'xlabelfontsize':12,          
            'ylabel':'PDOS',
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
    
    rcParams ={'axes.axisbelow': True,
               'axes.edgecolor': '#212121',
               'axes.facecolor': 'white',    
               'axes.labelcolor': '#212121',
               'axes.linewidth': 1.50,
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
               'ytick.minor.size': 3.5} # mpl.rcParams


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass 
    return False


            
def analyze_atom_object(string):
    import re
    #string = re.sub('type:|orbital:|atomlist:|\s', '', string)
    split_string = re.split(',',string)
    
    atomtype = split_string[0]
    atomorbital = split_string[1].lower()
    if atomorbital == 'all' or atomorbital == 'none':
        atomorbital = None
    
    atomlist = []
    for term in split_string[2:]:
        if '-' in term:
            start = int(term.split('-')[0])
            end = int(term.split('-')[-1])
            [atomlist.append(_) for _ in range(start, end+1)]
        else:
            if term.isdigit():
                term = int(term)
                atomlist.append(term)
            else: #color code
                #print (state.nxobj)
                state.cmap.insert(state.nxobj, term)
                
    #remove duplicates in  atomlist       
    atomlist = list(dict.fromkeys(atomlist))
    atomtype, atomorbital, atomlist
    return atomtype, atomorbital, atomlist
    
     
def read_input(filename='ini.dos'):
    with open(filename, 'r') as txt:
        for line in txt:
            #Skip lines start with: '#', '!', '\s', '\n'
            skip = any([line.startswith(_) for _ in ['#', '!', '\s', '\n']])
            if not skip:
                line = re.sub('(?<!\'|\")#.*', '', line)
                line = re.sub('\s|\'|\"', '', line)
                line = re.split('=', line)
                keyword = line[0].upper()
                
                try:
                    arguments = line[1]
                except:
                    continue
                print(keyword, ':', arguments)
                
                if keyword == 'SYSTEM':
                    if arguments != '':
                        state.jobname = arguments + '-' 
                        
                elif keyword == 'FILENAME':
                    state.dosfile = arguments
                                        
                elif keyword == 'EFERMI':
                    state.efermi = float(arguments)
                    
                elif keyword == 'ATOM_OBJECT':
                    atomtype, atomorbital, atomlist = analyze_atom_object(arguments)
                    state.kwargs[atomtype] = (atomorbital, atomlist)
                    state.atomobj = True
                    state.nxobj +=1
                    state.obj.append(analyze_atom_object(arguments))
                    
                elif keyword == 'MOLE_OBJECT':
                    arguments = re.split('[:\+]', arguments)
                    moletype = arguments[0].strip()
                    kwargs = OrderedDict()
                    state.tmp = []
                    for idx in range(1, len(arguments)):
                        atomtype, atomorbital, atomlist = analyze_atom_object(arguments[idx])
                        kwargs[atomtype] = (atomorbital, atomlist)
                        state.tmp.append((atomtype,(atomorbital, atomlist)))
                        state.kwargs_mole[moletype]=kwargs
                     
                    import copy
                    state.kwargs2 = copy.copy(state.kwargs)
                    state.kwargs2[moletype] =('mole', state.tmp)
                    state.moleobj = True
                    state.nxobj +=1
                    
                elif keyword == 'UNDECOMPOSED': 
                    if arguments.capitalize() == 'True':
                        state.undecomposed = True
                    else:
                        state.undecomposed = False  
  
                elif  keyword == 'SPLITDOS':
                    if arguments.capitalize() == 'True':
                        state.splitdos = True  
                    else:
                        state.splitdos = False 
                        
                elif  keyword == 'SAVEPDOS':
                    if arguments.capitalize() == 'True':
                        state.savepdos = True
                    else:
                        state.savepdos = False

                elif keyword == 'ANALDBC':
                    if arguments.capitalize() == 'True':
                        state.dbc = True 
                    else:
                        state.dbc = False
                        
                elif keyword == 'DBCORBITAL':
                    state.dbcorbital = arguments.lower()

                elif keyword == 'DBCINTEGRANGE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    if arguments[0].lower()[0] == 'd': #default
                        state.dbcintegrange  = None
                    else:
                        state.dbcintegrange = (float(arguments[0]), float(arguments[1]))
                        
                elif  keyword == 'PLOTPDOS':
                    if arguments.capitalize() == 'True':
                         state.plotpdos = True
                    else:
                        state.plotpdos = False
                        
                elif keyword == 'PLOTERANGE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    if arguments[0].lower()[0] == 'd': #default
                        state.erange =  None
                    else:
                        state.erange = (float(arguments[0]), float(arguments[1]))
                        
                elif keyword == 'PLOTDOSRANGE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    if arguments[0].lower()[0] == 'd': #default
                        state.dosrange =  None
                    else:
                        state.dosrange = (float(arguments[0]), float(arguments[1]))
                        
                elif keyword == 'SMOOTH':
                    if arguments.capitalize() == 'True':
                         state.smooth = True
                    else:
                        state.smooth = False
                        
                elif  keyword == 'FIGURESIZE' or keyword == 'FIGSIZE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    state.figParams['figsize'] = (float(arguments[0]), float(arguments[1])) 
                    
                elif keyword == 'FIGUREDPI' or keyword == 'FIGDPI' :
                    state.figParams['figdpi'] = int(arguments)
                        
                elif keyword == 'FILLAREA' or keyword == 'FILL': 
                    if arguments.capitalize() == 'True':
                        state.figParams['fill'] = True
                    else:
                        state.figParams['fill'] = False
        
                elif keyword == 'LINEWIDTH':
                    state.figParams['linewidth'] = float(arguments)
                elif keyword == 'LINEALPHA':
                    state.figParams['linealpha'] = float(arguments)    
                
                elif keyword == 'FILLALPHA':
                    state.figParams['fillalpha'] = float(arguments)
                    
                elif keyword == 'FONTSIZE':
                    fontsize = float(arguments)
                    for key in state.figParams.keys():
                        if 'fontsize' in key:
                            state.figParams[key] = fontsize
                    
                elif keyword == 'CMAP':
                    cmap = [term.strip(',').strip('[').strip(']') for term in arguments.split(',')]
                    state.figParams['cmap'] = cmap  +  state.figParams['cmap']
                    
                elif '[' and ']' in keyword: 
                    xparams = keyword.split('[')[0]
                    key = keyword.split('[')[1].strip(']').lower()
                    values = arguments.split(',')
                    
                 
                    if len(values) == 1: 
                        value = values[0]
                        if is_number(value):
                                value = float(value)
                        elif value == 'True':
                            value = True
                        elif value == 'False':
                            value = False        
                    else:
                        value = [float(i) for i in values]
                        pass
                            
                    if xparams == 'RCPARAMS':
                        state.rcParams[key] = value
                    else:
                        state.figParams[key] = value
 
                else:
                    raise ValueError('KEYWORD: %s is not supported' %(keyword))                   
                    
                    
                
                """
                elif keyword == 'STACK' or keyword == 'FIGSTACK' :
                    state.figParams['figstack'] = arguments.lower()
                    
                elif keyword == 'RCPARAMS' or  keyword == 'FIGPARAMS': 
                    arguments = arguments.strip(')').strip('(').split(',')
                    for argument in arguments:
                        key = argument.split(':')[0]
                        value = argument.split(':')[1]
                        if is_number(value):
                            value = float(value)
                        else:
                            if value == 'True':
                                value = True
                            elif value == 'False':
                                value = False
                        if keyword == 'RCPARAMS':
                            state.rcParams[key] = value
                        else:
                            state.figParams[key] = value

                elif  keyword == 'NEAT':
                    if arguments.capitalize() == 'True':
                        state.figParams['axhlinewidth']=0 
                        state.figParams['axvlinewidth']=0 
                        state.figParams['xlabel']=''
                        state.figParams['ylabel']=''
                        state.figParams['ticklabelleft']=False
                        state.figParams['tickleft']=False
                        state.figParams['ticklabelbottom']=False
                        state.figParams['tickbottom']=False
                        state.figParams['showlegend']=False
                        state.rcParams['axes.grid']=False
                        state.rcParams['axes.spines.left']=False
                        state.rcParams['axes.spines.right']=False
                        state.rcParams['axes.spines.top']=False
                        state.rcParams['axes.spines.bottom']=False
                        
                elif  keyword == 'ONLY_BOTTOM':
                    if arguments.capitalize() == 'True':
                        state.figParams['fill']=False
                        state.figParams['linewidth']=0 
                        state.figParams['axhlinewidth']=0 
                        state.figParams['axvlinewidth']=0 
                        state.figParams['xlabel']=''
                        state.figParams['ylabel']=''
                        state.figParams['ticklabelleft']=False
                        state.figParams['tickleft']=False
                        state.figParams['ticklabelbottom']=True
                        state.figParams['tickbottom']=True
                        state.figParams['showlegend']=False
                        state.rcParams['axes.grid']=False
                        state.rcParams['axes.spines.left']=False
                        state.rcParams['axes.spines.right']=False
                        state.rcParams['axes.spines.top']=False
                        state.rcParams['axes.spines.bottom']=True

                elif  keyword == 'ONLY_LEFT':
                    if arguments.capitalize() == 'True':
                        state.figParams['fill']=False
                        state.figParams['linewidth']=0 
                        state.figParams['axhlinewidth']=0 
                        state.figParams['axvlinewidth']=0 
                        state.figParams['xlabel']=''
                        state.figParams['ylabel']=''
                        state.figParams['ticklabelleft']=True
                        state.figParams['tickleft']=True
                        state.figParams['ticklabelbottom']=False
                        state.figParams['tickbottom']=False
                        state.figParams['showlegend']=False
                        state.rcParams['axes.grid']=False
                        state.rcParams['axes.spines.left']=True
                        state.rcParams['axes.spines.right']=False
                        state.rcParams['axes.spines.top']=False
                        state.rcParams['axes.spines.bottom']=False
                """



            
# Main Function   
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        inp = sys.argv[1]
    else:
        inp = 'pdos.ini'
    
    
    print("o=======================================================o")
    print("            Reading DOSCAR and pdos.ini files            ")
    print("o=======================================================o")
    
    state = State()
    read_input(filename=inp)

    doscar = dio.Doscar(filename=state.dosfile, efermi=state.efermi, 
                        undecomposed=state.undecomposed,)
    
    if state.splitdos:
        print("o=======================================================o")
        print("                 Splitting pDOS                          ")
        print("o=======================================================o")
        
        doscar.split_doscar()
        print("DOSCAR file has been splitted")    
    
    if state.atomobj:
        df_atomobj = doscar.get_pdos(state.kwargs)
   
    if state.moleobj:
        dfobjs = []
        for key in state.kwargs_mole.keys():
            kwargs = state.kwargs_mole[key]
            _pdos = doscar.get_pdos(kwargs)
            
            if doscar.ispin:
                columns = ['energy', key + '_up', key + '_down']
            else:
                columns = ['energy', key]
                
            _pdos = doscar.to_tdos(_pdos, columns)
            dfobjs.append(_pdos)
            
        df_moleobj = dfobjs[0]   
        for _pdos in dfobjs[1:]:
            df_moleobj = doscar.concat_atomic_pdos([df_moleobj, _pdos]) 
    
    if state.atomobj and state.moleobj:
        doscar.pdos = doscar.concat_atomic_pdos([df_atomobj, df_moleobj]) 
    elif state.atomobj and not state.moleobj:
        doscar.pdos = df_atomobj
    elif not state.atomobj and state.moleobj:
        doscar.pdos = df_moleobj
        

    if state.dbc:
        print("o=======================================================o")
        print("                 Analyzing %s-band center                " %(state.dbcorbital))
        print("o=======================================================o")
        doscar.calculate_dbc(orbital=state.dbcorbital, erange=state.dbcintegrange) 
        dbc = doscar.dbc
        dbw = doscar.dbw
        dbe = doscar.dbe
        ned = doscar.ned
        
        print("%s-band center:   %.2f (eV), %.2f (eV)" %(state.dbcorbital, dbc[0], dbc[1]))
        print("%s-band width:    %.2f (eV), %.2f (eV)" %(state.dbcorbital, dbw[0], dbw[1]))
        print("%s-band edge:    %.2f (eV), %.2f (eV)" %(state.dbcorbital, dbe[0], dbe[1]))
        print("%s-band filling:  %.2f (no.)" %(state.dbcorbital, ned))
        print("\n")
        
        #fmt0 = '%-15s %15.2f eV'
        txt = open(state.jobname+'dbc.dat', 'w')
        txt.write("%s-band center:   %.2f (eV), %.2f (eV)\n" %(state.dbcorbital, dbc[0], dbc[1]))
        txt.write("%s-band width:    %.2f (eV), %.2f (eV)\n" %(state.dbcorbital, dbw[0], dbw[1]))
        txt.write("%s-band edge:    %.2f (eV), %.2f (eV)\n" %(state.dbcorbital, dbe[0], dbe[1]))
        txt.write("%s-band filling:  %.2f (no.)\n" %(state.dbcorbital, ned))
        txt.close()
            
    if state.savepdos:
        doscar.save_dos(doscar.pdos, state.jobname+'pdos.dat')
        print("pDOS data have been saved to %s file" %(state.jobname+'pdos.dat'))
        
    if state.plotpdos:
        print("o=======================================================o")
        print("                 Plotting pDOS                          ")
        print("o=======================================================o")

        
        doscar.plot_pdos(filename=state.jobname+'pdos.png', 
                         erange=state.erange, 
                         dosrange=state.dosrange,
                         smooth=state.smooth, 
                         figParams=state.figParams,
                         rcParams=state.rcParams)
        print('pDOS figure has been saved to %s' %(state.jobname+'pdos.png'))


        

         

    
    
