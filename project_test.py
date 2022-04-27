#!/usr/bin/env python3

######################################
#    AraSim simple event reader      #
#        Dennis H. calderon          #
#    calderon-madera.1@osu.edu       #
######################################

#######################################################
"""
=======================
##project_test.py##
======================
Author: Dennis H. Calderon
Email: calderon-madera.1@osu.edu
Date: Nomveber 02, 2021
Modified: March 24, 2022
=======================
Descripiton: 
This PYTHON script takes two sets of AraSim output .root files. For each set, it makes a cut for triggered events, pulls variables, and makes histograms comparing the two. 

This scrpit was make for a comparison of antennas for the ARA Bicone (vpol) and an evolved antenna using GENETIS (vpol). This current verion is comparing Direct & Refracted/Reflected Events using variables (theta_rec, rec_ang, reflect_ang) for each simulation run.
=======================
Usage:
python project.py <source> [options] <source_2>
<source_1> is where the ROOT file from your AraSim output
<source_2> is path where the other ROOT file to compare
<source_3> is path where the other ROOT file to compare
<source_4> is path where the other ROOT file to compare
<source_5> is path where the other ROOT file to compare
<source_6> is path where the other ROOT file to compare.
=======================
Options:
[-s2, -s3, -s4, -s5, -s6]  tells program that you are putting in anoter source of simulation files.
=======================
example:
python all_vars.py ../output_files/AraOut.Bicone.run{0..9}.root -s2 ../output_files/AraOut.GENETIS.run{0..9}.root
=======================
"""

#######################################################
import timeit
start = timeit.default_timer()
#######################################################
print("\n")
print('\033[1;37m#\033[0;0m'*50)
print("Now running \033[1;4;5;31mproject_test.py\033[0;0m!")
print('\033[1;37m#\033[0;0m'*50)
print('\n')
##########################################
print("\033[1;37mPlease wait patiently...\033[0;0m")
print('Importing libraries...')

##########################################
#System libraries
#import sys
import argparse
#import csv
#import types
#import os
import warnings
warnings.filterwarnings("ignore")
print('...')

#PyRoot libraries
import ROOT
#from ROOT import TCanvas, TGraph
#from ROOT import gROOT
from ROOT import gInterpreter, gSystem
#from ROOT import TChain, TSelector, TTree
from ROOT import TChain
print('...')

#Python libraries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
print('...')
##########################################

#####
#AraSim specific headers needed
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Position.h"')#"/users/PAS0654/dcalderon/AraSim/Position.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Report.h"')#"/users/PAS0654/dcalderon/AraSim/Report.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Detector.h"')#"/users/PAS0654/dcalderon/AraSim/Detector.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Settings.h"')#"/users/PAS0654/dcalderon/AraSim/Settings.h"')

gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/ara_build/lib/libAra.so')#'/users/PAS0654/dcalderon/AraSim/libAra.so') 

##########################################
# We want to give an output file as an input. This checks that we have a fle to read
parser = argparse.ArgumentParser(
        description='Read AraSim file and produce some plots. Can also compare two AraSim files.')
parser.add_argument("source_1", help = "Path to the AraSim file you want to use.", nargs='+')
parser.add_argument("--source_2", "-s2", help = "Path to another AraSim file you want to comprare to.", nargs='+')
parser.add_argument("--source_3", "-s3", help = "Path to another AraSim file you want to comprare to.", nargs='+')
parser.add_argument("--source_4", "-s4", help = "Path to another AraSim file you want to comprare to.", nargs='+')
parser.add_argument("--source_5", "-s5", help = "Path to another AraSim file you want to comprare to.", nargs='+')
parser.add_argument("--source_6", "-s6", help = "Path to another AraSim file you want to comprare to.", nargs='+')

g = parser.parse_args()

#print(g)
#print('#'*28)

# source_name = g.source[0].split('.')[1]
# source_2_name = g.source_2[0].split('.')[1]
# 
# print(source_name)
# print(source_2_name)
# print('#'*28)
##########################################
'''
can put this inside as well
'''

#Making a dictionary of the parsed arguments
source_dict = g.__dict__
#Deleting empty arguments from dictionary
source_dict = {k:v for k, v in source_dict.items() if v != None}
print('#'*28)
print(source_dict)
print('#'*28)
print('\n')


########################
##Variables needed
########################
energy = np.power(10,18)
earth_depth = 6359632.4
core_x = 10000.0
core_y = 10000.0
#stations[i].strings[j].antennas[k].GetX() << " : " <<

# Alocal_x = np.array(Aposnu_x) - 10000.0
# Alocal_y = np.array(Aposnu_y) - 10000.0

# Blocal_x = np.array(Bposnu_x) - 10000.0
# Blocal_y = np.array(Bposnu_y) - 10000.0

# Adisplace = np.sqrt(Alocal_x**2 + Alocal_y**2)
# Bdisplace = np.sqrt(Blocal_x**2 + Blocal_y**2)

# Adepth = 6359632.4 - np.array(Aposnu_z)
# Bdepth = 6359632.4 - np.array(Bposnu_z)



##################################
###Loop over Evenets
##################################
##########################################

print('#'*28)
print("Now lets do the loop")
print("Please wait patiently...")
print('...')
print('\n')

data_dict = {}
for i in range(len(source_dict.keys())):
        
        #General info for each simulation set
        print('#'*50)
        
        #setting trees
        var_dict = {}
        #list of all variable names
        var = ['trigg', 'weight', 'posnu_x', 'posnu_y', 'posnu_z',
               'rec_ang_0', 'theta_rec_0', 'reflect_ang_0',
               'dist_0', 'arrival_time_0', 'reflection_0', 
               'l_att_0', 'view_ang_0', 'launch_ang_0',
               'rec_ang_1', 'theta_rec_1', 'reflect_ang_1',
               'dist_1', 'arrival_time_1', 'reflection_1', 
               'l_att_1', 'view_ang_1', 'launch_ang_1',
               'current', 'flavor', 'elast',
               'nnu_theta', 'nnu_phi', 'ShowerEnergy',
               'depth', 'distance']
        
        #loop for making dictionary of variables and empty list
        for x in var:
                var_dict['{0}'.format(x)] = []
                #print('yes')
                # print(var_dict)

        SimTree = [] #sets SimTree and makes empty list
        SimTree = TChain("AraTree2") #Which tree I will be chaining
        for line in list(source_dict.values())[i]: #for every filename in my list
                SimTree.AddFile(line)
        reportPtr = ROOT.Report()#report pointer
        eventPtr = ROOT.Event()#event pointe
        #detectorPtr = ROOT.Detector()
        #can also add more pointers if needed
        #print(reportPtr)
        #print(SimTree)
        SimTree.SetBranchAddress("report", ROOT.AddressOf(reportPtr))
        SimTree.SetBranchAddress("event", ROOT.AddressOf(eventPtr))
        #SimTree.SetBranchAddress("detector", ROOT.AddressOf(detectorPtr))
        
        #basic info of data
        totalEvents = SimTree.GetEntries()
        key = []
        key =  list(source_dict)[i]
        
        print('\033[1;37m{0}\033[0;0m'.format(key))
        print('Total Events: {0}'.format(totalEvents))
        print('#'*50)

        #Now we loop over all the events 
        for j in range(totalEvents):
                #print(j)
                SimTree.GetEntry(j)
                
                #Selecting only triggered events and a weight between 0 and 1
                if (reportPtr.stations[0].Global_Pass > 0) and (eventPtr.Nu_Interaction[0].weight >= 0 and eventPtr.Nu_Interaction[0].weight <= 1):
                        #print(j)
                        #print(key)
                        trigg = j
                        var_dict['trigg'].append(j)
                        #If value is seen in both antennas (Top Vpol and Bot Vpol) then we take an average of two
                        #var_dict['trigg'].append(j)
                        try:                                                                 
                                #interaction position in ice
                                posnu_x = eventPtr.Nu_Interaction[0].posnu.GetX()
                                posnu_y = eventPtr.Nu_Interaction[0].posnu.GetY()
                                posnu_z = eventPtr.Nu_Interaction[0].posnu.GetZ()
                                
                                #Getting angle of received signal in antenna
                                #Direct solutioins
                                rec_ang_0 = ((reportPtr.stations[0].strings[1].antennas[0].rec_ang[0] + 
                                             reportPtr.stations[0].strings[1].antennas[2].rec_ang[0])/2.0)
                                reflect_ang_0 = ((reportPtr.stations[0].strings[1].antennas[0].reflect_ang[0] +
                                                 reportPtr.stations[0].strings[1].antennas[2].reflect_ang[0])/2.0)
                                theta_rec_0 = ((reportPtr.stations[0].strings[1].antennas[0].theta_rec[0] +
                                               reportPtr.stations[0].strings[1].antennas[2].theta_rec[0])/2.0)
                                
                                dist_0 = reportPtr.stations[0].strings[1].antennas[0].Dist[0]
                                arrival_time_0 = reportPtr.stations[0].strings[1].antennas[0].arrival_time[0] 
                                reflection_0 = reportPtr.stations[0].strings[1].antennas[0].reflection[0]
                                l_att_0 = reportPtr.stations[0].strings[1].antennas[0].L_att[0]
                                
                                view_ang_0 = reportPtr.stations[0].strings[1].antennas[0].view_ang[0]
                                launch_ang_0 = reportPtr.stations[0].strings[1].antennas[0].launch_ang[0]

                                #Refracted/Reflected solutions
                                rec_ang_1 = ((reportPtr.stations[0].strings[1].antennas[0].rec_ang[1] +
                                             reportPtr.stations[0].strings[1].antennas[2].rec_ang[1])/2.0)
                                reflect_ang_1 = ((reportPtr.stations[0].strings[1].antennas[0].reflect_ang[1] +
                                                 reportPtr.stations[0].strings[1].antennas[2].reflect_ang[1])/2.0)
                                theta_rec_1 = ((reportPtr.stations[0].strings[1].antennas[0].theta_rec[1] +
                                               reportPtr.stations[0].strings[1].antennas[2].theta_rec[1])/2.0)
                                
                                dist_1 = reportPtr.stations[0].strings[1].antennas[0].Dist[1]
                                arrival_time_1 = reportPtr.stations[0].strings[1].antennas[0].arrival_time[1] 
                                reflection_1 = reportPtr.stations[0].strings[1].antennas[0].reflection[1]
                                l_att_1 = reportPtr.stations[0].strings[1].antennas[0].L_att[1]
                                
                                view_ang_1 = reportPtr.stations[0].strings[1].antennas[0].view_ang[1]
                                launch_ang_1 = reportPtr.stations[0].strings[1].antennas[0].launch_ang[1]
                                
                                #incomeing neutrino info
                                nnu_theta = eventPtr.Nu_Interaction[0].nnu.Theta()
                                nnu_phi = eventPtr.Nu_Interaction[0].nnu.Phi()
                                
                                current = eventPtr.Nu_Interaction[0].currentint
                                flavor = eventPtr.nuflavorint
                                elast = eventPtr.Nu_Interaction[0].elast_y
                                
                                #weight
                                weight = eventPtr.Nu_Interaction[0].weight
                                                
                                if current == 1 and flavor == 1:
                                        ShowerEnergy = energy                                        
                                else:
                                        ShowerEnergy = energy * elast
                                
                                depth = posnu_z - earth_depth
                                distance =  ((posnu_x - core_x)**2 + (posnu_y - core_y)**2 )**(0.5)
                                #detectorPtr.stations[0].strings[1].antennas[0].GetX()
                                

                                all_var = [trigg, weight, posnu_x, posnu_y, posnu_z,
                                       rec_ang_0, theta_rec_0, reflect_ang_0,
                                       dist_0, arrival_time_0, reflection_0, 
                                       l_att_0, view_ang_0, launch_ang_0,
                                       rec_ang_1, theta_rec_1, reflect_ang_1,
                                       dist_1, arrival_time_1, reflection_1, 
                                       l_att_1, view_ang_1, launch_ang_1,
                                       current, flavor, elast,
                                       nnu_theta, nnu_phi, ShowerEnergy,
                                       depth, distance]
                               
                                for k in range(1,len(all_var)):
                                        var_dict['{0}'.format(var[k])].append(all_var[k])
                                print(j)
                                        
                        except IndexError:
                                
                                #Both antennas didn't see a signal, so we try with index 0 (Bot Vpol)
                                try: 
                                        
                                        #interaction position in ice
                                        posnu_x = eventPtr.Nu_Interaction[0].posnu.GetX()
                                        posnu_y = eventPtr.Nu_Interaction[0].posnu.GetY()
                                        posnu_z = eventPtr.Nu_Interaction[0].posnu.GetZ()
                                        
                                        #angles seen by antenna
                                        rec_ang_0 = reportPtr.stations[0].strings[1].antennas[0].rec_ang[0]
                                        theta_rec_0 = reportPtr.stations[0].strings[1].antennas[0].theta_rec[0]
                                        reflect_ang_0 = reportPtr.stations[0].strings[1].antennas[0].reflect_ang[0]
                                        
                                        dist_0 = reportPtr.stations[0].strings[1].antennas[0].Dist[0]
                                        arrival_time_0 = reportPtr.stations[0].strings[1].antennas[0].arrival_time[0] 
                                        reflection_0 = reportPtr.stations[0].strings[1].antennas[0].reflection[0]
                                        l_att_0 = reportPtr.stations[0].strings[1].antennas[0].L_att[0]
                                        
                                        view_ang_0 = reportPtr.stations[0].strings[1].antennas[0].view_ang[0]
                                        launch_ang_0 = reportPtr.stations[0].strings[1].antennas[0].launch_ang[0]
                                       
                                        rec_ang_1 = reportPtr.stations[0].strings[1].antennas[0].rec_ang[1]
                                        theta_rec_1 = reportPtr.stations[0].strings[1].antennas[0].theta_rec[1]
                                        reflect_ang_1 = reportPtr.stations[0].strings[1].antennas[0].reflect_ang[1]
                                        
                                        #other info 
                                        
                                        dist_1 = reportPtr.stations[0].strings[1].antennas[0].Dist[1]
                                        arrival_time_1 = reportPtr.stations[0].strings[1].antennas[0].arrival_time[1] 
                                        reflection_1 = reportPtr.stations[0].strings[1].antennas[0].reflection[1]
                                        l_att_1 = reportPtr.stations[0].strings[1].antennas[0].L_att[1]
                                        
                                        view_ang_1 = reportPtr.stations[0].strings[1].antennas[0].view_ang[1]
                                        launch_ang_1 = reportPtr.stations[0].strings[1].antennas[0].launch_ang[1]       
                                        
                                        #incomeing neutrino info
                                        nnu_theta = eventPtr.Nu_Interaction[0].nnu.Theta()
                                        nnu_phi = eventPtr.Nu_Interaction[0].nnu.Phi()
                                        
                                        current = eventPtr.Nu_Interaction[0].currentint
                                        flavor = eventPtr.nuflavorint
                                        elast = eventPtr.Nu_Interaction[0].elast_y
                                        
                                        #weight
                                        weight = eventPtr.Nu_Interaction[0].weight
                                                
                                        if current == 1 and flavor == 1:
                                                ShowerEnergy = energy                                        
                                        else:
                                                ShowerEnergy = energy * elast
                                                
                                        depth = posnu_z - earth_depth
                                        distance = ((posnu_x - core_x)**2 + (posnu_y - core_y)**2 )**(0.5)
                                                                                        
                                        all_var = [trigg, weight, posnu_x, posnu_y, posnu_z,
                                                   rec_ang_0, theta_rec_0, reflect_ang_0,
                                                   dist_0, arrival_time_0, reflection_0, 
                                                   l_att_0, view_ang_0, launch_ang_0,
                                                   rec_ang_1, theta_rec_1, reflect_ang_1,
                                                   dist_1, arrival_time_1, reflection_1, 
                                                   l_att_1, view_ang_1, launch_ang_1,
                                                   current, flavor, elast,
                                                   nnu_theta, nnu_phi, ShowerEnergy,
                                                   depth, distance]
                                                
                                        for k in range(1,len(all_var)):
                                                var_dict['{0}'.format(var[k])].append(all_var[k])
                                                
                                        print(str(j)+" only has Bot Vpol signal")

                                except IndexError:
                                        try: #Have this here because not always that both antenna see
                                                
                                                #interaction position in ice
                                                posnu_x = eventPtr.Nu_Interaction[0].posnu.GetX()
                                                posnu_y = eventPtr.Nu_Interaction[0].posnu.GetY()
                                                posnu_z = eventPtr.Nu_Interaction[0].posnu.GetZ()
                                                
                                                #angles seen by antenna
                                                rec_ang_0 = reportPtr.stations[0].strings[1].antennas[2].rec_ang[0]
                                                theta_rec_0 = reportPtr.stations[0].strings[1].antennas[2].theta_rec[0]
                                                reflect_ang_0 = reportPtr.stations[0].strings[1].antennas[2].reflect_ang[0]
                                                
                                                rec_ang_1 = reportPtr.stations[0].strings[1].antennas[2].rec_ang[1]
                                                theta_rec_1 = reportPtr.stations[0].strings[1].antennas[2].theta_rec[1]
                                                reflect_ang_1 = reportPtr.stations[0].strings[1].antennas[2].reflect_ang[1]
                                                
                                                #other info 
                                                dist_0 = reportPtr.stations[0].strings[1].antennas[2].Dist[0]
                                                arrival_time_0 = reportPtr.stations[0].strings[1].antennas[2].arrival_time[0] 
                                                reflection_0 = reportPtr.stations[0].strings[1].antennas[2].reflection[0]
                                                l_att_0 = reportPtr.stations[0].strings[1].antennas[2].L_att[0]
                                                
                                                view_ang_0 = reportPtr.stations[0].strings[1].antennas[2].view_ang[0]
                                                launch_ang_0 = reportPtr.stations[0].strings[1].antennas[2].launch_ang[0]
                                                
                                                dist_1 = reportPtr.stations[0].strings[1].antennas[2].Dist[1]
                                                arrival_time_1 = reportPtr.stations[0].strings[1].antennas[2].arrival_time[1] 
                                                reflection_1 = reportPtr.stations[0].strings[1].antennas[2].reflection[1]
                                                l_att_1 = reportPtr.stations[0].strings[1].antennas[2].L_att[1]
                                                
                                                view_ang_1 = reportPtr.stations[0].strings[1].antennas[2].view_ang[1]
                                                launch_ang_1 = reportPtr.stations[0].strings[1].antennas[2].launch_ang[1]       
                                                
                                                #incomeing neutrino info
                                                nnu_theta = eventPtr.Nu_Interaction[0].nnu.Theta()
                                                nnu_phi = eventPtr.Nu_Interaction[0].nnu.Phi()
                                                
                                                current = eventPtr.Nu_Interaction[0].currentint
                                                flavor = eventPtr.nuflavorint
                                                elast = eventPtr.Nu_Interaction[0].elast_y
                                                
                                                #weight
                                                weight = eventPtr.Nu_Interaction[0].weight
                                                
                                                
                                                if current == 1 and flavor == 1:
                                                        ShowerEnergy = energy                                        
                                                else:
                                                        ShowerEnergy = energy * elast
                                                        
                                                depth = posnu_z - earth_depth
                                                distance = ((posnu_x - core_x)**2 + (posnu_y - core_y)**2 )**(0.5)
                                                
                                                all_var = [trigg, weight, posnu_x, posnu_y, posnu_z,
                                                           rec_ang_0, theta_rec_0, reflect_ang_0,
                                                           dist_0, arrival_time_0, reflection_0, 
                                                           l_att_0, view_ang_0, launch_ang_0,
                                                           rec_ang_1, theta_rec_1, reflect_ang_1,
                                                           dist_1, arrival_time_1, reflection_1, 
                                                           l_att_1, view_ang_1, launch_ang_1,
                                                           current, flavor, elast,
                                                           nnu_theta, nnu_phi, ShowerEnergy,
                                                           depth, distance]
                                                        
                                                for k in range(1,len(all_var)):
                                                        var_dict['{0}'.format(var[k])].append(all_var[k])
                                                                
                                                
                                                print(str(j)+" only has Top Vpol signal")                                                             
                                        except IndexError:
                                                print("Event "+str(j)+" has no signal in either Top or Bot Vpol")
                                                #exit()
                                                continue
                                                
        
        #end of loop                                                    
        data_dict['{0}'.format(list(source_dict.keys())[i])] = var_dict
        print("#"*28)
        print('\n')

        
#print(data_dict.keys())
print('\n')
print("We have now looped over alll events and selected only triggered events")
print("Now we can let the fun begin...")
print('#'*50)
print('\n')

#print(data_dict['source_1']['distance'])
#exit()
#######################################
###Plots
#######################################
print('#'*50)
print("Now lets make some pots!")
print('#'*50)

source_names = list(data_dict.keys())

##
w = 2.0
binsize = np.linspace(-1.0, 1.0, 41)
bindepth = 20
bindistance = np.linspace(0,4000, 21)

bin_cos = np.linspace(-1.0, 1.0, 41)
bin_dist = np.linspace(0,4000, 41)
binsize = np.linspace(-1.0, 1.0, 41)
bindepth = 20
bindistance = np.linspace(0,4000, 41)

'''
make histogram plot
for cos(rec_ang_0) and cos(rec_ang_1)
'''


##Setting up legends 
colors = ['r','b','g','c','m','y']


custom_lines_style = [Line2D([0], [0], color='k', ls='-'),
                      Line2D([0], [0], color='k', ls='--')]
###
#Making legends
custom_lines_color = []
for i in range(len(source_names)):
        custom_lines_color.append(Line2D([0], [0], color=colors[i], lw=4))
custom_lines_color.append(Line2D([0], [0], color='k', ls ='-'))
custom_lines_color.append(Line2D([0], [0], color='k', ls ='--'))

legend_names = list(data_dict.keys())
legend_names.append('Direct')
legend_names.append('Refracted')


custom_legend = []
for i in range(len(source_names)):
        custom_legend.append(Line2D([0], [0], color=colors[i], lw=4))

#colors[i], lw=4))
#new_custom = custom_lines_color.append(custom_lines_style)
#new_custom = custom_lines_color.append(custom_lines_style)
#print(new_custom)

legend_names = list(data_dict.keys())
#legend_names = []
#legend_names.append(source_names)
legend_names.append('Direct')
legend_names.append('Refracted')
#print(legend_names)
#print(source_names)
#print(custom_lines_color)
#exit()

#Variable arrays for plotting
hist_vars = ['rec_ang','theta_rec','view_ang','launch_ang','reflect_ang',
             'nnu_theta', 'nnu_phi',
             'dist', 'ShowerEnergy', 'depth', 'distance', 'flavor', 'elast', 'weight']
bins = [bin_cos, bin_cos, bin_cos, bin_cos, bin_cos, bindistance]
ang_strings = ['ang', 'theta', 'phi']
# var = ['trigg', 'weight', 'posnu_x', 'posnu_y', 'posnu_z',
#                'rec_ang_0', 'theta_rec_0', 'reflect_ang_0',
#                'dist_0', 'arrival_time_0', 'reflection_0', 
#                'l_att_0', 'view_ang_0', 'launch_ang_0',
#                'rec_ang_1', 'theta_rec_1', 'reflect_ang_1',
#                'dist_1', 'arrival_time_1', 'reflection_1', 
#                'l_att_1', 'view_ang_1', 'launch_ang_1',
#                'current', 'flavor', 'elast',
#                'nnu_theta', 'nnu_phi', 'ShowerEnergy',
#                'depth', 'distance']
        
'''
Below the first for statement below is my plotting function for histograms


'''

def hist_maker(hist_var):
        #print(hist_var)
        for i in range(len(source_names)):
                #print("Plotting...")
                #print(source_names[i])
                #print("...")
                try:    
                        if 'ang' in hist_var or 'theta' in hist_var or 'phi' in hist_var:
                                plt.hist(np.cos(data_dict[source_names[i]]['{0}_0'.format(hist_var)]), 
                                         weights=data_dict[source_names[i]]['weight'],bins=bin_cos, density=False, 
                                         histtype='step', color=colors[i], ls='-', label=str(source_names[i])+' direct')
                                plt.hist(np.cos(data_dict[source_names[i]]['{0}_1'.format(hist_var)]), 
                                         weights=data_dict[source_names[i]]['weight'], bins=bin_cos, density=False, 
                                         histtype='step', color=colors[i], ls='--', label=str(source_names[i])+' refracted')
                                plt.xlabel("Cos({0})".format(hist_var), fontsize=12)
                                
                        else:
                                plt.hist(data_dict[source_names[i]]['{0}_0'.format(hist_var)], 
                                         weights=data_dict[source_names[i]]['weight'],bins=bindistance, density=False, 
                                         histtype='step', color=colors[i], ls='-', label=str(source_names[i])+' direct')
                                plt.hist(data_dict[source_names[i]]['{0}_1'.format(hist_var)], 
                                         weights=data_dict[source_names[i]]['weight'], bins=bindistance, density=False, 
                                         histtype='step', color=colors[i], ls='--', label=str(source_names[i])+' refracted')
                                plt.xlabel("{0}".format(hist_var), fontsize=12)
                        
                                
                        legend_1 = plt.legend(custom_lines_color, legend_names, loc='best')
                        
                        plt.ylabel("Events", fontsize=12)
                        plt.grid(linestyle='--')
                        plt.gca().add_artist(legend_1)
                        plt.tight_layout()
                    
                        plt.savefig('test_plots/Hist_{0}_0_{0}_1_.png'.format(hist_var),dpi=300)
  
                except:
                        
                        if 'ang' in hist_var or 'theta' in hist_var or 'phi' in hist_var:
                                plt.hist(np.cos(data_dict[source_names[i]]['{0}'.format(hist_var)]), 
                                         weights=data_dict[source_names[i]]['weight'], bins=bin_cos, density=False, 
                                         histtype='step', color=colors[i], ls='-', label=str(source_names[i]))
                                plt.xlabel("Cos({0})".format(hist_var), fontsize=12)
                                legend_2 = plt.legend(custom_legend, source_names, loc='upper left')
                                
                        elif 'weight' in hist_var:
                                plt.hist(data_dict[source_names[i]]['{0}'.format(hist_var)], 
                                         log=True, density=False, 
                                         histtype='step', color=colors[i], ls='-', label=str(source_names[i]))#, bins =40)
                                plt.xlabel("{0}".format(hist_var), fontsize=12)
                                legend_2 = plt.legend(custom_legend, source_names, loc='upper center')

                        elif 'ShowerEnergy' in hist_var:
                                plt.hist(data_dict[source_names[i]]['{0}'.format(hist_var)],
                                         density=False, weights=data_dict[source_names[i]]['weight'],
                                         histtype='step', log=True, 
                                         color=colors[i], ls='-', label=str(source_names[i]))#, bins= )
                                plt.xlabel("{0}".format(hist_var), fontsize=12)
                                legend_2 = plt.legend(custom_legend, source_names, loc='upper left')

                        elif 'depth' in hist_var or 'distance' in hist_var:
                                plt.hist(data_dict[source_names[i]]['{0}'.format(hist_var)],
                                         density=False, weights=data_dict[source_names[i]]['weight'],
                                         histtype='step', 
                                         color=colors[i], ls='-', label=str(source_names[i]), bins= 40)
                                plt.xlabel("{0}".format(hist_var), fontsize=12)
                                legend_2 = plt.legend(custom_legend, source_names, loc='upper left')
                        
                        else:
                                plt.hist(data_dict[source_names[i]]['{0}'.format(hist_var)], 
                                         weights=data_dict[source_names[i]]['weight'],density=False, 
                                         histtype='step', color=colors[i], ls='-', label=str(source_names[i]))#, bins= )
                                plt.xlabel("{0}".format(hist_var), fontsize=12)
                                legend_2 = plt.legend(custom_legend, source_names, loc='best')
        
                        plt.ylabel("Events", fontsize=12)
                        plt.grid(linestyle='--')
                        plt.gca().add_artist(legend_2)
                        plt.tight_layout() 

                        plt.savefig('test_plots/Hist_{0}.png'.format(hist_var),dpi=300)
                                        
  
scatter_vars = ['distance', 'depth', 'rec_ang_0']
# print(scatter_vars[0])
# print(data_dict[source_names[0]][scatter_vars[0]])
# print(source_names)
# exit()
###
def scatter_maker(var1, var2):
        #print(var1+' '+var2)
        for i in range(len(source_names)):
                #print("Plotting...")
                #print(source_names[i])
                #print("...")
                # if 'ang' in scatter_vars or 'theta' in hist_vars or 'phi' in hist_vars:
                plt.figure(100+i, figsize=(8,6))
                plt.scatter(data_dict[source_names[i]]['{0}'.format(var1)], 
                            data_dict[source_names[i]]['{0}'.format(var2)], 
                            s=1.0, alpha=0.25, color=colors[i], label=str(source_names[i]))
                                
                plt.xlabel("{0}".format(var1), fontsize=12)
                plt.ylabel("{0}".format(var2), fontsize=12)
                
                # legend_3 = plt.legend([Line2D([0], [0], color=colors[i], lw=4)], source_names[i])#, loc='best')
                # plt.gca().add_artist(legend_3)
                
                plt.legend()
                plt.grid(linestyle='--')
                plt.tight_layout()
                plt.savefig('test_plots/Scatter_{2}_{0}_{1}_.png'.format(var1, var2, source_names[i]), dpi=300)
        
        

# plt.scatter(Adist_0, -Adepth, s=1.0, alpha=0.25, color='r', label='BICONE')
'''
can save same plot multiple times
'''
print("Histograms!")
for i in range(len(hist_vars)):
        plt.figure(i, figsize=(8,6))
        print("Plotting...")
        hist_maker(hist_vars[i])
        plt.clf()

print("Done!")
 

# print(custom_legend)
# print(custom_legend[0])
print("Scatter Plots!")
scatter_maker(scatter_vars[0], scatter_vars[1])
print("Done!")

'''
for i in range(len(source_names)):
scatter_maker(scatter_vars[0], scatter_vars[1]

''



# #############3
# ### 2D histograms
# #######

# bin_cos = np.linspace(-1.0, 1.0, 81)
# bin_dist = np.linspace(0,4000, 81)
# binsize = np.linspace(-1.0, 1.0, 81)
# bindepth = 20
# bindistance = np.linspace(0,4000, 81)

# print("2D histograms!!")

# ###Weighted
# #BICONE
# plt.figure(3001, figsize=(8,6))
# weight_AH_0 = plt.hist2d(Adist_0, np.cos(Arec_ang_0), bins=(bin_dist,bin_cos), weights=Aweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('BICONE: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_BICONE_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

# #GENETIS
# plt.figure(3003, figsize=(8,6))
# weight_BH_0 = plt.hist2d(Bdist_0, np.cos(Brec_ang_0), bins=(bin_dist,bin_cos), weights=Bweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('GENETIS: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_GENETIS_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

# weight_diff_0= weight_BH_0[0]- weight_AH_0[0]
# plt.figure(3005, figsize=(8,6))
# plt.pcolormesh(bin_dist, bin_cos, weight_diff_0.T, cmap='bwr')
# #plt.hist2d(diff, bins=(bin_dist,bin_cos))#, weights=Bweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('DIFF: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_DIFF_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

#D_vars =['distance', 'depth']


'''
function for scatter plot
'''

# for j in range(len(hist_vars)):
#         print("{0}".format(hist_vars[j]))
#         print("Plotting...")
#         plt.figure(j, figsize=(8,6))
#         #plt.title('', fontsize=16, fontweight="bold")
#         for i in range(len(source_names)):
#                 #print("Plotting...")
#                 #print(source_names[i])
#                 #print("...")
#                 try:
#                         if 'ang' in hist_vars[j] or 'theta' in hist_vars[j] or 'phi' in hist_vars[j]:
#                                 plt.hist(np.cos(data_dict[source_names[i]]['{0}_0'.format(hist_vars[j])]), 
#                                          weights=data_dict[source_names[i]]['weight'],bins=bin_cos, density=False, 
#                                          histtype='step', color=colors[i], ls='-', label=str(source_names[i])+' direct')
#                                 plt.hist(np.cos(data_dict[source_names[i]]['{0}_1'.format(hist_vars[j])]), 
#                                          weights=data_dict[source_names[i]]['weight'], bins=bin_cos, density=False, 
#                                          histtype='step', color=colors[i], ls='--', label=str(source_names[i])+' refracted')
#                                 plt.xlabel("Cos({0})".format(hist_vars[j]), fontsize=12)
                        
#                         else:
#                                 plt.hist(data_dict[source_names[i]]['{0}_0'.format(hist_vars[j])], 
#                                          weights=data_dict[source_names[i]]['weight'],bins=bindistance, density=False, 
#                                          histtype='step', color=colors[i], ls='-', label=str(source_names[i])+' direct')
#                                 plt.hist(data_dict[source_names[i]]['{0}_1'.format(hist_vars[j])], 
#                                          weights=data_dict[source_names[i]]['weight'], bins=bindistance, density=False, 
#                                          histtype='step', color=colors[i], ls='--', label=str(source_names[i])+' refracted')
#                                 plt.xlabel("{0}".format(hist_vars[j]), fontsize=12)
                                
#                         legend_1 = plt.legend(custom_lines_color, legend_names, loc='best')
                                
#                         plt.ylabel("Events", fontsize=12)
#                         plt.grid(linestyle='--')
#                         plt.gca().add_artist(legend_1)
#                         plt.tight_layout()
#                         plt.savefig('test_plots/Hist_{0}_0_{0}_1_.png'.format(hist_vars[j]),dpi=300)

#                 except:
#                         if 'ang' in hist_vars[j] or 'theta' in hist_vars[j] or 'phi' in hist_vars[j]:
#                                 plt.hist(np.cos(data_dict[source_names[i]]['{0}'.format(hist_vars[j])]), 
#                                          weights=data_dict[source_names[i]]['weight'],bins=bin_cos, density=False, 
#                                          histtype='step', color=colors[i], ls='-', label=str(source_names[i]))
#                                 plt.xlabel("Cos({0})".format(hist_vars[j]), fontsize=12)
                        
#                         elif 'weight' in hist_vars[j]:
#                                 plt.hist(data_dict[source_names[i]]['{0}'.format(hist_vars[j])], 
#                                          bins=bindistance, density=False, 
#                                          histtype='step', color=colors[i], ls='-', label=str(source_names[i]))
#                                 plt.xlabel("{0}".format(hist_vars[j]), fontsize=12)
                                 
                                
#                         else:
#                                 plt.hist(data_dict[source_names[i]]['{0}'.format(hist_vars[j])], 
#                                          weights=data_dict[source_names[i]]['weight'],bins=bindistance, density=False, 
#                                          histtype='step', color=colors[i], ls='-', label=str(source_names[i]))
#                                 plt.xlabel("{0}".format(hist_vars[j]), fontsize=12)
                                
#                         legend_2 = plt.legend(custom_legend, source_names, loc='best')
                                
#                         plt.ylabel("Events", fontsize=12)
#                         plt.grid(linestyle='--')
#                         plt.gca().add_artist(legend_2)
#                         plt.tight_layout()
#                         plt.savefig('test_plots/Hist_{0}.png'.format(hist_vars[j]),dpi=300)


# print("Done!")

'''
scatter variables
distance
depth
thats all i need really
i guess i could put in the other ones if I wanted to 
woud use the same if statement and try and except block as above for that

'''



'''
Now need to make a scatter plot function
then make a 2d histogram function
maybe use np.histgoram
plot it individually 
but if more than two data sets, need to compare 
then use the colormesh

'''
#legend_1 = plt.legend(custom_lines_color, source_names, loc='best')
#legend_2 = plt.legend(custom_lines_style, ['Direct', 'Refracted'], loc='upper center')
'''
how to add other variablses to plot into the histogram loop? 
could add try and except block :
is this  bad knowingly some will not work? 
idea is that for nnu and shower energy varibales, then it wfll fail the abov
but go straight to the except block
then I can have one loop for histogram variables
above is pretty much a function

'''
# from matplotlib.lines import Line2D
# cmap = plt.cm.coolwarm
# custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
#                 Line2D([0], [0], color=cmap(.5), lw=4),
#                 Line2D([0], [0], color=cmap(1.), lw=4)]


# cmap = plt.cm.coolwarm
# custom_lines = [Line2D([0], [0], color='r', lw=4),
#                 Line2D([0], [0], color='b', lw=4),
#                 Line2D([0], [0], color='g', lw=4),
#                 Line2D([0], [0], color='c', lw=4)]

'''
make empty array for custom_lines
use for loop to append to empty array 
then use custom_lines for my legend
'''


'''
need to figure out how to add other legend
only have for direct (solid) & refracted(dashed)
'''
# # #fig, ax = plt.subplots()
# # #lines = ax.plot(data)

#plt.legend(loc='upper right')
# #best')


# # lines = plt.get_lines()
# # legend1 = plt.legend([lines[i] for i in [0,1,2]], ["algo1", "algo2", "algo3"], loc=1)
# # legend2 = plt.legend([lines[i] for i in [0,3,6]], parameters, loc=4)
# plt.add_artist(legend1)
# plt.add_artist(legend2)


'''
maybe make above a function
call '{0}_0'format(variable_array[i])
where variable_array has common var names
like rec_ang, dist, arrival_time, etc
then I can use that common name in the name of the plot filename
also use it for the tite possibly
also use it for all other labels
then if I make more of these types for other plots,
then I can just add to the loop
for bin size, 
i can change variale name before loop
'''        






#######################################
###General Info
#######################################
for i in range(len(source_names)):
        print('#'*28)
        print('\033[1;37m{0}\033[0;0m'.format(source_names[i]))
        print('#'*28)
        print('\033[4;37mEvents\033[0;0m')
        print('Triggered: \033[1;31m{0}\033[0;0m'.format(len(data_dict[source_names[i]]['trigg'])))
        print('Usable: \033[1;31m{0}\033[0;0m'.format(len(data_dict[source_names[i]]['weight'])))
        print('Weighted: \033[1;31m{0}\033[0;0m'.format(np.sum(data_dict[source_names[i]]['weight'])))
        print('#'*50)
        print('\n')
stop = timeit.default_timer()
print('Time: \033[1;31m{0}\033[0;0m'.format(stop - start))
exit()


 
#                  plt.hist(data_dict[source_names[i][k], weights = data_dict[source_names[i]['weight'], bins

# ## can give name as arguemnt for 
# '''
# give different arguements for plotting function
# maybe only od the plt.hist
# use larger function or loop
# '''
# def simple_plotmaker(var1):
#         plt.hist(var1)
#         name = str(var1)
#         print(name)
#         plt.savefig('test_plots/histogram.png')
#         return

# hist_vars = ['rec_ang_0', 'rec_ang_1']
# for k in hist_vars:
#         for i in range(len(data_dict.keys())):
#                 plt.hist(data_dict[source_names[i][k], weights = data_dict[source_names[i]['weight'], bins

# plt.legend()#loc='upper left')
# plt.grid(linestyle='--')
# plt.xlabel("Depth (m)", fontsize=12)
# plt.ylabel("Events", fontsize=12)
# plt.title('BICONE vs GENETIS', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/BICONE_GENETIS_Depth.png',dpi=300)




# plt.hist(Adepth, bins=bindepth, weights=Aweight, density=False, histtype='step', color='r', ls='-', label='BICONE')
#simple_plotmaker(data_dict['source']['rec_ang_0'])


# print('weight')
# print(len(data_dict['source']['weight']))
# print(len(data_dict['source']['trigg']))
# print(len(data_dict['source_2']['trigg']))


# stop = timeit.default_timer()
# print('Time: \033[1;31m{0}\033[0;0m'.format(stop - start))

# #extra
# from astropy.table import QTable, Table, Column
# from astropy import units as u

# table = Table(data_dict['source'][1:])
# print(table)



# print(locals().keys())
# print(globals().keys())

#######################################
###Plots
#######################################

###
w = 2.0
binsize = np.linspace(-1.0, 1.0, 41)
bindepth = 20
bindistance = np.linspace(0,4000, 21)

bin_cos = np.linspace(-1.0, 1.0, 41)
bin_dist = np.linspace(0,4000, 41)
binsize = np.linspace(-1.0, 1.0, 41)
bindepth = 20
bindistance = np.linspace(0,4000, 41)

print("Help is on on the way, dear!!")
###
# print("Histogram time...")
# #Depth
# plt.figure(101, figsize=(8,6))
# # (n_Areflect_ang_1, bins_Areflect_ang_0, patches_Areflect_ang_1) = plt.hist(np.cos(Areflect_ang_1), weights=Aweight, density=False, histtype='step', color='r', ls='--', label='BICONE Refracted')
# # (n_Breflect_ang_1, bins_Breflect_ang_1, patches_Breflect_ang_1) = plt.hist(np.cos(Breflect_ang_1), weights=Bweight, density=False, histtype='step', color='b', ls='--', label='GENETIS Refracted')


# plt.hist(Bdepth, bins=bindepth, weights=Bweight, density=False, histtype='step', color='b', ls='-', label='GENETIS')

# plt.legend()#loc='upper left')
# plt.grid(linestyle='--')
# plt.xlabel("Depth (m)", fontsize=12)
# plt.ylabel("Events", fontsize=12)
# plt.title('BICONE vs GENETIS', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/BICONE_GENETIS_Depth.png',dpi=300)

# print("plotting...")

# #Both
# fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True, sharex=True)
# axs[0].scatter(Adisplace, -Adepth, s=1.0, alpha=0.25, color='r', label='BICONE')

# # axs[0].legend()#loc='upper left')
# axs[0].grid(linestyle='--')
# axs[0].set_xlabel("XY Distance (m)", fontsize=8)
# axs[0].set_ylabel("Depth (m)", fontsize=8)

# axs[0].set_title('BICONE: XY Distance vs. Depth (m)', fontsize=10, fontweight="bold")

# axs[1].scatter(Bdisplace, -Bdepth, s=1.0, alpha=0.25, color='b', label='GENETIS')

# axs[1].grid(linestyle='--')
# axs[1].set_title('GENETIS: XY Distance vs. Depth (m)', fontsize=10, fontweight="bold")

# axs[1].set_xlabel("XY Distance (m)", fontsize=8)
# #axs[1].set_ylabel("Depth (m)", fontsize=8)
# plt.tight_layout()
# plt.savefig("test_plots/TEST_BOTH_SUBPLOTS_XY_DISTANCE_DEPTH.png", dpi=300)

# print("plotting...")

# #################################
# ##scatter 


# ####Dist_0
# ##Depth vs XY Distance
# plt.figure(401, figsize=(8,6))
# plt.scatter(Adist_0, -Adepth, s=1.0, alpha=0.25, color='r', label='BICONE')
# plt.legend()#loc='upper left')
# plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Depth (m)", fontsize=12)
# plt.title('BICONE: Depth vs. Distance', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/BICONE_depth_dist_0.png', dpi=300)

# print("plotting...")

# ###
# #Color Plots

# ## Dist vs theta_rec


# ##BICONE
# #theta_rec_0
# print("Praise the color plot!")

# plt.figure(901, figsize=(8,6))

# plt.scatter(np.array(Adist_0)[np.where(Adepth <= 200)], np.cos(np.array(Atheta_rec_0)[np.where(Adepth <= 200)]), s=1.0, alpha=0.25)#, color='r')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 200) & (Adepth <= 400) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 200) & (Adepth <= 400) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 400) & (Adepth <= 600) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 400) & (Adepth <= 600) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 600) & (Adepth <= 800) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 600) & (Adepth <= 800) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 800) & (Adepth <= 1000) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 800) & (Adepth <= 1000) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 1000) & (Adepth <= 1200) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 1000) & (Adepth <= 1200) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 1200) & (Adepth <= 1400) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 1200) & (Adepth <= 1400) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 1400) & (Adepth <= 1600) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 1400) & (Adepth <= 1600) )]), s=1.0, alpha=0.25)#, color='g')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 1600) & (Adepth <= 1800) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 1600) & (Adepth <= 1800) )]), s=1.0, alpha=0.25)#, color='b')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 1800) & (Adepth <= 2000) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 1800) & (Adepth <= 2000) )]), s=1.0, alpha=0.25)#, color='b')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 2000) & (Adepth <= 2200) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 2000) & (Adepth <= 2200) )]), s=1.0, alpha=0.25)#, color='c')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 2200) & (Adepth <= 2400) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 2200) & (Adepth <= 2400) )]), s=1.0, alpha=0.25)#, color='c')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 2400) & (Adepth <= 2600) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 2400) & (Adepth <= 2600) )]), s=1.0, alpha=0.25)#, color='m')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 2600) & (Adepth <= 2800) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 2600) & (Adepth <= 2800) )]), s=1.0, alpha=0.25)#, color='k')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 2800) & (Adepth <= 3000) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 2800) & (Adepth <= 3000) )]), s=1.0, alpha=0.25)#, color='y')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 3000) & (Adepth <= 3200) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 3000) & (Adepth <= 3200) )]), s=1.0, alpha=0.25)#, color='y')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 3200) & (Adepth <= 3400) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 3200) & (Adepth <= 3400) )]), s=1.0, alpha=0.25)#, color='orange')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 3400) & (Adepth <= 3600) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 3400) & (Adepth <= 3600) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 3600) & (Adepth <= 3800) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 3600) & (Adepth <= 3800) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 3800) & (Adepth <= 4000) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 3800) & (Adepth <= 4000) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 4000) & (Adepth <= 4200) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 4000) & (Adepth <= 4200) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 4200) & (Adepth <= 4400) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 4200) & (Adepth <= 4400) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 4400) & (Adepth <= 4600) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 4400) & (Adepth <= 4600) )]), s=1.0, alpha=0.25)#, color='indigo')
# plt.scatter(np.array(Adist_0)[np.where( (Adepth > 4600) & (Adepth <= 4800) )], np.cos(np.array(Atheta_rec_0)[np.where( (Adepth > 4600) & (Adepth <= 4800) )]), s=1.0, alpha=0.25)#, color='indigo')

# #plt.legend()#loc='upper left')
# plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(theta_rec_0)", fontsize=12)
# plt.title('BICONE: Distance vs. Cos(theta_rec_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/TEST_BICONE_Dist_theta_rec_0.png', dpi=300)

# print("plotting...")

# #############3
# ### 2D histograms
# #######

# bin_cos = np.linspace(-1.0, 1.0, 81)
# bin_dist = np.linspace(0,4000, 81)
# binsize = np.linspace(-1.0, 1.0, 81)
# bindepth = 20
# bindistance = np.linspace(0,4000, 81)

# print("2D histograms!!")

# ###Weighted
# #BICONE
# plt.figure(3001, figsize=(8,6))
# weight_AH_0 = plt.hist2d(Adist_0, np.cos(Arec_ang_0), bins=(bin_dist,bin_cos), weights=Aweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('BICONE: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_BICONE_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

# #GENETIS
# plt.figure(3003, figsize=(8,6))
# weight_BH_0 = plt.hist2d(Bdist_0, np.cos(Brec_ang_0), bins=(bin_dist,bin_cos), weights=Bweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('GENETIS: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_GENETIS_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

# weight_diff_0= weight_BH_0[0]- weight_AH_0[0]
# plt.figure(3005, figsize=(8,6))
# plt.pcolormesh(bin_dist, bin_cos, weight_diff_0.T, cmap='bwr')
# #plt.hist2d(diff, bins=(bin_dist,bin_cos))#, weights=Bweight)
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("Distance (m)", fontsize=12)
# plt.ylabel("Cos(rec_ang_0)", fontsize=12)
# plt.title('DIFF: Distance vs. Cos(rec_ang_0)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_DIFF_2DHist_dist_0_rec_ang_0.png', dpi=300)

# print("plotting...")

# ##Depth Distancc
# A_distance = np.sqrt(np.array(Adist_0)**2 - Adepth**2)
# B_distance = np.sqrt(np.array(Bdist_0)**2 - Bdepth**2)
# binsize = 80
# ###Scatter 

# plt.figure(6001, figsize=(8,6))
# plt.scatter(A_distance, -Adepth, s=1.0, alpha=0.25, color='r')
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("XY Distance (m)", fontsize=12)
# plt.ylabel("Depth (m)", fontsize=12)
# plt.title('BICONE: Distance vs. Depth (m)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_test_BICONE_scatter_distance_depth.png', dpi=300)


# plt.figure(6002, figsize=(8,6))
# plt.scatter(B_distance, -Bdepth, s=1.0, alpha=0.25, color='b')
# #plt.legend()#loc='upper left')
# #plt.grid(linestyle='--')
# plt.xlabel("XY Distance (m)", fontsize=12)
# plt.ylabel("Depth (m)", fontsize=12)
# plt.title('GENETIS: Distance vs. Depth (m)', fontsize=16, fontweight="bold")
# plt.tight_layout()
# plt.savefig('test_plots/weighted_test_GENETIS_scatter_distance_depth.png', dpi=300)

# from astropy.table import Table, Column
# from astropy.io import ascii
# A_table = Table()
# A_table['x'] = Aposnu_x
# A_table['y'] = Aposnu_y
# A_table['z'] = Aposnu_z

# A_table.write('table.csv', format='csv', overwrite=True)  

