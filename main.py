#/usr/bin/env python3
# -*- coding: utf-8 -*-

#from os import times
from datetime import date
from operator import truediv
from os import write
from pickle import FALSE, TRUE
import re
import pandas as pd
import numpy as np
from scipy import integrate
import csv

class reader:

    def gain_condition(file_out_log):
        f = open(file_out_log,mode = "r")
        lines = f.readlines()

        for i,line in enumerate(lines):

            if "dynamic parameters and flags:" in line:
                fragments = (lines[i + 3].split())[-1]
                ene_kin_ini = float((lines[i + 5].split())[-1])

                #print(fragments)

                nucleus = {}
                for fragment in range(1,int(fragments)+1):
                    key_n = "N" + str(fragment)
                    key_z = "Z" + str(fragment)
                    nucleus[key_n] = int((lines[i + 12 + fragment].split())[-2])
                    nucleus[key_z] = int((lines[i + 12 + fragment].split())[-1])

        f.close

        reducedmass = 938.0 * (float((nucleus["N1"] + nucleus["Z1"]) * (nucleus["N1"] + nucleus["Z1"])))\
            /(float((nucleus["N1"] + nucleus["Z1"]) + (nucleus["N1"] + nucleus["Z1"])))

        v_ini = np.sqrt(2*ene_kin_ini/reducedmass)
        ene_kin_ini = ene_kin_ini/sum(nucleus.values())

        conditions = list(nucleus.values())
        conditions.append(ene_kin_ini)
        conditions.append(v_ini)

        return conditions,reducedmass,nucleus

    def gain_plotrange(file_out_log):

        f = open(file_out_log,mode = "r")
        lines = f.readlines()

        #investigate plot style
        for i,line in enumerate(lines):

            if "Contour 5=      " in line:
                top_plot = i + 2

                range_plot_width = [lines[i + 3].find("i"),lines[i + 3].rfind("i")]
                #print(type(range_plot_width[0]))

        
            if "### Writeden:  F" in line:
                bottom_plot = i - 2

            if "top_plot" in locals() and "bottom_plot" in locals():
                range_plot_vertical = [0,bottom_plot - top_plot]
                #print(range_plot)
                break

        f.close

        #print()

        range_plot = [range_plot_vertical,range_plot_width]

        return range_plot

    def gain_plot(plot_range_vertical,file_out_log):

        #data plot has the plot graph at each time
        f = open(file_out_log,mode = "r")
        lines = f.readlines()
        plots = []

        for i,line in enumerate(lines):

            if "cpu time(sec)  =       0.0000" in line:
                plots.append(lines[i+8:i+8+plot_range_vertical[1]+1])

        plots = dict(zip(range(len(plots)),plots))

        return plots

    def gain_blank(plots,plot_range):

        blank_alltime = []

        for plot in plots.values():
            plot_listed = []

            for plot_line in plot:
                plot_listed.append(list(plot_line.rstrip('\n')))

            plot_pd = pd.DataFrame(plot_listed)
            plot_pd = plot_pd.drop(plot_range[0],axis=0)
            plot_pd = plot_pd.drop(list(range(0,plot_range[1][0] + 1)),axis=1)
            plot_pd = plot_pd.drop(plot_range[1][1],axis=1)

            range_nuclei_eachtime = []
            for i,line_vertical in plot_pd.iteritems():

                if "1" in line_vertical.tolist():
                    range_nuclei_eachtime.append(i)

            blank_eachtime = []
            for i in range(min(range_nuclei_eachtime),max(range_nuclei_eachtime)+1):
                if all([string == " " for string in plot_pd[i].tolist()]):
                    blank_eachtime.append(i)
            
            blank_alltime.append(blank_eachtime)

        
        return blank_alltime

class writer:
    def clear_summary():
        f = open("summary.csv",mode="w")
        writer_summary = csv.writer(f)

        header = ["N1","Z1","N2","Z2"]
        header.append("E_ini(MeV)")
        header.append("v_ini(c)")
        header.append("status")
        header.append("num")
        header.append("E_diss(MeV)")
        header.append("integral_Rdot^2(fm*c)")
        header.append("fric_ave(MeV/(fm*c))")
        #[header.append(i) for i in ["E_ini(MeV)","v_ini(c)","status","num","E_diss(MeV)","integral_Rdot^2(fm*c)","fric_av(MeV/(fm*c))"]]
        writer_summary.writerow(header)
        f.close

    def summary_write(reaction_status,fric_ave_result):
        f = open("summary.csv",mode="a")
        writer_summary = csv.writer(f)

        result = []

        if reaction_status["scattering"]:
            status = "scattering"
            result.extend(conditions)
            result.append(status)
            result.append(2)
            result.append("not calculated")
            writer_summary.writerow(result)
        elif reaction_status["fusion"]:
            status = "fusion"
            result.extend(conditions)
            result.append(status)
            result.append(1)
            result.append("not calculated")
            writer_summary.writerow(result)
        elif reaction_status["fission"]:
            status = "fission"
            result.extend(conditions)
            result.append(status)
            result.append(num_fragment_alltime[-1])
            [result.append(i) for i in fric_ave_result]
            writer_summary.writerow(result)
        
        f.close

    def clear_detail():
        f = open("summary_detail.txt",mode="w")

        str_start = '==============calculation start=============\n'
        f.write(str_start)
        f.close

    def detail_write(reaction_status,file_out_log,file_out_cm,conditions,reducedmass,nucleus,plots,fric_ave_detail,ene):

        f = open("summary_detail.txt",mode="a")
                
        #print(f_name_ene)
        #f_ene = open("")

        f.write("file_out_log = " + file_out_log + '\n')
        f.write("file_cm = " + file_out_cm + '\n')

        str_nucleus = ""

        for num_nucleus in range(len(list(nucleus.keys()))):
            str_nucleus = str_nucleus + list(nucleus.keys())[num_nucleus] + ":" + str(list(nucleus.values())[num_nucleus]) + ","

        str_nucleus = str_nucleus[:-1]

        f.write("nucleus = " + str_nucleus + '\n')

        f.write("reduced mass = " + str(reducedmass) + '\n')
        f.write("E_ini,v_ini = " + str(conditions[4]) + "," + str(conditions[5]) + '\n')
        f.write('\n')

        if reaction_status["scattering"]:
            f.write("status = scattering")
            f.write('\n')
            [f.write(plot_line) for plot_line in plots[0]]
            f.write('\n')
            [f.write(plot_line) for plot_line in plots[len(plots) - 1]] 
            f.write('\n')

        elif reaction_status["fusion"]:
            f.write("status = fusion")
            f.write('\n')
            [f.write(plot_line) for plot_line in plots[0]]
            f.write('\n')
            [f.write(plot_line) for plot_line in plots[len(plots) - 1]]
            f.write('\n')

        elif reaction_status["fission"]:
            f.write("status = fission")
            f.write('\n')
            f.write("fusion point : index , time, R = " + str(fric_ave_detail[0]) + ", " + str(fric_ave_detail[1]) + ", " + str(fric_ave_detail[2]) + '\n')
            [f.write(plot_line) for plot_line in plots[fric_ave_detail[0]]]
            f.write('\n')
            f.write("kinetic energy,potential energy = " + str(fric_ave_detail[3]) + ", " + str(fric_ave_detail[4]) + '\n')

            f.write('\n')

            f.write("fission point : index , time, R = " + str(fric_ave_detail[5]) + ", " + str(fric_ave_detail[6])  + ", " + str(fric_ave_detail[7]) + '\n')
            [f.write(plot_line) for plot_line in plots[fric_ave_detail[5]]]
            f.write('\n')
            f.write("kinetic energy,potential energy = " + str(fric_ave_detail[8]) + ", " + str(fric_ave_detail[9]) + '\n')
            f.write('\n')

            #make energyfile
            f_name_ene = "pot" + "_" + str(list(nucleus.values())[0]) + "-" + str(list(nucleus.values())[1]) + "to" + str(list(nucleus.values())[2]) + "-" + str(list(nucleus.values())[3]) + "_" + str(conditions[4]) + ".csv" 
            f_ene = open(f_name_ene, "w",newline='')
            header_ene = ["time(fm/c)","R(fm)","nuclear potential(MeV)"]
            writer = csv.writer(f_ene)
            writer.writerow(header_ene)
            writer.writerows(ene)
            f_ene.close()



        f.close

class reaction:
    def count_fragment(blank_alltime):

        num_fragment_alltime = []
        for blank_eachtime in blank_alltime:

            if len(blank_eachtime) == 0:
                num_fragment_eachtime = 1
                num_fragment_alltime.append(num_fragment_eachtime)

            else:
                num_fragment_eachtime = 2
                for i in range(len(blank_eachtime)-1):
                    if blank_eachtime[i+1] - blank_eachtime[i] != 1:
                        num_fragment_eachtime = num_fragment_eachtime + 1

                num_fragment_alltime.append(num_fragment_eachtime)
                
        return num_fragment_alltime
    
    def calc_fric_pot(file_out_cm,reducedmass,nucleus,num_fragment_alltime):

        datas = np.loadtxt(file_out_cm,dtype="float",usecols=[0,1],unpack=True)
        datas[1] = datas[1]*2

        #serch before fusion time and after fission time
        is_fragment_alltime = [num_fragment_eachtime == 1 for num_fragment_eachtime in num_fragment_alltime]
        index_bfu = is_fragment_alltime.index(True) - 1
        index_afi = is_fragment_alltime.index(True) + is_fragment_alltime[is_fragment_alltime.index(True):].index(False)

        #calculate dissipation energy
        ##calculate dR/dt
        cm_differential_first = (np.diff(datas[1],n = 1))/(np.diff(datas[0],n = 1))

        K_bfu = 0.5 * reducedmass * (cm_differential_first[index_bfu]**2)
        V_bfu = (float(nucleus["Z1"]) * float(nucleus["Z2"]) * (197.3/137.00) / datas[1][index_bfu])
        E_bfu = K_bfu + V_bfu

        K_afi = 0.5 * reducedmass * (cm_differential_first[index_afi]**2)
        V_afi = (float(nucleus["Z1"]) * float(nucleus["Z2"]) * (197.3/137.00) / datas[1][index_afi])
        E_afi = K_afi + V_afi

        E_diss = E_bfu - E_afi

        #calculate integral (dR/dt)^2
        #cm_diff_2 = cm_differential_first**2
        cm_differential_first_integrated = integrate.simps((cm_differential_first**2)[index_bfu+1:index_afi-1],datas[0][index_bfu+1:index_afi-1])

        #calculate average friction coefficient
        fric_ave = E_diss / cm_differential_first_integrated
        #print(fric_ave)

        fric_ave_result = [E_diss,cm_differential_first_integrated,fric_ave]
        fric_ave_detail = [index_bfu,datas[0][index_bfu],datas[1][index_bfu],K_bfu,V_bfu,index_afi,datas[0][index_afi],datas[1][index_afi],K_afi,V_afi]

        # calculated V_nucle in terms of the equation of motion
        cm_differential_second = (np.diff(cm_differential_first,n = 1))/(np.diff(datas[0][:-1],n = 1))
        V_differential_byt = (-fric_ave * cm_differential_first[:-1] - reducedmass * cm_differential_second) * cm_differential_first[:-1]

        V = []
        V.append(V_bfu)
        [V.append(V[i - 1] - (V_differential_byt[i +index_bfu - 1]) * (datas[0][i + index_bfu] - datas[0][i + index_bfu - 1])) for i in range(1,index_afi - index_bfu - 1)]

        V_ele = []
        [V_ele.append((float(nucleus["Z1"]) * float(nucleus["Z2"]) * (197.3/137.00) / datas[1][i])) for i in range(index_bfu,index_afi - 1)]

        #V_nucle = []
        V_nucle = [i - j for (i, j) in zip(V, V_ele)]
  
        ene = []
        ene = list(zip(datas[0][index_bfu+1:index_afi],datas[1][index_bfu+1:index_afi],V_nucle))

        #print(ene)
        #print(V_ele)
        #print(K)
        #print(V_nucle)
        #print(E_diss_all)
        #print((cm_differential_first**2)[index_bfu+1:index_bfu+3])
        #print(V)
        #print(V_ele)
        #print(len(V))
        #print(len(V_ele))

        return fric_ave_result,fric_ave_detail,ene

#=====main========

writer.clear_summary()
writer.clear_detail()


file_list = open("out_and_cm.txt",mode="r")

files = file_list.readlines()

for file in files:

    file_out_cm = file.split()[1]
    file_out_log = file.split()[0]
    print(file_out_cm)
    print(file_out_log)

    nucleus = {}
    conditions,reducedmass,nucleus = reader.gain_condition(file_out_log)

    plot_range = reader.gain_plotrange(file_out_log)
    plots = reader.gain_plot(plot_range[0],file_out_log)
    blank_alltime = reader.gain_blank(plots,plot_range)
    num_fragment_alltime = reaction.count_fragment(blank_alltime)
    fric_ave_result = []
    fric_ave_detail = []

    is_scattering = all([num_fragment_eachtime == 2 for num_fragment_eachtime in num_fragment_alltime])
    is_fusion = (num_fragment_alltime[-1] == 1)
    
    if not (is_scattering or is_fusion):
        is_fission = True
        fric_ave_result,fric_ave_detail,ene = reaction.calc_fric_pot(file_out_cm,reducedmass,nucleus,num_fragment_alltime)
    else:
        is_fission = False
        fric_ave_result = []
        fric_ave_detail = []
        ene = []

    reaction_status = {"scattering":is_scattering,"fusion":is_fusion,"fission":is_fission}

    writer.summary_write(reaction_status,fric_ave_result)
    writer.detail_write(reaction_status,file_out_log,file_out_cm,conditions,reducedmass,nucleus,plots,fric_ave_detail,ene)

    print(conditions)
    print(fric_ave_result)






