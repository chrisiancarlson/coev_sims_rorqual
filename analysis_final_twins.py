import numpy as np 
import matplotlib
import itertools
import sys
import json
import multiprocessing as mp
import os
#import scipy
current_directory = os.getcwd()
import re
import csv
import random
from matplotlib import pyplot as plt
import seaborn as sns


######## import full output of SLiM simulation: #################
def import_output_full(path):
   # input path to full output and a character identifier 
    
    # import full output and return two dictionaries. 
#First dictionary is called mut_dict.
    #for  each [mutation_key]: [effect, origin time, and number of genomes]
# Second dictionary is [individual_key]: [[x,y],[genome1ID,genome2ID],[[mutations on genome1],[mutations on genome2]] ]

    file=open(path,"r")
    space_dict={}
    mut_dict={}

    flagm = False
    flagi = False
    flagg=False
    for vals in file:
        new=vals.split()
        
        if new[0]=="Mutations:":
            flagm = True
            continue
        if flagm == True: 
            if new[0]=="Individuals:": 
                flagm=False
                flagi=True
                continue
            #output the mutation's effect, time it originated, number of genomes it's in, and its position
            mut_dict[new[0]]=[new[4],new[7],new[8],new[3]]
            #print(mut_dict)
            #sys.exit()
        if flagi == True: 
            if new[0]=="Genomes:": 
                flagi=False
                flagg=True
                continue
                #break
            ### get spatial positions of individuals
            x=float(new[4])
            y=float(new[5])
            id=new[0]
            # an individual's x y position, and the genomes that it has. 
            space_dict[id]=[[x,y],[new[2],new[3]],[]]
        
        if flagg==True: 
            if len(new)> 2:
                id=new[0]
                #print(id)
                # if an individual carries the right genomes, append the mutations that individual
                # carries to the space dictionary in the third position
                for ind in space_dict: 
                    if id==space_dict[ind][1][0] or id==space_dict[ind][1][1]:
                        space_dict[ind][2].append(new[2:])
    return space_dict,mut_dict


######## sample random individuals from each species: #################
# sample a radius around a focal individual 'num' times [ind1,ind2]
# sigma = deme radius, num = total deme number
def random_sampler(space_dict_m,space_dict_p,sigma,num):
    keys=[key for key in space_dict_m.keys()]
    
    sample_num=num//2 ###### center half the samples on moths
    choices=0
    samples_m=[]
    while choices < sample_num:
        choice=np.random.choice(keys,1,replace=False)
        
        x,y=space_dict_m[choice[0]][0]
        # make sure we're not too close to a spatial boundary b/c I don't want to deal with periodicity
        if x>0+sigma and x<10-sigma and y>0+sigma and y<10-sigma: 
            
            for plant in space_dict_p: 
            # make sure there will be a partner: 
                if (space_dict_p[plant][0][0]-x)**2+ (space_dict_p[plant][0][1]-y)**2<=sigma**2: 
                    samples_m.append(choice[0])
                    choices+=1
                    break

     ###### center the other half the samples on plants
    choices=0
    samples_p=[]
    keys2=[key for key in space_dict_p.keys()]
    while choices < sample_num:
        choice=np.random.choice(keys2,1,replace=False)
        x,y=space_dict_p[choice[0]][0]
        if x>0+sigma and x<10-sigma and y>0+sigma and y<10-sigma: 
            for moth in space_dict_m: 
            # make sure there will be a partner: 
                if (space_dict_m[moth][0][0]-x)**2+ (space_dict_m[moth][0][1]-y)**2<=sigma**2: 
                    samples_p.append(choice[0])
                    choices+=1
                    break

    positions=[]
    for sample in samples_m: 
        positions.append(space_dict_m[sample][0])
    
    for sample in samples_p: 
        positions.append(space_dict_p[sample][0]) 
    return positions # positions around which to sample


######## Find demes around focal individuals: #################

# return a list of the individuals around a focal individual [[inds in deme 1],[inds in deme 2]]
# positions = xy positions of focal individuals
def find_individuals(positions,space_dict,sigma): 
    deme_list=[[] for pos in positions]
    for i,pos in enumerate(positions): 
        # search for the individuals in that deme: 
        for ind in space_dict:
            indpos=space_dict[ind][0]
            # see if an individual is within radius sigma of the focal spatial position: 
            if (indpos[0]-pos[0])**2 + (indpos[1]-pos[1])**2<=sigma**2: 
                deme_list[i].append(ind)
    return deme_list

######## Find mean phenotype of deme: #################

# return a list of lists of the phenotypes in a given deme [ [zs in deme 1], [zs in deme 2]]
# deme_list = list of individuals that make up demes
def phenotypes(deme_list,space_dict,mut_dict):
    z_list=[[] for deme in deme_list]
    for i,deme in enumerate(deme_list): 
        # search for the individuals in that deme: 
        for ind in deme:     
            phenotype=0
            for mut_list in space_dict[ind][2]: 
                for mut in mut_list: 
                    phenotype+= float(mut_dict[mut][0])
            z_list[i].append(phenotype)
    return z_list

######## Calculate correlation between two vectors: #################

def calculate_correlation(x, y):
    """
    Calculates the Pearson correlation coefficient between two vectors.

    Args:
        x (list or tuple): The first vector.
        y (list or tuple): The second vector.

    Returns:
        float: The Pearson correlation coefficient, or None if an error occurs.
    """
    n = len(x)
    if n != len(y) or n <= 1:
        return None  # Vectors must have the same length and at least two elements
    
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x_squared = sum(xi ** 2 for xi in x)
    sum_y_squared = sum(yi ** 2 for yi in y)

    numerator = n * sum_xy - sum_x * sum_y
    denominator = ((n * sum_x_squared - sum_x ** 2) * (n * sum_y_squared - sum_y ** 2)) ** 0.5

    if denominator == 0:
      return None # To avoid division by zero

    return numerator / denominator


######## Calculate total mean phenotype: #################
# returns the total mean phenotype using space_dict and mut_dict
def calc_total_mean_phenotype(space_dict,mut_dict): 
    phenotype_list=[]
    for ind in space_dict: 
        phenotype=0 # reset phenotype
        for mut_list in space_dict[ind][2]: # list of list of mutations in each genome
            for mut in mut_list: # for each mutation in each list
                phenotype+= float(mut_dict[mut][0]) # add the mutation's value to the phenotype
        phenotype_list.append(phenotype)
    return np.mean(phenotype_list)

############## Calculate Zst: #################
# returns a vector of mean( Zs-Zt)^2 for each deme:     
#z_list = list of list of phenotypes, total_phenotype = global mean phenotype        
def calc_Zst(z_list,total_phenotype): 
    # for every deme: 
    all_difs=[]
    for deme_zs in z_list: 
        local_differences=[]
        # search for the individuals in that deme: 
        for phenotype in deme_zs: 
            local_differences.append((phenotype-total_phenotype)**2)
        all_difs.append(np.mean(local_differences))
    
    #print(len(all_difs))
    
    # return the list of deviations, and the mean deviation
    return all_difs 

# return the fitness of an interacting moth and plant  
# zm = moth phenotype, zp = plant phenotype, BM = fecundity benefit from mutualism
def calc_fitness(zm,zp,BM): 
    return BM*np.exp(((zm-zp)**2)*-1)


# an omnibus function that runs all of the phenotypic analysis 
def phenotypic_analysis(space_dict_m,mut_dict_m,space_dict_p,mut_dict_p,sigma,num): 
   # positions
    positions=random_sampler(space_dict_m,space_dict_p,sigma,num)
    # list of individuals in demes
    moth_demes=find_individuals(positions,space_dict_m,sigma)
    plant_demes=find_individuals(positions,space_dict_p,sigma)
    # list of deme phenotypes
    z_list_m=phenotypes(moth_demes,space_dict_m,mut_dict_m)
    z_list_p=phenotypes(plant_demes,space_dict_p,mut_dict_p)
    # total mean phenotypes
    z_t_m=calc_total_mean_phenotype(space_dict_m,mut_dict_m)
    z_t_p=calc_total_mean_phenotype(space_dict_p,mut_dict_p)
    # zsts: 
    zsts_m=calc_Zst(z_list_m,z_t_m)
    zsts_p=calc_Zst(z_list_p,z_t_p)

    # local means : 
    m_local_means=[np.mean(zs) for zs in z_list_m]
    p_local_means=[np.mean(zs) for zs in z_list_p]

    # spatial phenotypic correlation: 
    z_corr=calculate_correlation(m_local_means,p_local_means)

    # simulated local adaptation experiments: 
    home_list=[]
    global_list=[]
    foreign_list=[]
    for i,zms in enumerate(z_list_m): 
        for j,zps in enumerate(z_list_p): 
            for zm in zms: 
                for zp in zps: 
                    global_list.append(calc_fitness(zm,zp,40))
                    if i==j:
                        home_list.append(calc_fitness(zm,zp,40))
                    if i !=j: 
                        foreign_list.append(calc_fitness(zm,zp,40)) 
    
    # output summary stats: 
    output1=[len(space_dict_m),len(space_dict_p),len(mut_dict_m),len(mut_dict_p)]
    # output mean deme size: 
    output2= [np.mean([len(inds) for inds in moth_demes]),np.mean([len(inds) for inds in plant_demes]) ] 
    # output zsts and correlation: 
    output3=[np.mean(zsts_m),np.mean(zsts_p),z_corr]
    # output local adaptation results
    output4=[np.mean(home_list),np.mean(global_list),np.mean(foreign_list)]
    outs = [output1,output2,output3,output4]
    
    return outs,moth_demes,plant_demes

## an omnibus function that divides space into a bigxbig grid and calculates/records
## mutation effect, origin time, total frequency, spatial frequency, fst, and fis
def genomic_analysis(space_dict,mut_dict,big): 
    mut_list=[[],[],[],[],[],[]] # effect, time of origin,  total frequency, spatial range, fst,fis
    # A grid search for the spatial distribution of mutations
    ind_num=len(space_dict)
    neg_count=0
    for mut in mut_dict: 
        mut_list[0].append(abs(float(mut_dict[mut][0]))) # mutation effect
        mut_list[1].append(float(mut_dict[mut][1])) # origin time
        ft=float(mut_dict[mut][2])/(ind_num*2) # total mutation frequency
            
        mut_list[2].append(ft)
        # grid search for fst and fis: 
        f_s=[]
        f_is=[]
        delta = 10/big
        space_freq=0
        for xpos in range(1,big+1): 
            xpos=xpos/(big/10) # normalize to grid size
            for ypos in range(1,big+1): 
                ypos=ypos/(big/10) # normalize to grid size
                # set counters
                s_num=0# number of chromosomes 
                s_mut_num=0# number of chromosomes with mutation of interest
                flag=False # flag to track the first time a mutation occurs
                heterozygotes=0 # the number of homozygotes
                for ind in space_dict: 
                    # for every individual:
                    pos=space_dict[ind][0]
                    # if that individual is in the grid: 
                    if pos[0]<xpos and pos[0]>=xpos-delta and pos[1]<ypos and pos[1]>=ypos-delta: 
                        s_num+=2 # number of chromosomes
                         # number of chromosomes a mutation is on 
                        homozygote_flag=False # flag for homoozygote
                        # does an individual carry the relevant mutation?: 
                        local_num=0
                       
                        if mut in space_dict[ind][2][0]:
                            local_num+=1
                        if mut in space_dict[ind][2][1]:
                            local_num+=1
                        if flag==False and mut in space_dict[ind][2][0]: 
                            flag=True
                            space_freq+=1
                        if flag==False and mut in space_dict[ind][2][1]: 
                            flag=True
                            space_freq+=1
                        if local_num==1: 
                            heterozygotes+=1     
                        s_mut_num+=local_num
                        
                 
                if s_num >0 : # if a subpopulation had individuals, calcuate heterozygosity: 
                     # p*q
                    p=s_mut_num/s_num
                    q=1-p
                    f_12=heterozygotes/(s_num/2) # number of heterozygotes in subpop
                    
                    if p !=0 and p !=1: ##### if an allele is not absent or fixed
                        f_s.append(p*q)
                        #print(2*p*q)
                        #print(f_12)
                        #print(1-(f_12/(2*p*q)))
                        f_is.append(1-(f_12/(2*p*q)))
                        if p>1 or p <0: ## if this happened something went wrong
                            print('danger!')
                            #print(s_num)
                            #sys.exit()
                        
                    else: ## if an allele is fixed or absent: 
                        f_s.append(0)
                        f_is.append(1)           
        if ft>0 and ft<1: 
            fst = 1-(np.mean(f_s)/(ft*(1-ft)))
            fis=np.mean(f_is)
            #print(fis)
            if fst<0: 
                neg_count+=1
                fst=0
        else: 
            fst =0
        mut_list[3].append(space_freq)
        mut_list[4].append(fst)
        mut_list[5].append(fis)
    return mut_list

######## Function that formats output of genomic_analysis for heat maps
# returns 3 arrays one with mean, one with stdev, and one with number of mutations that fall in a bin
def new_heatmap_fn(num_bins,tot1,tot2,z_list,x_list,y_list,min1,min2):
    #print('indices',ind1,ind2)
    range1=tot1-min1
    range2=tot2-min2
    delta_f_x=range1/num_bins
    delta_f_y=range2/num_bins
    my_list2=[]
    for num1 in np.linspace(min1+delta_f_x,tot1,num_bins): 
        big_list=[]
        for num2 in np.linspace(min2+delta_f_y,tot2,num_bins): 
            big_list.append([])
        my_list2.append(big_list)
    
    # loop over the relevant lists: 
    index_count=-1
    for valx,valy in zip(x_list,y_list) :
        index_count+=1 
        for x,binx in enumerate(np.linspace(min1+delta_f_x,tot1,num_bins)):
            for y,biny in enumerate(np.linspace(min2+delta_f_y,tot2,num_bins)): 
                if valx<=binx and valx>binx-delta_f_x and valy<=biny and valy>biny-delta_f_y:
                    #if mut==mut2: print('made it! ')
                    my_list2[num_bins-y-1][x].append(abs(float(z_list[index_count])))
                    

    mean_plot2=[]    
    sigma_plot2=[]   
    num_plot=[]     
    for first in my_list2: 
        row=[]
        row_sigma=[]
        row_num=[]
        for second in first: 
            mean=0
            std=0
            if np.sum(second)!=0 and len(second)>1: 
                mean=np.mean(second)
                std =np.std(second)
            if np.sum(second)!=0 and len(second)==1: 
                mean=second[0]
                #std =np.std(second)
            numb=len(second)
            row_num.append(numb)
            row.append(mean)
            row_sigma.append(std)
        mean_plot2.append(row)
        sigma_plot2.append(row_sigma)
        num_plot.append(row_num)
    return mean_plot2,sigma_plot2,num_plot


## omnibus function that formats all output of genomic analysis for 
# heatmaps by applying new_heatmap_fn to the output of genomic_analysis
def array_format_genomic(mut_list,id_var1,id_var2,num):
    # mut_list has the form: [ [effect], [time of origin], [spatial range], [frequency],[fst]]
    r_max=max(mut_list[2])
    ft_max=max(mut_list[3])
    fst_max=max(mut_list[4])
    r_min=min(mut_list[2])
    ft_min=min(mut_list[3])
    fst_min=min(mut_list[4])
    #print(r_max,ft_max,fst_max,r_min,ft_min,fst_min)
    id_list=['tXr','fXr','tXfst']
    plot_list=['mean','sigma','num']
    mean_sigma_num_list_tr=new_heatmap_fn(num,11000,r_max,mut_list[0],mut_list[1],mut_list[2],0,r_min) # the portions of the mutation list to put on the x and y axis
    mean_sigma_num_list_fr=new_heatmap_fn(num,ft_max,r_max,mut_list[0],mut_list[3],mut_list[2],ft_min,r_min) # the portions of the mutation list to put on the x and y axis
    mean_sigma_num_list_tfst=new_heatmap_fn(num,11000,fst_max,mut_list[0],mut_list[1],mut_list[4],0,fst_min) # the portions of the mutation list to put on the x and y axis
    
    #print(mean_sigma_num_list_tfst[0])

    overall_list=[mean_sigma_num_list_tr,mean_sigma_num_list_fr,mean_sigma_num_list_tfst]

    for id,list_o_lists in zip(id_list,overall_list): 
        #print('should be 3',len(list_o_lists))
        for plot,f_list in zip(plot_list,list_o_lists): 
            #print('should be 5',len(f_list))
            #print('should be 5',len(f_list[0]))
            with open("{}_{}_space_{}_{}.csv".format(plot,id,str(id_var1),str(id_var2)), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_list)
            # with open("{}_{}_space_{}_{}.txt".format(plot,id,str(id_var1),str(id_var2)), 'w') as fp:
            #     json.dump(list, fp)

    max_params=['r_max',r_max,'ft_max',ft_max,'fst_max',fst_max,'r_min',r_min,'ft_min',ft_min,'fst_min',fst_min]
    return max_params



# function that returns the a dictionary of the spatial frequency of each 
# sampled mutation with the form: [mutation ID]: [list of frequencies]
# returns a second dictionary of mutations that occur in at least
# 10% of the demes
def make_mut_freq_dict(space_dict,demes,s_f): 
    # collect every mutation in every deme and the number of 
    # individuals in the deme
    deme_muts=[]
    size_list=[]
    muts_per_deme=[]
    for deme in demes:
        deme_muts_list=[]
        deme_size=0
        for ind in deme:
            deme_size+=1
            chromosomes=space_dict[ind][2]
            for chromosome in chromosomes: 
                for mut in chromosome:
                    deme_muts.append(mut)
                    deme_muts_list.append(mut)
        muts_per_deme.append(deme_muts_list)
        size_list.append(deme_size)
    #print(deme_muts)
    #sys.exit()
    mut_freq_dict={}
    f_mut_freq_dict={}
    deme_num=len(demes)
    # for every unique mutation record its frequency in every deme
    # record this in a dictionary with the form:
    ##  [mutation ID]: [list of mutation frequencies in demes]
    for mut in set(deme_muts): 
        space_num=0 # number of demes mut is in
        freq_list=[] # list of the mutation's frequencies
        for size,my_deme_muts in zip(size_list,muts_per_deme): 
            if mut in my_deme_muts: 
                space_num+=1 
                freq=my_deme_muts.count(mut)/(size*2) 
                freq_list.append(freq)
            else: 
                freq_list.append(0)
            
        mut_freq_dict[mut]=freq_list
        if space_num/len(demes)>=s_f: f_mut_freq_dict[mut]=freq_list
         
    return mut_freq_dict,f_mut_freq_dict

## spatial frequency of every mutation, input is a dictionary whose key is a mutation and 
## values are a list of spatial frequencies, samples is the number of demes, cutoff is a frequency
## cutoff that is user supplied. 
def s_afs(mut_freq_dict,samples,cutoff):
    afs=[]
    freq_mut_list=[]
    for mut in mut_freq_dict: 
        #print(mut_freq_dict[mut])
        occurrences=len([x for x in mut_freq_dict[mut] if x != 0])
        afs.append(occurrences)
        #print(occurrences)
        if occurrences/samples>=cutoff: 
            freq_mut_list.append(mut)
    return afs,freq_mut_list

# return a mutation's effect and origin time
def mutation_effects_time(mut_list,mut_dict): 
    effects=[]
    times=[]
    for mut in mut_list: 
        effects.append(abs(float(mut_dict[mut][0])))
        times.append(mut_dict[mut][1])
    return effects, times


### a function that breaks down ILD patterns: 
def ILD_stats(m_space_dict,p_space_dict,mut_dict_m,mut_dict_p,moth_demes,plant_demes,s_f): 
    # make dictionary of mutation frequency: 
    m_mut_freq_dict,m_f_mut_freq_dict=make_mut_freq_dict(m_space_dict,moth_demes,s_f)
    p_mut_freq_dict,p_f_mut_freq_dict=make_mut_freq_dict(p_space_dict,plant_demes,s_f)

    m_afs,m_freq_muts=s_afs(m_mut_freq_dict,len(moth_demes),s_f)
    p_afs,p_freq_muts=s_afs(p_mut_freq_dict,len(plant_demes),s_f)
   #print(m_freq_muts)
   #print(p_freq_muts)
    
    
    m_ef_freq,m_times_freq=mutation_effects_time(m_f_mut_freq_dict.keys(),mut_dict_m)
    p_ef_freq,p_times_freq=mutation_effects_time(p_f_mut_freq_dict.keys(),mut_dict_p)
    
    
    m_ef_tot,m_times_tot=mutation_effects_time(m_mut_freq_dict.keys(),mut_dict_m)
    p_ef_tot,p_times_tot=mutation_effects_time(p_mut_freq_dict.keys(),mut_dict_p)

    ILDS=[]
    freq_ILDS=[]
    for mutm,mutp in itertools.product(m_mut_freq_dict.keys(),p_mut_freq_dict.keys()):
        correlation = calculate_correlation(m_mut_freq_dict[mutm], p_mut_freq_dict[mutp])
        if correlation !=None:
            ILDS.append(correlation)
            if mutm in m_f_mut_freq_dict.keys() and mutp in p_f_mut_freq_dict.keys(): 
                freq_ILDS.append(correlation)

    mut_stuff=[m_ef_tot,p_ef_tot,m_ef_freq,p_ef_freq,m_afs,p_afs,m_times_tot,p_times_tot,m_times_freq,p_times_freq,ILDS,freq_ILDS]
    
    return mut_stuff

def pairwise_mean(list1,list2): 
    newlist=[]
    for l1,l2 in itertools.product(list1,list2): 
        newlist.append(np.mean(l1,l2))
    return newlist


################ put ILD stats in array format for heatmap #######################
## omnibus function that formats all output of ILD analysis for 
# heatmaps by applying new_heatmap_fn to the output of ILD_stats
# id_var1 is a text input of the dispersal rate
def array_format_ILD(mut_list,id_var1,num):
    ## take means of explanatory stats

    ef_mean=pairwise_mean(mut_list[0],mut_list[1])
    range_mean=pairwise_mean(mut_list[4],mut_list[5])
    times_mean=pairwise_mean(mut_list[6],mut_list[7])
    ILDS=mut_list[10]

    abs_ILDS=[abs(i) for i in ILDS]

    r_max=max(range_mean)
    t_max=max(times_mean)
    ef_max=max(ef_mean)
    ILD_max=max(ILDS)
    ILD_abs_max=max(abs_ILDS)

    r_min=min(range_mean)
    t_min=min(times_mean)
    ef_min=min(ef_mean)
    ILD_min=min(ILDS)
    ILD_abs_min=min(abs_ILDS)
    
    id_list=['tXrXILD','txexild','txildxe','txildaxe']
    plot_list=['mean','sigma','num']
    m_s_n_txrxild=new_heatmap_fn(num,11000,r_max,abs_ILDS,times_mean,range_mean,0,r_min) # the portions of the mutation list to put on the x and y axis
    m_s_n_list_txexild=new_heatmap_fn(num,11000,ef_max,abs_ILDS,times_mean,range_mean,0,ef_min) # the portions of the mutation list to put on the x and y axis
    m_s_n_txildxe=new_heatmap_fn(num,11000,ILD_max,ef_mean,times_mean,ILDS,0,ILD_min) # the portions of the mutation list to put on the x and y axis
    m_s_n_txildaxe=new_heatmap_fn(num,11000,ILD_abs_max,ef_mean,times_mean,abs_ILDS,0,0)
    #print(mean_sigma_num_list_tfst[0])

    overall_list=[m_s_n_txrxild,m_s_n_list_txexild,m_s_n_txildxe,m_s_n_txildaxe]

    for id,list_o_lists in zip(id_list,overall_list): 
        #print('should be 3',len(list_o_lists))
        for plot,f_list in zip(plot_list,list_o_lists): 
            #print('should be 5',len(f_list))
            #print('should be 5',len(f_list[0]))
            with open("{}_{}_space_ILD_{}.csv".format(plot,id,str(id_var1)), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_list)
            

    max_params=['r_max',r_max,'ef_max',ef_max,'ILD_max',ILD_max,'ILD_abs_max',ILD_abs_max,'r_min',r_min,'ef_min',ef_min,'ILD_min',ILD_min,'ILD_abs_min',ILD_abs_min]
    return max_params

######### a function that sorts lists #############
def sort_lists(list1,list2): 
    combined=zip(list1,list2)
    # Sort the combined list (it will sort by the first element of each tuple)
    sorted_combined = sorted(combined,reverse=True)
    sorted_list1, sorted_list2 = zip(*sorted_combined)
    return sorted_list1, sorted_list2

############ a function that filters a list to only have values that fall within 
# a set interval  effects = the list, interval = interval      ######## 

def filter_interval(effects,interval):
    labels=[]
    marker=True
    j=0
    for i in effects: 
        if marker==True: 
            labels.append(i)
            marker=False
        else: 
            labels.append("")
            j+=1
            if j==interval:
                marker=True
                j=0
    return labels

############# de-trended ILD heatmaps for a single replicate ###################
############## s_f is a frequency cutoff for the least frequent mutation's spatial prevalence #
def ILD_full_array(m_space_dict,p_space_dict,mut_dict_m,mut_dict_p,moth_demes,plant_demes,s_f,sigma):
    # make dictionary of mutation frequency: 
    m_mut_freq_dict,m_f_mut_freq_dict=make_mut_freq_dict(m_space_dict,moth_demes,s_f)
    p_mut_freq_dict,p_f_mut_freq_dict=make_mut_freq_dict(p_space_dict,plant_demes,s_f)

    m_keys=[key for key in m_mut_freq_dict.keys()]
    p_keys=[key for key in p_mut_freq_dict.keys()]

    m_effects = [abs(float(mut_dict_m[key][0])) for key in m_mut_freq_dict.keys()]
    p_effects = [abs(float(mut_dict_p[key][0])) for key in p_mut_freq_dict.keys()]

    # Zip them together, with the sorting criteria first
    sorted_effects_m, sorted_m_keys = sort_lists(m_effects, m_keys)
    # Unzip them back into two sorted lists
    sorted_effects_p, sorted_p_keys = sort_lists(p_effects, p_keys)

    ILDS=[]
    ILD_matrix=[]
  
    for mutm in sorted_m_keys:
        row=[]
        for mutp in sorted_p_keys:
            correlation = calculate_correlation(m_mut_freq_dict[mutm], p_mut_freq_dict[mutp])
            if correlation !=None:
                row.append(abs(correlation)) 
            else: 
                row.append(0)
        ILD_matrix.append(row)
   
    sorted_effects_m_r=[round(i,2)for i in sorted_effects_m ] 
    ylabels=filter_interval(sorted_effects_m_r,len(sorted_effects_m_r)//10)
    sorted_effects_p_r=[round(i,2)for i in sorted_effects_p ] 
    xlabels=filter_interval(sorted_effects_p_r,len(sorted_effects_m_r)//10)

    with open("ILD_full_effects_{}.csv".format(sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(ILD_matrix)
    
    with open("ILD_full_effects_xlab_{}.csv".format(sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(xlabels)
    
    with open("ILD_full_effects_ylab_{}.csv".format(sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(ylabels)


    fkeysm=[key for key in m_f_mut_freq_dict.keys()]
    fkeysp=[key for key in p_f_mut_freq_dict.keys()]

    f_m_effects = [abs(float(mut_dict_m[key][0])) for key in fkeysm]
    f_p_effects = [abs(float(mut_dict_p[key][0])) for key in fkeysp]

    # Zip them together, with the sorting criteria first
    f_sorted_effects_m, f_sorted_m_keys = sort_lists(f_m_effects, fkeysm)
    # Unzip them back into two sorted lists
    f_sorted_effects_p, f_sorted_p_keys = sort_lists(f_p_effects, fkeysp)

    f_ILD_matrix=[]
    
    for mutm in f_sorted_m_keys:
        row=[]
        for mutp in f_sorted_p_keys:
            correlation = calculate_correlation(m_mut_freq_dict[mutm], p_mut_freq_dict[mutp])
            if correlation !=None:
                row.append(abs(correlation))
            else: 
                row.append(0)
        f_ILD_matrix.append(row)

    sorted_effects_m_r=[round(i,2)for i in f_sorted_effects_m ] 
    f_ylabels=filter_interval(sorted_effects_m_r,len(sorted_effects_m_r)//10)

    sorted_effects_p_r=[round(i,2)for i in f_sorted_effects_p ] 
    f_xlabels=filter_interval(sorted_effects_p_r,len(sorted_effects_m_r)//10)

    with open("ILD_f_{}_full_effects_{}.csv".format(str(s_f), sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_ILD_matrix)
    
    with open("ILD_f_{}_effects_xlab_{}.csv".format(str(s_f),sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_xlabels)
    
    with open("ILD_f_{}_effects_ylab_{}.csv".format(str(s_f),sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_ylabels)

    f_m_times = [abs(float(mut_dict_m[key][1])) for key in fkeysm]
    f_p_times = [abs(float(mut_dict_p[key][1])) for key in fkeysp]

    # Zip them together, with the sorting criteria first
    f_sorted_times_m, f_sorted_m_keys = sort_lists(f_m_times, fkeysm)
    # Unzip them back into two sorted lists
    f_sorted_times_p, f_sorted_p_keys = sort_lists(f_p_times, fkeysp)

    f_ILD_matrix=[]
    for mutm in f_sorted_m_keys:
        row=[]
        for mutp in f_sorted_p_keys:
            correlation = calculate_correlation(m_mut_freq_dict[mutm], p_mut_freq_dict[mutp])
            if correlation !=None:
                row.append(abs(correlation))
            else: 
                row.append(0)
        f_ILD_matrix.append(row)

    sorted_times_m_r=[round(i,2)for i in f_sorted_times_m ] 
    ylabels=filter_interval(sorted_times_m_r,len(sorted_times_m_r)//10)

    sorted_times_p_r=[round(i,2)for i in f_sorted_times_p ] 
    xlabels=filter_interval(sorted_times_p_r,len(sorted_times_p_r)//10)

    with open("ILD_f_{}_full_times_{}.csv".format(str(s_f), sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(f_ILD_matrix)
    
    with open("ILD_f_{}_times_xlab_{}.csv".format(str(s_f),sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(xlabels)
    
    with open("ILD_f_{}_times_ylab_{}.csv".format(str(s_f),sigma), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(ylabels)


########## hybridization experiments #################
# impute parent id1, parent id2, and recombination rate
# returns list of offspring phenotypes and dictionary of offspring with form: 
# ID number : [[mutations on chromosome 1],[mutations on chromosome 2]]
def make_offspring(p1,p2,recomb_rate,genome_size,num_offspring,space_dict,mut_dict,F2,offspring_g_dict_F1,offspring_count):
    p= recomb_rate*genome_size
    recomb1=[np.random.binomial(genome_size,recomb_rate) for i in range(num_offspring)]
    recomb2=[np.random.binomial(genome_size,recomb_rate) for i in range(num_offspring)] 
    
    offspring_zs=[]
    offspring_g_dict={}
    #for i in range(num_offspring): offspring_g_dict[str(i)]=[]
    for offspring in range(num_offspring): 
        c_list=[]
        if True: 
            breakpoints1=[random.randint(0, genome_size) for i in range(recomb1[offspring])]
            breakpoints1.sort()
            
            c_choice=np.random.binomial(1,0.5)
            if F2==False: 
                muts_c=space_dict[p1][2][c_choice]
                muts_nc=space_dict[p1][2][1-c_choice]
            if F2==True: 
                muts_c=offspring_g_dict_F1[p1][c_choice]
                muts_nc=offspring_g_dict_F1[p1][1-c_choice]
                


            mut_list=[]
            if recomb1[offspring]==0: 
                c_list.append(muts_c)
            else: 
                for i,point in enumerate(breakpoints1): 
                    if i ==0: 
                        for mut in muts_c: 
                            if float(mut_dict[mut][3])<=point: 
                            
                                mut_list.append(mut)
                        if breakpoints1[-1]==point: 
                            for mut in muts_nc: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)
                    
                    if i >0 and i%2==1: 

                        for mut in muts_nc: 
                            position=float(mut_dict[mut][3])
                            if position>breakpoints1[i-1] and position<=point: 
                                mut_list.append(mut)
                        
                        if breakpoints1[-1]==point: 
                            for mut in muts_c: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)
                    
                    if i >0 and i%2==0: 
                        for mut in muts_c: 
                            position=float(mut_dict[mut][3])
                            if position>breakpoints1[i-1] and position<=point: 
                                mut_list.append(mut)
                        
                        if breakpoints1[-1]==point: 
                            for mut in muts_c: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)

                    
            
                c_list.append(mut_list)

            breakpoints2=[random.randint(0, genome_size) for i in range(recomb2[offspring])]
            breakpoints2.sort()
            c_choice=np.random.binomial(1,0.5)
            
            muts_c=space_dict[p2][2][c_choice]
            muts_nc=space_dict[p2][2][1-c_choice]
            mut_list=[]

            mut_list=[]
            if recomb2[offspring]==0: 
                c_list.append(muts_c) 
            else: 
                for i,point in enumerate(breakpoints2): 
                    if i ==0: 
                        for mut in muts_c: 
                            if float(mut_dict[mut][3])<=point: 
                            
                                mut_list.append(mut)
                        if breakpoints2[-1]==point: 
                            for mut in muts_nc: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)
                    
                    if i >0 and i%2==1: 

                        for mut in muts_nc: 
                            position=float(mut_dict[mut][3])
                            if position>breakpoints2[i-1] and position<=point: 
                                mut_list.append(mut)
                        
                        if breakpoints2[-1]==point: 
                            for mut in muts_c: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)
                    
                    if i >0 and i%2==0: 
                        for mut in muts_c: 
                            position=float(mut_dict[mut][3])
                            if position>breakpoints2[i-1] and position<=point: 
                                mut_list.append(mut)
                        
                        if breakpoints2[-1]==point: 
                            for mut in muts_c: 
                                if float(mut_dict[mut][3])>point: 
                                    mut_list.append(mut)

                c_list.append(mut_list)
                
            # recombination in parent 1: 
            #breakpoints=[random.randint(0, genome_size) for i in range(recomb[offspring])]
            #chrom_start=np.random.binomial(num_offspring,0.5)

        phenotype=0
        for mut_list in c_list:
            for mut in mut_list: 
                phenotype+= float(mut_dict[mut][0])

        offspring_zs.append(phenotype)
       
        offspring_g_dict[str(offspring_count)]=c_list
        offspring_count+=1
    return offspring_zs,offspring_g_dict


## sample a focal deme and a twin. num=number of twin demes, factor = number of sigmas 
## phenotypic twins must be from a focal deme. 
def sample_similar_demes(space_dict_m,mut_dict_m,space_dict_p,mut_dict_p,num,sigma,factor): 
    # list of xy positions to sample
    positions=random_sampler(space_dict_m,space_dict_p,sigma,num)
   # list of list of individuals in each deme [[inds in deme 1],[inds in deme 2]]
    deme_list_m=find_individuals(positions,space_dict_m,sigma)
    deme_list_p=find_individuals(positions,space_dict_p,sigma)
# list of list of phenotypes in each deme [[phenotypes in deme 1],[phenotypes in deme 2]]
    z_list_m= phenotypes(deme_list_m,space_dict_m,mut_dict_m)
    z_list_p= phenotypes(deme_list_p,space_dict_p,mut_dict_p)
    
    # find a phenotypic 'twin' to a focal deme: 
    twins=[]
    twin_zs=[]
    keys_m=[key for key in space_dict_m.keys()]
    for i,deme_z in enumerate(z_list_m): 
        z_bar_m=np.mean(deme_z)
        z_bar_p=np.mean(z_list_p[i])
        local_w=calc_fitness(z_bar_m,z_bar_p,40)
        # mean fitness of locals
        search=True
        
        while search == True: 
            choice=np.random.choice(keys_m,1,replace=False)
            phenotype=0
            for mut_list in space_dict_m[choice[0]][2]: 
                for mut in mut_list: 
                    phenotype+= float(mut_dict_m[mut][0])
            
            if calc_fitness(phenotype,z_bar_p,40) - local_w >= -0.5: 
                x,y=space_dict_m[choice[0]][0]

                # make sure we're not too close to a spatial boundary b/c I don't want to deal with periodicity
                if x>0+sigma and x<10-sigma and y>0+sigma and y<10-sigma: 
                    
                    # make sure we're nott too close to the focal subpopulation: 
                    if (x-positions[i][0])**2 + (y - positions[i][1])**2 > factor*sigma:
                        #print(x,y,positions[i])
                        
                        for plant in space_dict_p: 
                        # make sure there will be a partner: 
                            if (space_dict_p[plant][0][0]-x)**2+ (space_dict_p[plant][0][1]-y)**2<=sigma**2: 
                                search=False
                                twins.append(choice[0])
                                twin_zs.append(phenotype)
                                break
            
    
    F1_deme_zs=[]
    F2_deme_zs=[]
    for focal,twin in zip(deme_list_m,twins): 
        #parents=random.choices(focal, k=5)
        offspring_count=0
        offspring_g_dict={}
        F1_z_list=[]
        for parent in focal: 
            z_list_F1,offspring_g_dict_F1=make_offspring(parent,twin,10**-5,2*10**5,5,space_dict_m,mut_dict_m,False,[],offspring_count)
            offspring_g_dict.update(offspring_g_dict_F1)
            for z in z_list_F1: F1_z_list.append(z)
        F1_deme_zs.append(F1_z_list)
        offspring_ids=list(offspring_g_dict.keys())
        offspring_count=0
        F2_z_list=[]
        for parent in focal: 
            for offspring in offspring_ids: 
                z_list_F2,offspring_g_dict_F2=make_offspring(offspring,parent,10**-5,2*10**5,5,space_dict_m,mut_dict_m,True,offspring_g_dict,offspring_count)
                for z in z_list_F2: F2_z_list.append(z)
        F2_deme_zs.append(F2_z_list)

    

    moth_native_w=[]
    for deme_m,deme_p in zip(z_list_m,z_list_p): 
        for z_m in deme_m: 
            for z_p in deme_p: 
                moth_native_w.append(calc_fitness(z_m,z_p,40))


    moth_invader_w=[]
    for twin,deme_p in zip(twin_zs,z_list_p): 
        for z_p in deme_p: 
            moth_invader_w.append(calc_fitness(twin,z_p,40))

    
    F1_fitness=[]
    for deme_m,deme_p in zip(F1_deme_zs,z_list_p): 
        for z_m in deme_m: 
            for z_p in deme_p: 
                F1_fitness.append(calc_fitness(z_m,z_p,40))
    

    F2_fitness=[]
    for deme_m,deme_p in zip(F2_deme_zs,z_list_p): 
        for z_m in deme_m: 
            for z_p in deme_p: 
                F2_fitness.append(calc_fitness(z_m,z_p,40))

    return [np.mean(moth_native_w),np.mean(moth_invader_w),np.mean(F1_fitness),np.mean(F2_fitness)]


####### write a .csv file with a text name, and values for each row

def write_my_csv(name,values):
    with open(name, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerows(values)

def run_all_analysis():
    wd = os.getcwd()
    files=os.listdir(wd)
    pdict={}
    mdict={}
    ds=[]
    
    for file in files: 
        if file[-11:] == "11000_o.txt":
            vals=[float(s) for s in re.findall(r"-?\d+\.?\d*", file)]
            if file[0]=='p': 
                #print(vals[0:3])
                pdict[file]=[vals[0:3]]
            if file[0]=='m': 
                #print(vals[0:3])
                mdict[file]=[vals[0:3]]
                ds.append(vals[0])
                
                

    # list of means of genomic stats
    genom_means=[]
    ## dictionary for genomic stats per mutation: 
    dss=set(ds)
    d_dict_m={}
   
    for d in dss: 
        d_dict_m[str(d)]=[[],[],[],[],[],[]] # effect, time of origin,  total frequency, spatial range, fst,fis

    ef_dict={}
    # for d in dss: 
    #     for rep in range(1,11): 
    #         val=str(d)+"_"+str(rep)
    #         ef_dict[val]=[]


    ## phenotypic analysis list:
    z_dict={}
    for d in dss: 
        z_dict[str(d)]=[[],[],[],[]]

    ### dictionary for ILD analysis
    ILD_dict={}
    for d in dss: 
        ILD_dict[str(d)]=[[],[],[],[],[],[],[],[],[],[],[],[]]
    
    # list for twin study
    twin_list=[]

    # ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    # pool = mp.Pool(processes=ncpus)

    # data = itertools.product(mdict,pdict)
    genom=False
    z=True
    z_save=False
    ILD=False
    Twins=True
    for m,p in itertools.product(mdict,pdict): 
        if mdict[m]==pdict[p]:
            sigma_p,sigma_m,rep=pdict[p][0]
            space_dict_m,mut_dict_m = import_output_full("./{}".format(m))
            print(sigma_m)
            if genom==True and rep <=5: ## toggle genomic analysis
                num=20
                mut_list=genomic_analysis(space_dict_m,mut_dict_m,num) # effect, time of origin,  total frequency, spatial range, fst,fis
                
                # append summary stats to a master list: 
                for i,list in enumerate(mut_list): 
                    target=d_dict_m[str(sigma_m)][i]
                    for j in list: target.append(j)
                if rep<=3: 
                    effects=mut_list[0]
                    key=str(sigma_m)+"_"+str(rep)
                    ef_dict[key].append(effects)

                
                ### compute the means for each summary stat
                mean_muts=[np.mean(sublist) for sublist in mut_list]
                mean_muts.append(sigma_m)
                genom_means.append(mean_muts)
            print('getting closer! ')
            if z==True:  
                ## get plant data and run phenotypic analysis: 
                space_dict_p,mut_dict_p= import_output_full("./{}".format(p))
                num=200
                output,moth_demes,plant_demes=phenotypic_analysis(space_dict_m,mut_dict_m,space_dict_p,mut_dict_p,float(sigma_m),num)
                if z_save==True: 
                    for i,list in enumerate(output):
                        target=z_dict[str(sigma_m)][i]
                        for j in list: target.append(j)
            
            print('z analysis done! ')
            if ILD==True: # and rep <=5
                mut_stuff=ILD_stats(space_dict_m,space_dict_p,mut_dict_m,mut_dict_p,moth_demes,plant_demes,0.1)
                # ^ compute stats for ILD
                for i,list in enumerate(mut_stuff):
                    target=ILD_dict[str(sigma_m)][i]
                    for j in list: 
                        target.append(j)
                print("ILD part 1")
                if float(rep)==1:
                    ILD_full_array(space_dict_m,space_dict_p,mut_dict_m,mut_dict_p,moth_demes,plant_demes,0.1,sigma_m)
                print('hooray! ')
            if Twins==True and rep<=5: 
            
                num=20
                # conduct twin study: 
                output=sample_similar_demes(space_dict_m,mut_dict_m,space_dict_p,mut_dict_p,num,sigma_m,10)

                twin_list.append(output)

            
    
    if z_save==True: 
    ## save all dictionaries:
        with open("z_dict.json", "w") as file:
                json.dump(z_dict, file)
    if genom==True: 
        with open("d_dict_m.json", "w") as file:
                json.dump(d_dict_m, file)
    if ILD==True: 
        with open("ILD_dict.json", "w") as file:
                json.dump(ILD_dict, file)
    if genom==True: 
        with open("ef_dict.json", "w") as file:
                json.dump(ef_dict, file)


    ### format outputs with all mutations to make heatmaps
    g_param_list=[]
    ILD_param_list=[]
    for d in dss: 
        num=22
        g_params=array_format_genomic(d_dict_m[str(d)],str(d),'m',num)
        g_param_list.append(g_params)
        ILD_params=array_format_ILD(ILD_dict[str(d)],str(d),num)
        ILD_param_list.append(ILD_params)
    

    write_my_csv('g_params.csv',g_param_list)
    write_my_csv('ILD_params.csv',ILD_param_list)
    write_my_csv('twin_study_results.csv',twin_list)

if __name__ == "__main__":
    run_all_analysis()


     
    

    
    