# GWAS-BCCAD



# Download the Data Files

Please download the image features (global and local connectivity metrics) in the [Data folder](https://github.com/elssam/GWAS-BCCAD/tree/master/Data)

# Feature Preparing

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:51:33 2018

@author: samar
             Data Analysis of MRI Data and Gene Expression

This script is meant to prepare the data for analysis. It:
    1- Import the individual files and get their global features,
    2- Do some descriptive statistics and plots
    3- Import the genes data and integrate them with the brain features

Input:
    Individual files (global), gene expression files.
Output:
    Save two datasets and important figures

"""
#############################################################################################################
### Import modules/Libraries

import pandas 
import matplotlib.pylab as plt
import seaborn
#import scipy.stats as stats
import glob

### Functions:
def combine(directory,visit):
    ''' This function combines csv files of MRI features in a directory - files should start with 0 (or 1) - directory should be the name of the directory you sotred your data; which normaly contains the subjects features of a specific visit and for either case/controls'''
    infiles=glob.glob(str(directory)+'/*'+str(visit)+"*")#glob.glob is similar to "ls" in linux command line, it saves all the files in current directory in a list
    outfile=pandas.DataFrame()
    for i in infiles:
        infile=pandas.read_csv(i)
        outfile=outfile.append(infile,ignore_index=True)
    return(outfile)
def add_column(dataframe,column_name,column_values):
    '''This function appends a column with a fix values at the end of it. Pass as parameter: dataframe,column_name,column_value'''
    dataframe_result=dataframe[:]
    dataframe_result[str(column_name)]=pandas.Series([str(column_values)]*len(dataframe))
    return(dataframe_result)
def append_both(dataframe1,dataframe2):
    '''This function takes two dataframes and append them into one'''
    return(dataframe1.append(dataframe2,ignore_index=True))
def calculate_diff(dataframe,column_name):
    '''This function calculates the difference between the metric in baseline and in follow-up'''
    before=" "+str(column_name)+"_"+"x"
    after=" "+str(column_name)+"_"+"y"
    dataframe[str(column_name)]=abs(dataframe[before]-dataframe[after])
    return(dataframe)
#############################################################################################################
### Combine the AD and Control subjects as followup/basline 
*Folders names:*

case='awdataontheway/AD_FA_02_binary'
control='awdataontheway/CONTROL_FA_02_binary'
print("Starting to work on the local features...")
AD_infile_followup = combine(case,"followup")
AD_infile_baseline = combine(case,"baseline")
Control_infile_followup=combine(control,'followup')
Control_infile_baseline=combine(control,'baseline')
#check
print("Four dataframes are ready, each corresponds to a StatusvsVisit.\n\nChecking...")
print(AD_infile_followup.head(n=2))
print(Control_infile_baseline.tail(n=2))
print(AD_infile_baseline.columns)
#Creating the Visit column
AD_infile_followup = add_column(dataframe=AD_infile_followup,column_name="Visit",column_values="FOLLOWUP")
AD_infile_baseline = add_column(AD_infile_baseline,column_name="Visit",column_values="BASELINE")
Control_infile_followup=add_column(Control_infile_followup,column_name="Visit",column_values="FOLLOWUP")
Control_infile_baseline=add_column(Control_infile_baseline,column_name="Visit",column_values="BASELINE")
#There will be an error message. Ignore it.
#check
print(AD_infile_baseline.columns)
print(Control_infile_followup.head(n=2))
print(AD_infile_baseline.head(n=2))
#Descriptive (sample size)
print("\n\nTotal of AD subjects at follow-up is",len(AD_infile_followup))
print("Total of AD subjects at baseline is",len(AD_infile_baseline))
print("Total of Control subjects at follow-up is",len(Control_infile_followup))
print("Total of Control subjects at baseline is",len(Control_infile_baseline))
#Combining cases together, and controls together
AD_infile = append_both(AD_infile_followup,AD_infile_baseline)
Control_infile=append_both(Control_infile_followup,Control_infile_baseline)
#check
print("Cases combined, and controls combined check:")
print(Control_infile.tail(n=2))
print(AD_infile.head(n=2))
#Descriptive to chek the previous sample size (supposed to be double the above)
print("Total of Controls, for both visits is",len(Control_infile))
print("Total of AD, for both visits is",len(AD_infile))
#Creating the Status column
AD_infile= add_column(AD_infile,column_name="Status",column_values="AD")
Control_infile= add_column(Control_infile,column_name="Status",column_values="CONTROL")
#check
print("Column Status been added")
print(AD_infile.head(n=3))
print(Control_infile.tail(n=3))
print(AD_infile.columns)
print(Control_infile.columns)
#Combine both above to one MRI dataframe
DTI_infile=append_both(AD_infile,Control_infile)
#check
print("Checking the overall local features (still individuals are duplicated, everyone got two values, in follow-up, and baseline")
print(DTI_infile.tail(n=3))
print(DTI_infile.head(n=3))
print(DTI_infile.columns)
#If there is 'inf'
#MRI_infile = MRI_infile.replace(to_replace='inf',value='NaN' ,regex=True)
#Save/Import to csv file
DTI_infile.to_csv("Local_Features_Duplicated.csv")
#Reshaping the data
Reshapred_Data = pandas.DataFrame() #It was Combined_Data
print("Reshaping the data to aggregate Status together, having two columns for baseline and follow-up:")
CONTROLS= pandas.merge(Control_infile_baseline,Control_infile_followup,on="Subject")
CONTROLS=add_column(dataframe=CONTROLS,column_name="Status",column_values="CONTROL")
CASES= pandas.merge(AD_infile_baseline,AD_infile_followup,on="Subject")
CASES=add_column(dataframe=CASES,column_name="Status",column_values="AD")
print("Done reshaing.\n\nNew columns of the reshaped data are:\n\n",)
Reshapred_Data = CASES.append(CONTROLS,ignore_index=True)
Reshapred_Data.to_csv("Local_Features_Reshapred.csv")
#############################################################################################################
#Calculating the differencecs
NEW=pandas.read_csv("Local_Features_Reshapred.csv")
print(NEW.head())
calculate_diff(dataframe=NEW,column_name="transitivity ")
calculate_diff(dataframe=NEW,column_name="global_eff")
calculate_diff(dataframe=NEW,column_name="louvain")
calculate_diff(dataframe=NEW,column_name="char_path_len")
print(NEW.columns)
############NEW.to_csv("Local_Features_Differences.csv")
#############################################################################################################
#Plotting: Descriptive Statistics
#Plots of Featurs differences
#seaborn.set(rc={'figure.figsize':(11.7,8.27),'axes.facecolor':'white','figure.facecolor':'white'}) # 'figure.facecolor':'white'})

#seaborn.set(rc={'figure.figsize':(3,3)})
seaborn.set_style('ticks')
seaborn.pairplot(NEW.loc[:,'transitivity ':'char_path_len'],kind="reg")
plt.savefig("Distribution_Differences.png")
plt.clf()
seaborn.set_style('ticks')
seaborn.set(rc={'figure.figsize':(1,1)})
#Plot features with Status
seaborn.pairplot(NEW.loc[:,'Status':'char_path_len'],hue="Status", markers=["o", "D"])
plt.savefig("Distribution_Differences_Status.png")
plt.clf()

#Plot features with Visit
seaborn.pairplot(DTI_infile,hue="Visit", palette="Set1", markers=["o", "D"])
plt.savefig("Distribution_Differences_Visit.png")
#Box plots
plt.clf()
seaborn.set_style('ticks')
seaborn.set(rc={'figure.figsize':(5,5)})
seaborn.boxplot(data=NEW.loc[:,'transitivity ':'char_path_len'],palette="RdBu")
plt.savefig("Boxplot_Differences.png")
#Plotting two categorical in a boxplot    
plt.clf()

seaborn.set(rc={'figure.figsize':(5,5)})
for y in [" louvain"," global_eff"," transitivity "," char_path_len"]:
    plt.clf()
    seaborn.set_style('ticks')
    z=seaborn.boxplot(x="Status", y=y, hue="Visit",data=DTI_infile,palette="Set3",order=["CONTROL","AD"], hue_order=["BASELINE","FOLLOWUP"])#,size=(4,4))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if y == " transitivity ":
        y=y[:-1]
    plt.savefig("Boxplot_Differences_Status&Visit"+y[1:]+".png") #avoid the space
    
#I added this to save them individually
y =" louvain"
seaborn.set(rc={'figure.figsize':(5,5)})
seaborn.set_style('ticks')
#seaborn.set(rc={'figure.figsize':(5,5)})
z=seaborn.boxplot(x="Status", y=y, hue="Visit",data=DTI_infile,palette="Set3",order=["CONTROL","AD"], hue_order=["BASELINE","FOLLOWUP"])#,size=(4,4))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


#############################################################################################################
#Data Description
print(NEW.describe())
#############################################################################################################
#Adding MCI
#Here you do need the first lines of this proram.
mci='awdataontheway/EMCI_FA0.2_binary'
print("Starting to work on joining the MCI data...")
MCI_infile_followup = combine(mci,"followup")
MCI_infile_baseline = combine(mci,"baseline")
#check
print("Two dataframes are ready, each corresponds to a MCI Visit.\n\nChecking...")
print(MCI_infile_followup.head(n=2))
print(MCI_infile_baseline.columns)
#Creating the Visit column
MCI_infile_followup = add_column(dataframe=MCI_infile_followup,column_name="Visit",column_values="FOLLOWUP")
MCI_infile_baseline = add_column(MCI_infile_baseline,column_name="Visit",column_values="BASELINE")
#There will be an error message. Ignore it.
#check
print(MCI_infile_baseline.columns)
print(MCI_infile_followup.head(n=2))
#Descriptive (sample size)
print("\n\nTotal of MCI subjects at follow-up is",len(MCI_infile_followup))
print("Total of MCI subjects at baseline is",len(MCI_infile_baseline))
#Combining cases together, and controls together
MCI_infile=append_both(MCI_infile_followup,MCI_infile_baseline)
#check
print("MCI combined check:")
print(MCI_infile.tail(n=2))
#Descriptive to chek the previous sample size (supposed to be double the above)
print("Total of MCI, for both visits is",len(MCI_infile))
#Creating the Status column
MCI_infile= add_column(MCI_infile,column_name="Status",column_values="MCI")
#check
print("Column Status been added")
print(MCI_infile.head(n=3))
print(MCI_infile.columns)
#Combine both above to one MRI dataframe
DTI2_infile=append_both(DTI_infile,MCI_infile)
#check
print("Checking the overall local features (still individuals are duplicated, everyone got two values, in follow-up, and baseline")
print(DTI2_infile.tail(n=3))
print(DTI2_infile.head(n=3))
print(DTI2_infile.columns)
#Save/Import to csv file
#############DTI2_infile.to_csv("Local_Features_Duplicated_withMCI.csv")
#############################################################################################################
#Reshaping the data
print("Reshaping the data with MCI to aggregate Status together, having two columns for baseline and follow-up:")
MCIS= pandas.merge(MCI_infile_baseline,MCI_infile_followup,on="Subject")
MCIS=add_column(dataframe=MCIS,column_name="Status",column_values="MCI")
Reshapred_Data2 = Reshapred_Data.append(MCIS,ignore_index=True)
#########Reshapred_Data2.to_csv("Local_Features_Reshapred_withMCI.csv")
#############################################################################################################
#Calculating the differencecs
NEW2=pandas.read_csv("Local_Features_Reshapred_withMCI.csv")
print(NEW2.head())
calculate_diff(dataframe=NEW2,column_name="transitivity ")
calculate_diff(dataframe=NEW2,column_name="global_eff")
calculate_diff(dataframe=NEW2,column_name="louvain")
calculate_diff(dataframe=NEW2,column_name="char_path_len")
print(NEW2.columns)
##################NEW2.to_csv("Local_Features_Differences_withMCI.csv")
#############################################################################################################


import scipy.stats as stats
x=['Status', 'transitivity ', 'global_eff', 'louvain', 'char_path_len']
#       'LPR1_2', 'SPON1', 'SPON1_2', 'SPON1_3', 'SPON1_4','louvain_x','louvain_y','louvain','global_eff_x','global_eff_y','global_eff', 'transitivity_x', 'transitivity_y', 'transitivity']
CCASE=NEW.Status=="AD"
CCONTROL=NEW.Status=='CONTROL'
for i in x:
 result = stats.mannwhitneyu(NEW[i][CCASE],NEW[i][CCONTROL])#lternative='greater')
#    #if result.pvalue < 0.05:
 print(i,":\t",result)#spearmanr(NEW[i],NEW[j]))
#
#TTEST
X=[' char_path_len_x',' global_eff_x', ' transitivity _x', ' louvain_y',' char_path_len_y', ' global_eff_y', ' transitivity _y']
for i in X:#x[1:]:#X:
    result = stats.ttest_ind(NEW[i][CCASE],NEW[i][CCONTROL],equal_var=False)
#    #if result.pvalue < 0.05:
    print(i,":\t",result)#spearmanr(NEW[i],NEW[j]))
#
## In[ ]:
print('char')
print(stats.ttest_rel(NEW[' char_path_len_x'], NEW[' char_path_len_y'], axis=0, nan_policy='propagate'))
print('global')
print(stats.ttest_rel(NEW[' global_eff_x'], NEW[' global_eff_y'], axis=0, nan_policy='propagate'))
print('trans')
print(stats.ttest_rel(NEW[' transitivity _x'], NEW[' transitivity _y'], axis=0, nan_policy='propagate'))
print('louvain')
print(stats.ttest_rel(NEW[' louvain_x'], NEW[' louvain_y'], axis=0, nan_policy='propagate'))
#

#WIlcoxon
print('char')
print(stats.wilcoxon(NEW[' char_path_len_x'], NEW[' char_path_len_y'], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW[' global_eff_x'], NEW[' global_eff_y'], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW[' transitivity _x'], NEW[' transitivity _y'],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW[' louvain_x'], NEW[' louvain_y'], zero_method='wilcox', correction=False))
#same but separately for cases, controls and mci




#For MCI
NEW2.columns
CCCASES=NEW2.Status=="AD"
CCCONTROLS=NEW2.Status=='CONTROL'
MMMCI=NEW2.Status=='MCI'


sum(CCCASES),sum(CCCONTROLS),sum(MMMCI)
#WIlcoxon with MCI
print("Only all three groups")
print('char')
print(stats.wilcoxon(NEW2[' char_path_len_x'], NEW2[' char_path_len_y'], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[' global_eff_x'], NEW2[' global_eff_y'], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[' transitivity _x'], NEW2[' transitivity _y'],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[' louvain_x'], NEW2[' louvain_y'], zero_method='wilcox', correction=False))
#only cases
print("Only cases")
print('char')
print(stats.wilcoxon(NEW2[CCCASES][' char_path_len_x'], NEW2[CCCASES][' char_path_len_y'], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[CCCASES][' global_eff_x'], NEW2[CCCASES][' global_eff_y'], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[CCCASES][' transitivity _x'], NEW2[CCCASES][' transitivity _y'],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[CCCASES][' louvain_x'], NEW2[CCCASES][' louvain_y'], zero_method='wilcox', correction=False))
#only controls
print("Only controls")
print('char')
print(stats.wilcoxon(NEW2[CCCONTROLS][' char_path_len_x'], NEW2[CCCONTROLS][' char_path_len_y'], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[CCCONTROLS][' global_eff_x'], NEW2[CCCONTROLS][' global_eff_y'], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[CCCONTROLS][' transitivity _x'], NEW2[CCCONTROLS][' transitivity _y'],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[CCCONTROLS][' louvain_x'], NEW2[CCCONTROLS][' louvain_y'], zero_method='wilcox', correction=False))
#only mci
print("Only mci")
print('char')
print(stats.wilcoxon(NEW2[MMMCI][' char_path_len_x'], NEW2[MMMCI][' char_path_len_y'], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[MMMCI][' global_eff_x'], NEW2[MMMCI][' global_eff_y'], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[MMMCI][' transitivity _x'], NEW2[MMMCI][' transitivity _y'],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[MMMCI][' louvain_x'], NEW2[MMMCI][' louvain_y'], zero_method='wilcox', correction=False))




#trying
#aanova non para baseline
print('char')
print(stats.kruskal(NEW2[CCCASES][' char_path_len_x'], NEW2[CCCONTROLS][' char_path_len_x'], NEW2[MMMCI][' char_path_len_x']))
print('global')
print(stats.kruskal(NEW2[CCCASES][' global_eff_x'], NEW2[CCCONTROLS][' global_eff_x'], NEW2[MMMCI][' global_eff_x']))
print('trans')
print(stats.kruskal(NEW2[CCCASES][' transitivity _x'], NEW2[CCCONTROLS][' transitivity _x'], NEW2[MMMCI][' transitivity _x']))
print('louvain')
print(stats.kruskal(NEW2[CCCASES][' louvain_x'], NEW2[CCCONTROLS][' louvain_x'],NEW2[MMMCI][' louvain_x']))
#anova non para followup
print('char')
print(stats.kruskal(NEW2[CCCASES][' char_path_len_y'], NEW2[CCCONTROLS][' char_path_len_y'], NEW2[MMMCI][' char_path_len_y']))
print('global')
print(stats.kruskal(NEW2[CCCASES][' global_eff_y'], NEW2[CCCONTROLS][' global_eff_y'], NEW2[MMMCI][' global_eff_y']))
print('trans')
print(stats.kruskal(NEW2[CCCASES][' transitivity _y'], NEW2[CCCONTROLS][' transitivity _y'], NEW2[MMMCI][' transitivity _y']))
print('louvain')
print(stats.kruskal(NEW2[CCCASES][' louvain_y'], NEW2[CCCONTROLS][' louvain_y'],NEW2[MMMCI][' louvain_y']))
#anova non para diff
print('char')
print(stats.kruskal(NEW2[CCCASES]['char_path_len'], NEW2[CCCONTROLS]['char_path_len'], NEW2[MMMCI]['char_path_len']))
print('global')
print(stats.kruskal(NEW2[CCCASES]['global_eff'], NEW2[CCCONTROLS]['global_eff'], NEW2[MMMCI]['global_eff']))
print('trans')
print(stats.kruskal(NEW2[CCCASES]['transitivity '], NEW2[CCCONTROLS]['transitivity '], NEW2[MMMCI]['transitivity ']))
print('louvain')
print(stats.kruskal(NEW2[CCCASES]['louvain'], NEW2[CCCONTROLS]['louvain'],NEW2[MMMCI]['louvain']))
#plot for the three boxplots
for y in [" louvain"," global_eff"," transitivity "," char_path_len"]:
    plt.clf()
    seaborn.set(style='ticks',font_scale=1.5)
    z=seaborn.boxplot(x="Visit", y=y, hue="Status",data=DTI2_infile,palette="Set3",order=["BASELINE","FOLLOWUP"], hue_order=["CONTROL","MCI","AD"])#,size=(4,4))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if y == " transitivity ":
        y=y[:-1]
    plt.savefig("Boxplot_withMCI"+y[1:]+"_exchanged.png") #avoid the space
#plot distrubution the three
##############################################################################################
#saving individually
 y = " louvain"
 y = " global_eff"
 y = " transitivity "
 y = " char_path_len"
    plt.clf()
    seaborn.set(style='ticks',font_scale=1.5)
    z=seaborn.boxplot(x="Visit", y=y, hue="Status",data=DTI2_infile,palette="Set3",order=["BASELINE","FOLLOWUP"], hue_order=["CONTROL","MCI","AD"])#,size=(4,4))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if y == " transitivity ":
        y=y[:-1]
##############################################################################################
#plot distribution again











#seaborn.set(rc={'figure.figsize':(3,3)})
seaborn.set(style='ticks',font_scale=1.5)
seaborn.pairplot(NEW2.loc[:,'transitivity ':'char_path_len'],kind="reg")
plt.savefig("Distribution_Differences_three.png")
plt.clf()
seaborn.set(style="ticks", color_codes=True,font_scale=1.5)
seaborn.set(rc={'figure.figsize':(1,1)})
#Plot features with Status
HUE=NEW2.loc[:,'Status':'char_path_len']
import seaborn as sns; sns.set(style="ticks", color_codes=True,font_scale=1.5)
g=seaborn.pairplot(HUE,hue="Status", markers=["o", "s","D"],diag_kind='kde',palette="muted")
plt.savefig("Distribution_Differences_Status_three.png")
plt.clf()

#Plot features with Visit
seaborn.pairplot(DTI2_infile,hue="Visit", markers=["o", "D"],diag_kind='kde',palette="muted")
plt.savefig("Distribution_Differences_Visit_three.png")
#Box plots
plt.clf()
seaborn.set_style('ticks')
seaborn.set(rc={'figure.figsize':(5,5)})
seaborn.boxplot(data=NEW.loc[:,'transitivity ':'char_path_len'],palette="RdBu")
plt.savefig("Boxplot_Differences.png")
#Plotting two categorical in a boxplot    
plt.clf()





#testing per group

#for i in ['CCCASES','CCCONTROLS','MMMCI']:
#CCCASES
print("The Wilcoxon test results among AD:")
print('char')
print(stats.wilcoxon(NEW2[' char_path_len_x'][CCCASES], NEW2[' char_path_len_y'][CCCASES], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[' global_eff_x'][CCCASES], NEW2[' global_eff_y'][CCCASES], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[' transitivity _x'][CCCASES], NEW2[' transitivity _y'][CCCASES],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[' louvain_x'][CCCASES], NEW2[' louvain_y'][CCCASES], zero_method='wilcox', correction=False))
#CCCONTROLS
print("The Wilcoxon test results among controls:")
print('char')
print(stats.wilcoxon(NEW2[' char_path_len_x'][CCCONTROLS], NEW2[' char_path_len_y'][CCCONTROLS], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[' global_eff_x'][CCCONTROLS], NEW2[' global_eff_y'][CCCONTROLS], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[' transitivity _x'][CCCONTROLS], NEW2[' transitivity _y'][CCCONTROLS],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[' louvain_x'][CCCONTROLS], NEW2[' louvain_y'][CCCONTROLS], zero_method='wilcox', correction=False))
#MMMCI
print("The Wilcoxon test results among MCI:")
print('char')
print(stats.wilcoxon(NEW2[' char_path_len_x'][MMMCI], NEW2[' char_path_len_y'][MMMCI], zero_method='wilcox', correction=False))
print('global')
print(stats.wilcoxon(NEW2[' global_eff_x'][MMMCI], NEW2[' global_eff_y'][MMMCI], zero_method='wilcox', correction=False))
print('trans')
print(stats.wilcoxon(NEW2[' transitivity _x'][MMMCI], NEW2[' transitivity _y'][MMMCI],zero_method='wilcox', correction=False))
print('louvain')
print(stats.wilcoxon(NEW2[' louvain_x'][MMMCI], NEW2[' louvain_y'][MMMCI], zero_method='wilcox', correction=False))

#checking:
i=CCCASES
NEW2[i]


#plotting for wilcoxon:
for y in [" louvain"," global_eff"," transitivity "," char_path_len"]:
    plt.clf()
    z=seaborn.boxplot(x="Status", y=y, hue="Visit",data=NEW2,palette="Set3",order=["CONTROL","MCI","AD"], hue_order=["BASELINE","FOLLOWUP"])
    if y == " transitivity ":
        y=y[:-1]
    plt.savefig("Boxplot_Differences_Status_all&Visit"+y[1:]+".png") #avoid the space
