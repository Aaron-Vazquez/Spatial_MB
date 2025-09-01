# Load Libraries
print('Processing...')
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import pickle as pickle
import xgboost as xgb
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib import pyplot
import os
import shap



### Consensus genes 5 Samples
d = {'SHH_Genes': ["EYA1","ISLR2","GABRD","PARM1","PPP1R1A","NTRK2","DPYSL3",
"GALNT9","GRIK1","TRH","LHX4","RGMB","EPHA8","NKX6-1","GFRA2","ENC1","RELN","KIF21B","SPATS2L"],
'Grp3_Genes':["OTX2","PABPC1","EEF2","UBC","EYS","GNB1","GPM6A","UBA52","TAFA4",
"CRX","SFRP1","EIF4A2","STMN2","CCNI","EIF4B"],
'Grp4_Genes':["ISLR2","MAP1B","RGMB","ACTB","PPP1R1A","BOK","APOE","TMSB4X",
"UCHL1","NNAT","MARCKSL1","GABRA6","ZIC1","HIST1H4D","UBE2QL1","TRH","LMO4","PARM1"]}

os.chdir("./")
for j in range(5):
    model = pickle.load(open('./Results/ML/Model/model_'+str(j)+'.pickle','rb'))
    X_train = pickle.load(open('./Results/ML/Splited_data/X_train_'+str(j)+'.pickle','rb'))
    print('Shap plots in progress..')
    explainer = shap.TreeExplainer(model)
#    shap_values = pd.read_csv("./Results/ML/SHAP_values/SHAP_Grupo4_cln_0.csv")
    shap_values = shap.TreeExplainer(model).shap_values(X_train, check_additivity=False)
#### Grp4 Genes
    fig2, axs = plt.subplots(nrows=int(np.ceil(len(d["Grp4_Genes"])/4)), ncols=4,
                             figsize=(10, 6), constrained_layout=True)
    axs = axs.ravel()
    for i, gene in enumerate(d['Grp4_Genes']):
        shap.dependence_plot(gene, shap_values, X_train, ax=axs[i], 
                             show = False, interaction_index = None)
    
    [fig2.axes[-1].remove() for i in range(np.ceil(len(d["Grp4_Genes"])/4)*4 - len(d["Grp4_Genes"]))]

    fig2.savefig('./Results/ML/figures/Shap_dependence_plot_Grp4_'+str(j)+'.png')

#### SHH Genes
    fig2, axs = plt.subplots(nrows=np.ceil(len(d["SHH_Genes"])/4), ncols=4,
                                 figsize=(10, 6), constrained_layout=True)
    axs = axs.ravel()
    for i, gene in enumerate(d['SHH_Genes']):
        shap.dependence_plot(gene, shap_values[1], X_train, ax=axs[i], 
        show = False, interaction_index = None)
    
    [fig2.axes[-1].remove() for i in range(np.ceil(len(d["SHH_Genes"])/4)*4 - len(d["SHH_Genes"]))]

    fig2.savefig('./Results/ML/figures/Shap_dependence_plot_SHH_'+str(j)+'.png')

#### Grp3 Genes
    fig2, axs = plt.subplots(nrows=np.ceil(len(d["Grp3_Genes"])/4), ncols=4,
                                 figsize=(10, 6), constrained_layout=True)
    axs = axs.ravel()
    for i, gene in enumerate(d['Grp3_Genes']):
        shap.dependence_plot(gene, shap_values[2], X_train, ax=axs[i], 
        show = False, interaction_index = None)
    
    [fig2.axes[-1].remove() for i in range(np.ceil(len(d["Grp3_Genes"])/4)*4 - len(d["Grp3_Genes"]))]

    fig2.savefig('./Results/ML/figures/Shap_dependence_plot_Grp3_'+str(j)+'.png')




   
