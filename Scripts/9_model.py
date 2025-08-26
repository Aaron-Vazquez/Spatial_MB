# Load Libraries
print('Processing...')
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import pickle5 as pickle
import xgboost as xgb
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib import pyplot
import os
import shap

#######################################################################################
# Preparation of Matrix 
# Load data (X and Y)

os.chdir("./")

### raw data
X = pd.read_csv('./Results/ML/Data_counts.csv', 
                header = 0, index_col = 0);
#X = np.transpose(X)

# 0 -> G4
# 1 -> SHH
# 2 -> G3
Y = pd.read_csv('./Results/ML/classes.csv',
                header = 0);

#Y = pd.read_csv('./Results/ML/labels.tsv',
#                sep='\t', header=None)
### sctranformed data
#X_sct = pd.read_csv('./Data/MCTS_sct.csv', 
#                header = 0, index_col = 0);

### scaled data (lop1p(RPKM))
#X_sld = pd.read_csv('./Data/MCTS_scaled.csv', 
#                header = 0, index_col = 0);

random = [0,1,2,3,12]
###############################################################################################################
for j in range(5):

#Split the data and Save X, y train and test.
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.25, 
                                                        random_state = random[j],
                                                        stratify = Y);
    
    
    ## Data separation
    #Original data
    # Y.plot.hist(grid=True, bins=3, rwidth=0.9)
    # Y_train.plot.hist(grid=True, bins=3, rwidth=0.9)
    # Y_test.plot.hist(grid=True, bins=3, rwidth=0.9)
    
    
    os.system('mkdir .  /Results/ML/Splited_data')
    pickle.dump(X_train, open("Results/ML/Splited_data/X_train_"+str(j)+".pickle", 'wb'), protocol=4);
    pickle.dump(X_test, open("Results/ML/Splited_data/X_test_"+str(j)+".pickle", 'wb'), protocol=4);
    
    ################################################################################################################
    # Constructing the XGBoost model. 
    print('Constructing model..')
    
    D_train = xgb.DMatrix(X_train, label=Y_train);
    D_test = xgb.DMatrix(X_test, label=Y_test);
    
    
    #parameters
    params = {
    #    'max_depth': 4,
        'eta': 0.2,
        'eval_metric':'auc',
        'objective': 'multi:softprob',  
        #### binary:hinge just for two classes
        'seed':90,
        #### number of classes
        'num_class': 3,
        'n_gpus': 0
    }
    
    steps = 10000
    
    model = xgb.train(params, D_train, steps)
    #Prepare the format
    # model_bytearray = model.save()[4:]
    # def myfun(self=None):
    #     return model_bytearray
    
    # model.save_raw = myfun
    
    ##################################################################################################################
    # Save model
    os.system('mkdir Results/ML/Model')
    pickle_out = open("Results/ML/Model/model_"+str(j)+".pickle","wb")
    pickle.dump(model, pickle_out)
    pickle_out.close()
    
    ###################################################################################################################
    # Assessment of the model. Confusion matrix
    # Performance of the model evaluating the test dataset.
    pred= model.predict(D_test)
    pred_data=pd.DataFrame(pred)
    pred_data.to_csv('Results/ML/figures/prediction_'+str(j)+'.csv', index=True)
    predictions = [np.round(value) for value in pred]
    
    def classification(value):
        if sum(value)==0 or sum(value)>1:
            val = 3
        else:
            val = np.where(value==1)[0][0]
        return val
    
    Pred_vals = [classification(value) for value in predictions]
    print(classification_report(Y_test, Pred_vals, 
                                zero_division = "warn"))
    cm = confusion_matrix(Y_test, Pred_vals)
    cm
    
    cmap = plt.get_cmap('Blues')
    
    os.system('mkdir Results/ML/figures')
    
    def plot_confusion_matrix(cm, classes, normalized=True, cmap=cmap): # 'bone'
        plt.figure(figsize=[7, 6])
        norm_cm = cm
        if normalized:
            norm_cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            where_are_NaNs = np.isnan(norm_cm)
            norm_cm[where_are_NaNs] = 0
    #        print(norm_cm)
            sns.heatmap(norm_cm, annot=cm, fmt='g', xticklabels=classes, yticklabels=classes, cmap=cmap)
            pyplot.savefig('Results/ML/figures/confusematrix_'+str(j)+'.png',dpi=800)
    
    # os.system('mkdir figures')
    
    # # I will save the figure data to make the figure in the paper.
    # # training and test performance
    # pickle_out = open("figures/cm.pickle","wb")
    # pickle.dump(cm, pickle_out)
    # pickle_out.close()
    
    plot_confusion_matrix(cm, ['Grupo4', 'SHH','Grupo3','Undefined'])
    # pyplot.savefig('figures/confusematrix_raw.png',dpi=800) 
    print('Confuse matrix saved!')
    
    
    ####################################################################################################################
    # Figures of the paper.
    
    # import matplotlib.pylab as plt
    
    # SHAP analysis
    print('Shap plots in progress..')
 #   explainer = shap.TreeExplainer(model)
    shap_values = shap.TreeExplainer(model).shap_values(X_train)
    ### for multiclass only bar available
    shap.summary_plot(shap_values, X_train, plot_type="bar",max_display = 30)
    pyplot.savefig('Results/ML/figures/shap_'+str(j)+'.png',dpi=800)
    
    ############# 
    ### ordering shap values
    ## a = np.mean(np.abs(shap_values[1]), axis = 0)
    ## b = np.argsort(a)
    ## X.columns[np.flip(b)]
    os.system('mkdir Results/ML/SHAP_values')
    letters = ['Grupo4', 'SHH','Grupo3']
    for i in range(3):
        a = np.mean(np.abs(shap_values[i]), axis = 0)
        b = np.argsort(a)
        d = {'Gene':X.columns[np.flip(b)], 
              'Shap':a[np.flip(b)]}
        df = pd.DataFrame(data=d)
        df.to_csv("SHAP_values/SHAP_"+letters[i]+"_"+str(j)+".csv")
    #    print(X.columns[np.flip(b)][0:10])
        print(df)
    ############
# # Save shap_values
# os.system('mkdir pickle_shap')
# pickle.dump(shap_values, open("pickle_shap/shap_values.pickle", 'wb'), protocol=4);

# X_display = X

# plt.clf()
# shap.summary_plot(shap_values, X)
# pyplot.savefig('figures/shap_summary1.png',format='png', dpi=800)

# plt.clf()
# shap.summary_plot(shap_values, X, plot_type="bar")
# pyplot.savefig('figures/shap_summary1_bar.png',format='png', dpi=800, bbox_inches='tight')
# print('Shap plots saved!!')
