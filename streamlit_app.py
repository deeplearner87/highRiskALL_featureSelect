#Run 'mapping_metadata_v1.ipynb' to create the metadata before using this notebook

#This notebook takes as input the following files:

#1. DRP data
#2. Clinical metadata
#3. Transcriptomics data
#4. Proteomics data


import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
import gdown

import warnings
warnings.filterwarnings('ignore')

st.write("""
# High risk ALL - DRP data analysis!
""")

if 'stage' not in st.session_state:
    st.session_state.stage = 0

def set_stage(stage):
    st.session_state.stage = stage
    

import os
 
dir = 'https://raw.githubusercontent.com/deeplearner87/highRiskALL_featureSelect/main/'
#os.chdir(dir)


# Drug response data
def create_groups(df, drugOfInterest):
    df = df.loc[df['drug']==drugOfInterest]
    # Convert column to numeric if not already
    df = df.copy()
    df['susceptibility_logAUC'] = pd.to_numeric(df['susceptibility_logAUC'], errors='coerce')
    
    # Define conditions for classification
    conditions = [
        (df['susceptibility_logAUC'] < 0.2),
        (df['susceptibility_logAUC'] >= 0.2) & (df['susceptibility_logAUC'] <= 0.75),
        (df['susceptibility_logAUC'] > 0.75)
    ]
    
    # Define values for each condition
    values = ['Sensitive', 'Intermediate', 'Resistant']
    
    # Use np.select to create the 'Class' column safely
    df = df.copy()  # Ensure we're working with a copy to avoid the warning
    df.loc[:, 'Class'] = np.select(conditions, values, default='Unknown')
    df.index = df['Labeling proteomics']
    df.drop(columns=['Labeling proteomics', 'drug'], inplace=True)
    return df

#Mapping Transcriptomics, Proteomics and DRP data against clinical metadata
def mapping_omicsandDRP2metadata(drugOfInterest):
    drp = pd.read_csv(dir+'Rank_Drugs_cleaned_only25_drugs_10122024.csv', header=0)
    drp['Labeling proteomics'] = drp['Labeling proteomics'].astype(str)
    drp.loc[:, 'Labeling proteomics'] = 'S' + drp['Labeling proteomics']
    #Removing rows corresponding to the contaminated sample '128'
    drp = drp.loc[drp['Labeling proteomics']!='S128']
    drp.drop(columns=['sample_id'], inplace=True)
    #Drop duplicates and keep only the first entry - not necessary any more as data has been cleaned already!
    #drp = drp.drop_duplicates(subset=['Labeling proteomics', 'drug'], keep='first')
    #Filtering for 25 drugs - not necessary any more as data has been filtered already!
    #drp = drp.loc[drp['drug'].isin(drugs_of_interest)]
    ##Check how many drugs each sample is treated with
    #drp.groupby('Labeling proteomics')['drug'].nunique()
    ##Check how many samples were treated with the 25-drugs panel
    #drp.groupby('drug')['Labeling proteomics'].nunique()
    drug_protein_df = create_groups(drp, drugOfInterest)
    drug_protein_df.index.names = ['Sample_ID']

    #Loading clinical metadata
    metadata = pd.read_csv(dir+'cleaned_metadata.csv')
    drug_df = drug_protein_df.reset_index()
    joined_df = metadata.merge(drug_df, how='inner', left_on='Sample ID Proteomics', right_on='Sample_ID')
    
    B_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe']== 'B-ALL', ['Sample ID Submitted', 'Sample ID Proteomics', 'Diagnosis/Relapse']]
    T_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe'] == 'T-ALL', ['Sample ID Submitted', 'Sample ID Proteomics', 'Diagnosis/Relapse']]
    
    #Loading the protein data
    #file_url = "https://drive.google.com/uc?id=1zOfgyP2ks6BnQRXttYfb7gn7qthooHPt"
    #Read the CSV file from Google Drive
    #protein = pd.read_csv(file_url, header=0, sep='\t', low_memory=False)
    
    #protein = pd.read_csv(dir+'Proteome_Atleast1validvalue_ImputedGD.txt', header=0, sep='\t', low_memory=False)
    #protein = protein.iloc[5:,:]
    protein_copy = protein.copy()
    protein.index = protein['Protein ID']
        
    protein = protein.iloc[:,0:127]
        
    T_ALL_protein_df = protein[protein.columns.intersection(T_ALL_samples['Sample ID Proteomics'])].T
    B_ALL_protein_df = protein[protein.columns.intersection(B_ALL_samples['Sample ID Proteomics'])].T

    #Loading Transcriptomics data
    #file_url = "https://drive.google.com/uc?id=1QmLl_ohlBm10Pd-kPy14POp39QbFrOLN"
    #Read the CSV file from Google Drive
    #rna = pd.read_csv(file_url, index_col=0)

    #rna = pd.read_csv(dir+'High-Risk-ALL_rna_preprocessed_protein_coding_genes.csv', index_col=0)

    B_ALL_rna_df = rna.loc[B_ALL_samples['Sample ID Submitted']]
    T_ALL_rna_df = rna.loc[T_ALL_samples['Sample ID Submitted']]
    
    drug_rna_df = joined_df
    drug_rna_df.index = drug_rna_df['Sample ID Submitted']
    drug_rna_df.index.names = ['Sample_ID']
    drug_rna_df.drop(columns=['Sample ID Submitted', 'Remarks (Dibyendu)', 'Sample ID Proteomics',
           'Immunophenoytpe', 'Diagnosis/Relapse', 'Sex', 'Age', 'ZNS',
           'Pred.resp.', 'Risiko MRD', 'Risk group', 'Calculated site of relapse',
           'Timepoint of relapse', 'matched pairs', 'Cytogenetics', 'Unnamed: 13',
           'Sample_ID'], inplace=True)
    return [drug_rna_df, drug_protein_df, T_ALL_rna_df, T_ALL_protein_df, B_ALL_rna_df, B_ALL_protein_df]

# Feature Selection
def preSelectFeatures(X, y, threshold, exp_name):
    import os
    X['Target'] = y
    corr_mat = pd.DataFrame(X.corr()['Target'])
    features = corr_mat.index[abs(corr_mat['Target']) >= threshold].tolist()   #consider both positive and negative correlations >=0.3 and <=-0.3
    return features[:-1]

def protein2gene(df, cols):
    #Loading the protein data
    #file_url = "https://drive.google.com/uc?id=1zOfgyP2ks6BnQRXttYfb7gn7qthooHPt"
    #Read the CSV file from Google Drive
    #protein = pd.read_csv(file_url, header=0, sep='\t', low_memory=False)
    
    #protein = pd.read_csv(dir+'Proteome_Atleast1validvalue_ImputedGD.txt', header=0, sep='\t', low_memory=False)
    #protein = protein.iloc[5:,:]
    protein_copy = protein.copy()
    protein.index = protein['Protein ID']
    protein = protein.iloc[:,0:127]
    protein2gene_mapping =  protein_copy[['Protein ID', 'Gene']]
    df = df[cols]
    genes = protein2gene_mapping.loc[protein2gene_mapping['Protein ID'].isin(df.columns), 'Gene']
    df.columns = genes
    df.columns = df.columns.astype(str)
    return df

def evaluateClassifiers(X, y):
    from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc
    from sklearn.model_selection import KFold, cross_validate, cross_val_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier
    from xgboost.sklearn import XGBClassifier
    from sklearn.svm import SVC
    from sklearn.linear_model import Lasso
    
    kfold = KFold(n_splits=10, shuffle=True, random_state=42) 
    models = [LogisticRegression(solver='liblinear', random_state=42), DecisionTreeClassifier(random_state=42), 
              RandomForestClassifier(random_state=42), SVC(kernel='linear', random_state=42), 
              XGBClassifier(random_state=42), Lasso(alpha=0.00001, random_state=42)]
    names = ["LR", "DT", "RF", "SVM", "XGB", "Lasso"]

    X = X.loc[:,~X.columns.duplicated()]
    for model, name in zip(models, names):
        print()
        print(name)
        for score in ["accuracy", "f1_weighted", "roc_auc", "r2", "neg_mean_absolute_error", "normalized_mutual_info_score", "neg_root_mean_squared_error", "explained_variance"]:
            results = cross_val_score(model, X.values, y, cv = kfold, scoring = score)
            #print(score,': {:.2f}'.format(results.mean()))


def importancePlot(feat_imp, exp_name):
    import matplotlib.pyplot as plt
    import os
    
    fig, ax = plt.subplots(figsize=(10,8))
    #feat_imp.plot.bar(yerr=std, ax=ax)
    feat_imp.plot.bar()
    ax.set_title("Feature importance_"+exp_name)
    ax.set_ylabel("Score")
    ax.set_xlabel('Gene')
    filename = exp_name+'_feature_importance_based_on_DRP.pdf'
    #plt.savefig(os.path.join(dir, 'Results/')+filename, dpi = 300, format = 'pdf', bbox_inches="tight")

def differentialPlot(df, conditions, exp_name):
    import scanpy as sc
    import anndata
    import os
    from scipy import stats
    import matplotlib.pyplot as plt
    cols = df.columns.astype(str)
    samp = df.index
    X = pd.DataFrame(np.array(df, dtype=float))
    X.columns = cols
    X.index = samp
    X = stats.zscore(X)
    ad = anndata.AnnData(X)
    ad.obs = pd.DataFrame(conditions, columns=['class'])
    ad.var_names = X.columns
    ad.var_names_make_unique()
    filename = exp_name+'_heatmap_based_on_DRP.pdf'
    with plt.rc_context():
        ax = sc.pl.heatmap(ad, ad.var_names, groupby='class', swap_axes=True, show_gene_labels=True, cmap="PiYG_r", show=False)
        ax['heatmap_ax'].set_ylabel("Gene")
        st.pyplot(plt.gcf())
        #plt.savefig(os.path.join(dir,'Results/')+filename, dpi = 300, format = 'pdf', bbox_inches="tight")

def majorityVoting(lof): #Input the list of features
    from collections import Counter
    cnt=Counter()
    for x in lof:
        cnt+= Counter(x)
    features = [k for k,v in dict(cnt).items() if v>=2] #selecting a feature if 2 or more classifiers agree
    return features

def selectFeatures(df, classifiers, exp_name, n):
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.tree import DecisionTreeClassifier
    from xgboost.sklearn import XGBClassifier
    from sklearn.svm import SVC
    from sklearn.linear_model import Lasso
    
    X = df[0]
    y = df[1]

    X = X.loc[:,~X.columns.duplicated()]
    
    #Train-test split
    #X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)
    
    feat_imp=[]
    acc = []
    if 'LR' in classifiers:
        #Train a Logistic Regression (LR) model
        lr = LogisticRegression(solver='liblinear', random_state=42)
        lr.fit(X, y)
        feat_imp_lr = pd.Series(lr.coef_[0], index=X.columns, name='LR')
        feat_imp_lr = feat_imp_lr.sort_values(ascending=False)
        feat_imp.append(feat_imp_lr.index)
        acc.append(lr.score(X.values, y))
        #importancePlot(feat_imp_lr, exp_name+'_LR') #Plot feature importances
    
    if 'DT' in classifiers:
        #Train a Decision Tree (DT) model
        dt = DecisionTreeClassifier(random_state=42)
        dt.fit(X, y)
        feat_imp_dt = pd.Series(dt.feature_importances_, index=X.columns, name='DT')
        feat_imp_dt = feat_imp_dt.sort_values(ascending=False)
        feat_imp.append(feat_imp_dt.index)
        acc.append(dt.score(X.values, y))
        #importancePlot(feat_imp_dt, exp_name+'_DT') #Plot feature importances
    
    if 'RF' in classifiers:
        #Train a Random Forest (RF) model  
        rf = RandomForestClassifier(random_state=42)
        rf.fit(X, y)
        feat_imp_rf = pd.Series(rf.feature_importances_, index=X.columns, name='RF')
        feat_imp_rf = feat_imp_rf.sort_values(ascending=False)
        feat_imp.append(feat_imp_rf.index)
        acc.append(rf.score(X.values, y))
        #importancePlot(feat_imp_rf, exp_name+'_RF') #Plot feature importances
        
    if 'SVC' in classifiers:
        #Train a SVM classifier model
        svc = SVC(kernel='linear', random_state=42)
        svc.fit(X, y)
        feat_imp_svc = pd.Series(svc.coef_[0], index=X.columns, name='SVC')
        feat_imp_svc = feat_imp_svc.sort_values(ascending=False)
        feat_imp.append(feat_imp_svc.index)
        acc.append(svc.score(X.values, y))
        #importancePlot(feat_imp_svc, exp_name+'_SVC') #Plot feature importances

    if 'XGB' in classifiers:
        #Train a XGBoost classifier model
        xgb = XGBClassifier(random_state=42)
        xgb.fit(X, y)
        feat_imp_xgb = pd.Series(xgb.feature_importances_, index=X.columns, name='XGB')
        feat_imp_xgb = feat_imp_xgb.sort_values(ascending=False)
        feat_imp.append(feat_imp_xgb.index)
        acc.append(xgb.score(X.values, y))
        #importancePlot(feat_imp_xgb, exp_name+'_XGB') #Plot feature importances

    if 'Lasso' in classifiers:
        #Train a Lasso model
        lasso = Lasso(alpha=0.00001, random_state=42)
        lasso.fit(X, y)
        feat_imp_lasso = pd.Series(lasso.coef_, index=X.columns, name='Lasso')
        feat_imp_lasso = feat_imp_lasso.sort_values(ascending=False)
        feat_imp.append(feat_imp_lasso.index)
        acc.append(lasso.score(X.values, y))
        #importancePlot(feat_imp_lasso, exp_name+'_Lasso') #Plot feature importances

    #print(acc)
    #print(f"Average training accuracy: {np.mean(acc):.2f}")
    selFeatures = majorityVoting(feat_imp)
    return selFeatures[0:n]

def classify(data, drug_data, exp_name, classifiers, num_features, threshold, omics_type):
    from sklearn.preprocessing import LabelEncoder

    drug_data = drug_data.loc[data.index] #select samples based on data i.e., (T-ALL or B-ALL)
    drug_data = drug_data.loc[drug_data['Class']!='Intermediate'] #filter out all sample with 'Intermediate' class
    data = data.loc[drug_data.index] #filter out samples corresponding to the 'Intermediate' class from the data
    labels = drug_data['Class']

    le = LabelEncoder()
    y = le.fit_transform(np.ravel(labels))
    cols = data.columns.astype(str)
    samples = data.index
    X = pd.DataFrame(np.array(data, dtype=float))
    X.columns = cols
    X.index = samples
    #print('Shape: ', X.shape)
    #preselect features based on correlation with target variable from entire protein expression data
    feat = preSelectFeatures(X, y, threshold, exp_name)
    if feat:
        X = X[feat]
    if omics_type == 'Proteomics':
        print('{} proteins were found to have significant positive or negative correlation with the annotations.'.format(len(feat)))
        X = protein2gene(X, X.columns)
    elif omics_type == 'Transcriptomics':
        print('{} genes were found to have significant positive or negative correlation with the annotations.'.format(len(feat)))
        #X = ensemblID2Gene(X, X.columns)
    #ans=input('Would you like to continue? (Y/N)')
    #if ans=='n' or ans=='N':
    #    return
    #evaluateClassifiers(X,y)
    selFeatures = selectFeatures([X,y], classifiers, exp_name, num_features)
    X = X[selFeatures]
    differentialPlot(X, labels.values, exp_name)
    X['Drug_Class']=labels
    return selFeatures

omics_type = st.selectbox('Select omics-type: ', ['Proteomics', 'Transcriptomics'])
if omics_type=='Transcriptomics':
    uploaded_file = st.file_uploader("Upload transcriptomics data")
    if uploaded_file:
        rna = pd.read_csv(uploaded_file, index_col=0)
        st.write("Uploaded transcriptomics data:")
        st.dataframe(rna)
elif omics_type=='Proteomics':
    uploaded_file = st.file_uploader("Upload proteomics data")
    if uploaded_file:
        protein = pd.read_csv(uploaded_file, header=0, sep='\t', low_memory=False)
        protein = protein.iloc[5:,:]
        st.write("Uploaded proteomics data:")
        st.dataframe(protein)
    
cell_type = st.selectbox('Select cell-type: ', ['T-ALL', 'B-ALL'])
drugs_of_interest = ['Idarubicin', 'Dasatinib', 'Ponatinib', 'Venetoclax', 'Navitoclax', 'Doxorubicin', 'Birinapant', 'Bortezomib', 'CB-103', 'Dexamethasone', 'Cytarabine', 'Etoposide', 'Methotrexate', 'Selinexor', 'Vincristine', 'Nilotinib', 'Temsirolimus', 'Bosutinib', 'Panobinostat', 'Trametinib', 'Ruxolitinib', 'Dinaciclib', 'A1331852', 'S-63845', 'Nelarabine']
drugOfInterest = st.selectbox('Select drug', options=[opt.strip() for opt in drugs_of_interest])
num_features = st.slider('Select number of genes you want to select',1, 100, 50)
threshold = st.slider('Select threshold for correlation-based feature pre-selection', 0.00, 1.00, 0.55) #threshold for correlation-based preselection
classifiers = st.multiselect('Select models - You may choose multiple among the following: [Logistic Regression, Decision Tree Classifier, Random Forest Classifier, Support Vector Machine Classifer, XG Boost Classifier and Lasso Regression]', ['LR', 'DT', 'RF', 'SVC', 'XGB', 'Lasso'])
#st.write(classifiers)

if uploaded_file is not None:
    data = mapping_omicsandDRP2metadata(drugOfInterest)
    
    if omics_type == 'Transcriptomics':
        drug_data = data[0]
    elif omics_type == 'Proteomics':
        drug_data = data[1]
    
    if cell_type == 'T-ALL':
        if omics_type == 'Transcriptomics':
            data = data[2]
        elif omics_type == 'Proteomics':
            data = data[3]
    elif cell_type == 'B-ALL':
        if omics_type == 'Transcriptomics':
            data = data[4]
        elif omics_type == 'Proteomics':
            data = data[5]
    
    analyze = st.button('Analyze', on_click=set_stage, args=(1,))
    if analyze:
        if len(classifiers) < 2:
            st.write('Please select at least 2 classifiers')
        else:
            #st.write(st.session_state)
            exp_name = cell_type+'_'+omics_type+'_'+drugOfInterest+'_'
            selFeatures = classify(data, drug_data, exp_name, classifiers, num_features, threshold, omics_type)
else:
    st.warning("Please upload a file to proceed!")
