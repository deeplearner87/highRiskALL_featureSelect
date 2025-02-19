#Run 'mapping_metadata_v2.ipynb' to create the metadata before using this notebook

#This notebook takes as input the following files:

#1. DRP data
#2. Clinical metadata
#3. Transcriptomics data
#4. Proteomics data
#5. ProteinID2Gene mapping


import os
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
#import gdown
import requests

import warnings
warnings.filterwarnings('ignore')

st.write("""
# High risk ALL - DRP data analysis!
""")

if 'stage' not in st.session_state:
    st.session_state.stage = 0
    
def set_stage(stage):
    st.session_state.stage = stage
    
# Drug response data
def create_groups(df, drugOfInterest):
    df = df.loc[df['drug']==drugOfInterest]
    #Convert column to numeric if not already
    df = df.copy()
    df['susceptibility_logAUC'] = pd.to_numeric(df['susceptibility_logAUC'], errors='coerce')
    
    #Define conditions for classification
    conditions = [
        (df['susceptibility_logAUC'] < 0.2),
        (df['susceptibility_logAUC'] >= 0.2) & (df['susceptibility_logAUC'] <= 0.75),
        (df['susceptibility_logAUC'] > 0.75)
    ]
    
    #Define values for each condition
    values = ['Sensitive', 'Intermediate', 'Resistant']
    
    #Use np.select to create the 'Class' column safely
    df = df.copy()  # Ensure we're working with a copy to avoid the warning
    df.loc[:, 'Class'] = np.select(conditions, values, default='Unknown')
    df.index = df['Labeling proteomics']
    df.drop(columns=['Labeling proteomics', 'drug'], inplace=True)
    return df

#Mapping Proteomics and DRP data against clinical metadata
def mapping_Proteomics_DRP_to_metadata(drugOfInterest):
    drp_data_url = st.secrets["data_links"]["drp_data"]
    #Download the file
    response = requests.get(drp_data_url)
    if response.status_code == 200:
        with open("drp.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    drp = pd.read_csv("drp.csv", header=0)
    drp['Labeling proteomics'] = drp['Labeling proteomics'].astype(str)
    drp.loc[:, 'Labeling proteomics'] = 'S' + drp['Labeling proteomics']
    #Removing rows corresponding to the contaminated sample '126'
    drp = drp.loc[drp['Labeling proteomics']!='S126']
    drp.drop(columns=['sample_id'], inplace=True)
    drug_protein_df = create_groups(drp, drugOfInterest)
    drug_protein_df.index.names = ['Sample_ID']
    """
    #Loading clinical metadata
    clinical_metadata_url1 = st.secrets["data_links"]["clinical_metadata_KR"]
    #Download the file
    response = requests.get(clinical_metadata_url)
    
    if response.status_code == 200:
        with open("metadata.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    metadata = pd.read_csv("metadata.csv", header=0)
    """
    
    clinical_metadata_url1 = st.secrets["data_links"]["clinical_metadata_KR"]
    response = requests.get(clinical_metadata_url1)

    if response.status_code == 200:
        with open("metadata.csv", "wb") as f:
            f.write(response.content)
            st.write("Download successful.")
    else:
            st.write(f"Failed to download file. Status code: {response.status_code}")
    
    
    """
    try:
        metadata = pd.read_excel("metadata.xlsx", engine="openpyxl")
        st.write(metadata.head())  # Print first few rows to confirm it's loaded
    except Exception as e:
        st.write(f"Error reading Excel file: {e}")
    """
    metadata['Sample ID Proteomics'] = metadata['Sample ID Proteomics'].astype('str')
    metadata['Sample ID Proteomics'] = 'S'+metadata['Sample ID Proteomics']
    metadata['Immunophenoytpe']=metadata['Immunophenoytpe'].astype('str')
    metadata.drop(metadata[metadata['Sample ID Proteomics'] == 'S126'].index, inplace=True) #dropping the contaminated sample
    metadata.loc[metadata['Immunophenoytpe'].isin(['T-ALL', 'T-ALL ', 'T-LBL/T-ALL']), 'Immunophenoytpe'] = 'T-ALL'
    metadata.loc[metadata['Sample ID Proteomics']=='S108', 'Diagnosis/Relapse'] = 'Diagnosis' #information collected from protein expression data
    
    drug_df = drug_protein_df.reset_index()
    joined_df = metadata.merge(drug_df, how='inner', left_on='Sample ID Proteomics', right_on='Sample_ID')
    
    B_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe']== 'B-ALL', ['Sample ID Proteomics', 'Diagnosis/Relapse']]
    T_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe'] == 'T-ALL', ['Sample ID Proteomics', 'Diagnosis/Relapse']]
    
    #Loading the protein data
    #protein_vsn_url = st.secrets["data_links"]["protein_vsn"]
    #protein_no_vsn_url = st.secrets["data_links"]["protein_no_vsn"]
    protein = st.secrets["data_links"]["protein"]

    #Download the file
    response = requests.get(protein)
    if response.status_code == 200:
        with open("protein.csv", "wb") as f:
            f.write(response.content)
    
    #Read the file
    protein = pd.read_csv("protein.csv", header=0, sep='\t')
    
    protein = protein.iloc[5:,:]
    protein_copy = protein.copy()
    protein.index = protein['Protein ID']
    #protein.head()
    protein2gene_mapping = protein_copy[['Protein ID', 'Gene']]
    T_ALL_protein_df = protein[protein.columns.intersection(T_ALL_samples['Sample ID Proteomics'])].T
    B_ALL_protein_df = protein[protein.columns.intersection(B_ALL_samples['Sample ID Proteomics'])].T
    
    return [drug_protein_df, T_ALL_protein_df, B_ALL_protein_df]

    
def mapping_Transcriptomics_DRP_to_metadata(drugOfInterest):
    drp_data_url = st.secrets["data_links"]["drp_data"]
    #Download the file
    response = requests.get(drp_data_url)
    if response.status_code == 200:
        with open("drp.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    drp = pd.read_csv("drp.csv")
    drp['Labeling proteomics'] = drp['Labeling proteomics'].astype(str)
    drp.loc[:, 'Labeling proteomics'] = 'S' + drp['Labeling proteomics']
    #Removing rows corresponding to the contaminated sample '126'
    drp = drp.loc[drp['Labeling proteomics']!='S126']
    drp.drop(columns=['sample_id'], inplace=True)
    drug_protein_df = create_groups(drp, drugOfInterest)
    drug_protein_df.index.names = ['Sample_ID']
    
    #Loading clinical metadata
    clinical_metadata_url = st.secrets["data_links"]["clinical_metadata"]
    #Download the file
    response = requests.get(clinical_metadata_url)
    if response.status_code == 200:
        with open("metadata.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    metadata = pd.read_csv("metadata.csv", header=0)
    
    drug_df = drug_protein_df.reset_index()
    joined_df = metadata.merge(drug_df, how='inner', left_on='Protein_Sample_ID', right_on='Sample_ID')
    
    B_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe']== 'B-ALL', ['RNA_Sample_ID_Available', 'Protein_Sample_ID', 'Diagnosis/Relapse']]
    T_ALL_samples = joined_df.loc[joined_df['Immunophenoytpe'] == 'T-ALL', ['RNA_Sample_ID_Available', 'Protein_Sample_ID', 'Diagnosis/Relapse']]
    
    #Loading the protein data
    #protein_vsn_url = st.secrets["data_links"]["protein_vsn"]
    #protein_no_vsn_url = st.secrets["data_links"]["protein_no_vsn"]
    protein = st.secrets["data_links"]["protein"]

    #Download the file
    response = requests.get(protein)
    if response.status_code == 200:
        with open("protein.csv", "wb") as f:
            f.write(response.content)
    
    #Read the file
    protein = pd.read_csv("protein.csv", header=0, sep='\t')
    
    protein = protein.iloc[5:,:]
    protein_copy = protein.copy()
    protein.index = protein['Protein ID']
    #protein.head()
    T_ALL_protein_df = protein[protein.columns.intersection(T_ALL_samples['Protein_Sample_ID'])].T
    B_ALL_protein_df = protein[protein.columns.intersection(B_ALL_samples['Protein_Sample_ID'])].T
    
    #Loading Transcriptomics data
    rna_url = st.secrets["data_links"]["rna"]
    #Download the file
    response = requests.get(rna_url)
    if response.status_code == 200:
        with open("rna.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    rna = pd.read_csv("rna.csv", index_col=0)
    B_ALL_rna_df = rna.loc[B_ALL_samples['RNA_Sample_ID_Available']]
    T_ALL_rna_df = rna.loc[T_ALL_samples['RNA_Sample_ID_Available']]
    
    drug_rna_df = joined_df
    drug_rna_df.index = drug_rna_df['RNA_Sample_ID_Available']
    drug_rna_df.index.names = ['Sample_ID']
    drug_rna_df.drop(columns=['Sample_ID_Submitted', 'RNA_Sample_ID_Available',
       'Duplicate_RNA_Sample_Available?', 'Protein_Sample_ID', 'WES_Sample_ID',
       'Immunophenoytpe', 'Diagnosis/Relapse', 'Sex', 'Age', 'ZNS',
       'Pred.resp.', 'Risiko MRD', 'Risk group', 'Calculated site of relapse',
       'Timepoint of relapse', 'matched pairs', 'Cytogenetics', 'Remarks',
       'Sample_ID'], inplace=True)
    return [drug_rna_df, drug_protein_df, T_ALL_rna_df, T_ALL_protein_df, B_ALL_rna_df, B_ALL_protein_df]



def selectFeaturesByWilcoxon(X, y, alpha, exp_name):
    from scipy.stats import ranksums
    from statsmodels.stats.multitest import multipletests
    # Ensure y is a pandas Series
    y = pd.Series(y, name="Target", index=X.index)

    # Check and align X and y
    if not X.index.equals(y.index):
        raise ValueError("Indices of X and y do not match. Align them before running this function.")

    # Prepare result storage
    results = []

    # Perform Wilcoxon rank-sum test for each feature
    for feature in X.columns:
        group_0 = X.loc[y == 0, feature]  # Resistant group
        group_1 = X.loc[y == 1, feature]  # Sensitive group
        
        # Check if there are enough data points in both groups
        if len(group_0) > 1 and len(group_1) > 1:
            stat, p_value = ranksums(group_0, group_1)
            results.append((feature, stat, p_value))
        else:
            results.append((feature, np.nan, np.nan))  # Append NaN if test cannot be run

    # Convert results to DataFrame
    results_df = pd.DataFrame(results, columns=["Feature", "Statistic", "P_Value"])

    # Adjust p-values for multiple testing
    results_df["Adjusted_P_Value"] = multipletests(results_df["P_Value"], method="fdr_bh")[1]

    # Select features with significant adjusted p-values
    significant_features = results_df[results_df["P_Value"] <= alpha]['Feature'].tolist()

    # Save results to CSV
    output_path = os.path.join("Results", f"{exp_name}_wilcoxon_test_results.csv")
    #os.makedirs("Results", exist_ok=True)  # Ensure directory exists
    results_df.to_csv(output_path, index=False)

    return significant_features


def preSelectFeatures(X, y, threshold, exp_name):
    # Ensure y is a pandas Series
    y = pd.Series(y, name="Target", index=X.index)

    # Check and handle issues in X and y
    if not X.index.equals(pd.Index(y.index)):
        raise ValueError("Indices of X and y do not match. Align them before running this function.")

    # Drop non-numeric or constant columns
    X = X.apply(pd.to_numeric, errors="coerce").loc[:, X.nunique() > 1]

    # Drop rows with NaN in X or y
    X = X.dropna()
    y = y.loc[X.index]

    # Compute correlation of each feature with the target
    corr_mat = X.corrwith(y)

    # Create a DataFrame to save correlations
    corr_df = corr_mat.reset_index()
    corr_df.columns = ["Feature", "Correlation"]
    #print(corr_df)
    
    # Ensure 'Correlation' column is numeric and drop NaNs
    corr_df["Correlation"] = pd.to_numeric(corr_df["Correlation"], errors="coerce")
    corr_df = corr_df.dropna(subset=["Correlation"])

    # Save correlation matrix to CSV
    output_path = os.path.join("Results", f"{exp_name}_correlation_with_target_DRP.csv")
    corr_df.to_csv(output_path, index=False)

    # Select features with absolute correlation above the threshold
    selected_features = corr_df.loc[abs(corr_df["Correlation"]) >= threshold, "Feature"].tolist()
    return selected_features


def protein2gene(df, cols):
    protein2gene_url = st.secrets["data_links"]["proteinIDtoGene"]
    #Download the file
    response = requests.get(protein2gene_url)
    if response.status_code == 200:
        with open("protein2gene.csv", "wb") as f:
            f.write(response.content)
    #Read the file
    protein2gene_mapping = pd.read_csv("protein2gene.csv")
    protein2gene_mapping['Protein.ID'] = protein2gene_mapping['Protein.ID'].str.strip()  #Strip whitespaces
    cols = [col.strip() for col in cols] #Strip whitespaces
    genes = protein2gene_mapping.loc[protein2gene_mapping['Protein.ID'].isin(cols), 'Gene']
    df = df[cols]
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
            print(score,': {:.2f}'.format(results.mean()))

def importancePlot(feat_imp, exp_name):
    import matplotlib.pyplot as plt
    
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
    ad.var_names = cols
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
    features = [k for k,v in dict(cnt).items() if v>2] #selecting a feature if 3 or more classifiers agree
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
"""
def classify(data, drug_data, exp_name, classifiers, num_features, threshold, omics_type):
    from sklearn.preprocessing import LabelEncoder

    drug_data = drug_data.loc[data.index] #select samples based on data i.e., (T-ALL or B-ALL)
    drug_data = drug_data.loc[drug_data['Class']!='Intermediate'] #filter out all sample with 'Intermediate' class
    data = data.loc[drug_data.index] #filter out samples corresponding to the 'Intermediate' class from the data
    labels = drug_data['Class']

    #label = pd.DataFrame(drug_data['Class'])
    #label = label[label['Class']!='Intermediate']
    #drug_data = drug_data.loc[drug_data['Class']!='Intermediate']
    #if omics_type=='Transcriptomics':
    #    samples = label
    #samples = data.index.intersection(label.index) #extracting sample IDs for drug classes 'Sensitive' and 'Resistant'
    #print(label)
    #print(samples)
    #X = data.loc[samples]
    #label = label.loc[samples]
    
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
    #X.to_csv(os.path.join(dir,'Results/')+exp_name+'DRP_ML_selFeatures_with_annotations.csv')
    return selFeatures
"""
def classify_proteomics(data, drug_data, exp_name, drugOfInterest, classifiers, num_features, threshold):
    from sklearn.preprocessing import LabelEncoder
    drug_data = drug_data.loc[data.index] #select samples based on data i.e., (T-ALL or B-ALL)
    #print(drug_data.shape)
    drug_data = drug_data.loc[drug_data['Class']!='Intermediate'] #filter out all sample with 'Intermediate' class
    data = data.loc[drug_data.index] #filter out samples corresponding to the 'Intermediate' class from the data
    labels = drug_data['Class']
    
    le = LabelEncoder()
    y = le.fit_transform(labels)
    cols = data.columns.astype(str)
    samples = data.index
    X = pd.DataFrame(np.array(data, dtype=float))
    
    X.columns = cols
    X.index = samples

    #preselect features based on correlation with target variable from entire protein expression data
    #feat = preSelectFeatures(X, y, threshold, exp_name)
    #print('{} proteins were found to have significant positive or negative correlation with the annotations.'.format(len(feat)))
    feat = selectFeaturesByWilcoxon(X, y, threshold, exp_name)
    print('{} proteins were found to be differentially expressed between the two groups.'.format(len(feat)))
    X = X[feat]
    #print('{} proteins were found to have significant positive or negative correlation with the annotations.'.format(len(feat)))
    X = protein2gene(X, X.columns)
    #X = protein2gene(X, X.columns, protein2gene_mapping)
    #evaluateClassifiers(X,y)
    selFeatures = selectFeatures([X,y], classifiers, exp_name, num_features)
    X = X[selFeatures]
    differentialPlot(X, labels.values, exp_name)
    X['Drug_Class']=labels
    #X.to_csv(os.path.join('Results/')+exp_name+'DRP_ML_selFeatures_with_annotations.csv')
    return selFeatures

def classify_transcriptomics(data, drug_data, exp_name, drugOfInterest, classifiers, num_features, threshold, preselect_choice):

    drug_data = drug_data.loc[data.index] #select samples based on data i.e., (T-ALL or B-ALL)
    #print(drug_data.shape)
    drug_data = drug_data.loc[drug_data['Class']!='Intermediate'] #filter out all sample with 'Intermediate' class
    data = data.loc[drug_data.index] #filter out samples corresponding to the 'Intermediate' class from the data
    labels = drug_data['Class']
    
    le = LabelEncoder()
    y = le.fit_transform(labels)
    cols = data.columns.astype(str)
    samples = data.index
    X = pd.DataFrame(np.array(data, dtype=float))
    
    X.columns = cols
    X.index = samples

    if preselect_choice == 'Correlation':
        #preselect features based on correlation with target variable from entire protein expression data
        feat = preSelectFeatures(X, y, threshold, exp_name)
        print('{} proteins were found to have significant positive or negative correlation with the annotations.'.format(len(feat)))
    elif preselect_choice == 'Wilcoxon':
        #preselect features based on non-parametric test (Wilcoxon)
        feat = selectFeaturesByWilcoxon(X, y, threshold, exp_name)
        print('{} proteins were found to be differentially expressed between the two groups.'.format(len(feat)))
    X = X[feat]
    #X = protein2gene(X, X.columns)
    #X = protein2gene(X, X.columns, protein2gene_mapping)
    #evaluateClassifiers(X,y)
    selFeatures = selectFeatures([X,y], classifiers, exp_name, num_features)
    X = X[selFeatures]
    differentialPlot(X, labels.values, exp_name)
    X['Drug_Class']=labels
    #X.to_csv(os.path.join('Results/')+exp_name+'DRP_ML_selFeatures_with_annotations.csv')
    return selFeatures



#if omics_type=='Transcriptomics':
#    uploaded_file = st.file_uploader("Upload transcriptomics data")
#    if uploaded_file:
#        rna = pd.read_csv(uploaded_file, index_col=0)
#        st.write("Uploaded transcriptomics data:")
#        st.dataframe(rna)
#elif omics_type=='Proteomics':
#    uploaded_file = st.file_uploader("Upload proteomics data")
#    if uploaded_file:
#        protein = pd.read_csv(uploaded_file, header=0, sep='\t', low_memory=False)
#        protein = protein.iloc[5:,:]
#        st.write("Uploaded proteomics data:")
#        st.dataframe(protein)

omics_type = st.selectbox('Select omics-type: ', ['Proteomics', 'Transcriptomics'])
cell_type = st.selectbox('Select cell-type: ', ['T-ALL', 'B-ALL'])
drugs_of_interest = ['Idarubicin', 'Dasatinib', 'Ponatinib', 'Venetoclax', 'Navitoclax', 'Doxorubicin', 'Birinapant', 'Bortezomib', 'CB-103', 'Dexamethasone', 'Cytarabine', 'Etoposide', 'Methotrexate', 'Selinexor', 'Vincristine', 'Nilotinib', 'Temsirolimus', 'Bosutinib', 'Panobinostat', 'Trametinib', 'Ruxolitinib', 'Dinaciclib', 'A1331852', 'S-63845', 'Nelarabine']
drugOfInterest = st.selectbox('Select drug', options=[opt.strip() for opt in drugs_of_interest])
num_features = st.slider('Select number of genes you want to select',1, 100, 50)
preselect_choice = st.selectbox('Select choice for dimension reduction: ', ['Wilcoxon', 'Correlation'])
threshold = st.slider('Select threshold for correlation-based feature pre-selection', 0.00, 1.00, 0.55) #threshold for correlation-based preselection
classifiers = st.multiselect('Select models - You may choose multiple among the following: [Logistic Regression, Decision Tree Classifier, Random Forest Classifier, Support Vector Machine Classifer, XG Boost Classifier and Lasso Regression]', ['LR', 'DT', 'RF', 'SVC', 'XGB', 'Lasso'])
#st.write(classifiers)

#if uploaded_file is not None:
#else:
    #st.warning("Please upload a file to proceed!")

if omics_type == 'Proteomics':
    data = mapping_Proteomics_DRP_to_metadata(drugOfInterest)
    if cell_type == 'T-ALL':
        drug_protein_df = data[1]
    elif cell_type == 'B-ALL':
        drug_protein_df = data[2]
    drug_data = data[0]
    analyze = st.button('Analyze', on_click=set_stage, args=(1,))
    if analyze:
        if len(classifiers) < 2:
            st.write('Please select at least 2 classifiers')
        else:
            #st.write(st.session_state)
            exp_name = cell_type+'_'+omics_type+'_'+drugOfInterest+'_'
            selFeatures = classify_proteomics(drug_protein_df, drug_data, exp_name, drugOfInterest, classifiers, num_features, threshold)

elif omics_type == 'Transcriptomics':
    data = mapping_Transcriptomics_DRP_to_metadata(drugOfInterest)
    drug_data = data[0]
    if cell_type == 'T-ALL':
        data = data[2]
    elif cell_type == 'B-ALL':
        data = data[4]
    preselect_choice = 'Wilcoxon'
    analyze = st.button('Analyze', on_click=set_stage, args=(1,))
    if analyze:
        if len(classifiers) < 2:
            st.write('Please select at least 2 classifiers')
        else:
            #st.write(st.session_state)
            exp_name = cell_type+'_'+omics_type+'_'+drugOfInterest+'_'
            selFeatures = classify_transcriptomics(data, drug_data, exp_name, drugOfInterest, classifiers, num_features, threshold, preselect_choice)
    
"""    
elif omics_type == 'Transcriptomics':

data = mapping_Proteomics_DRP_to_metadata(drugOfInterest)
    
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
"""
