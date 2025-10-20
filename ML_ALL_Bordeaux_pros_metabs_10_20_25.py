import numpy as np
import pandas as pd
from xgboost import XGBClassifier, plot_importance
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score,  RandomizedSearchCV, LeaveOneOut, cross_validate, GridSearchCV
from sklearn.metrics import balanced_accuracy_score, classification_report, confusion_matrix, roc_auc_score, matthews_corrcoef, make_scorer
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC, LinearSVC
import  umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def nestedCV(pipe,param_dist,cv_inner, cv_outer,X,y):
    '''
    Runs NestedCV to ensure that we are getting unbiased estimates of the performance of the models. 

    #----------
    #Inputs
    #-----------
    pipe - pipeline that will be run through GridSearchCV and cross validated
    param_grid - distribution of parameters that will be considered in the analysis
    cv_inner - the schema/folds that will be used for internal cross validation
    cv_outer - the schema/folds that will be used for exernal/outer cross validation
    X - the predictors
    y - the outcomes to train against.

    #----------
    # Outputs
    #----------
    scores - this outputs a dataframe with stored scores for Matthew's Correlation Coefficient, Balanced Accuracy Score and AUROC

    '''
    # Randomized search that will run INSIDE each outer fold
    search = GridSearchCV(
        estimator=pipe,                 # your Pipeline
        param_grid=param_dist,
        scoring="roc_auc",    #outer scores computed separately below
        cv=cv_inner,                       # inner CV for tuning
        n_jobs=64,
        verbose=1,
        refit=True
    )

    scoring = {
        "ba": make_scorer(balanced_accuracy_score),
        "auc": "roc_auc",
        "mcc": make_scorer(matthews_corrcoef)
    }

    scores = cross_validate(search, X, y, cv=cv_outer, scoring=scoring, n_jobs=64)
    
    return scores

# ---------------------------------------------
# 1. Reading in the data
# ---------------------------------------------

# Example:
X = pd.read_csv("/VOSNE/Fischbein/bhislop/metabolomics/Merged_path_added.csv")#read in the data

# ---------------------------------------------
# 2. Filter out the post data, and make sure the data only has Aneurysms for now
# ---------------------------------------------
#
X_metabs = pd.read_csv("/Users/bradyhislop/Library/CloudStorage/Box-Box/Brady Hislop's Files/HA_studies/Metabolomics/Analysis_spreadsheets/merged/Merged_path_added.csv")#read in the data
X_pros = pd.read_csv("/Users/bradyhislop/Library/CloudStorage/Box-Box/Brady Hislop's Files/HA_studies/Olink_Proteomics/Analysis_spreadsheets/proteins_plsda_10_16.csv")#read in the data

X_metabs = X_metabs.loc[
    (X_metabs["prePost"]=="pre")&(X_metabs["Surgery"]=="Aneurysm"),
]
X_metabs = X_metabs.drop(["prePost","Surgery"],axis=1)
# scaler = StandardScaler()

# X_metabs.iloc[:,8:] = scaler.fit_transform(X_metabs.iloc[:,8:])
X_metabs_scaling_cols = X_metabs.iloc[:,8:].columns
passthrough_cols = X_pros.iloc[:,8:].columns

X_pros = X_pros.loc[
    X_pros["Grades"].isin(["Mild","Mild-to-moderate","Moderate","Moderate-to-severe"]),
]

X_metabs["Grades"].replace({
    "Moderate-to-severe": "Moderate",
    "Mild-to-moderate": "Mild"
}, inplace=True)

X_pros["Grades"].replace({
    "Moderate-to-severe": "Moderate",
    "Mild-to-moderate": "Mild"
}, inplace=True)


X_pros = X_pros.set_index('SampleID')
X_metabs = X_metabs.set_index('SampleID')


# X_metabs.join(X_pros,on="Sex",how="inner")

X = pd.merge(X_metabs, X_pros.drop("Surgery",axis=1), on=["Age","Sex","Diameter","Location","Syndrome","Grades","Size"], how="inner")
# X_grades.to_csv("/Users/bradyhislop/Downloads/SurgeryWhy.csv")
y = X["Grades"]

X_grades = X.iloc[:,list(range(0,2))+list(range(6,7))+list(range(3,4))+list(range(4,5))+list(range(7,X.shape[1]))]
categorical_cols = ["Sex","Syndrome","Location","Size"]
X_grades[categorical_cols] = X_grades[categorical_cols].astype("category")
scaler = StandardScaler()
X_grades["Age"] = scaler.fit_transform(X_grades[["Age"]])
X_grades = pd.get_dummies(X_grades,columns=categorical_cols,drop_first=True)
le = LabelEncoder()
y=le.fit_transform(y)
X_metabs_scaling_cols = X_metabs.iloc[:,8:].columns
passthrough_cols = ["Age","Sex_Male","Syndrome_MFS","Syndrome_NS","Location_Root","Size_>5"] + list(X_pros.iloc[:,8:].columns)

# ---------------------------------------------
# 2.1 Set columns that will be passed through the standard scaler without transformation.
# ---------------------------------------------
#
# passthrough_cols = ["Age","Sex_Male","Syndrome_MFS","Syndrome_NS","Location_Root","Size_>5"]

# ---------------------------------------------
# 2.2 Set up the Column Transformer
# ---------------------------------------------
#
column_transformer = ColumnTransformer([
 ('scale',StandardScaler(), X_metabs_scaling_cols),
 ('pass','passthrough',passthrough_cols)
])


# ---------------------------------------------
# 3. Set up an outcomes Pandas DataFrame from which we can save all the model results.
# ---------------------------------------------
#
outcomes = pd.DataFrame(np.zeros((16,6)),columns=["BA_mean","BA_std","AUC_mean","AUC_std","MCC_mean","MCC_std"],index=['RFC','XGB','LGR','SVC','RFC_PCA','RFC_UMAP',
                                                                                                                       'RFC_tSNE','XGB_PCA','XGB_UMAP','XGB_tSNE',
                                                                                                                       'LGR_PCA','LGR_UMAP','LGR_tSNE',
                                                                                                                       'SVC_PCA','SVC_UMAP','SVC_tSNE'])


# ---------------------------------------------
# 4. Train Random Forest Classifier (no pre-variable selection)
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('model', RandomForestClassifier(
    random_state=42,
    class_weight="balanced_subsample",  # good for class imbalance
    n_jobs=64))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "model__n_estimators": np.linspace(300, 1200, 10, dtype=int),
    "model__max_depth": [None, 6, 8, 10, 12, 16],
    "model__min_samples_split": [2, 5, 10],
    "model__min_samples_leaf": [1, 2, 4],
    "model__max_features": ["sqrt", "log2", None],
    "model__bootstrap": [True, False]
}
scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)


#Save the outcomes from the nested CV
outcomes.loc["RFC"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 5. Train XGB Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('model', XGBClassifier(
    eval_metric="auc",
    tree_method="hist",
    enable_categorical=True,
    random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "model__n_estimators":      np.linspace(200, 1200, 6, dtype=int),
    "model__learning_rate":     np.logspace(np.log10(0.005), np.log10(0.2), 8),
    "model__max_depth":         [3, 4, 5, 6, 7, 8],
    "model__subsample":         np.linspace(0.1,0.8,8),
    }

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["XGB"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 6. Nested Cross Validation of Logistic Regression Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('model', LogisticRegression(penalty="elasticnet",solver="saga",class_weight="balanced",max_iter=5000,random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "model__C": np.logspace(-3, 1, 10),        # inverse regularization strength (smaller = stronger)
    "model__l1_ratio": np.linspace(0, 1, 6)    # 0=L2, 1=L1, values in-between = elastic net
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["LGR"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 7. Nested cross validation analysis of Support Vector Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('model', SVC(probability=True))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    'model__C': np.logspace(-3, 3, 20),
    'model__gamma': np.logspace(-4, 1, 20),
    'model__kernel': ['rbf', 'poly', 'sigmoid','linear']
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["SVC"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 8. Nested cross validation analysis of Random Forest Classifier with PCA before to reduce dimensionality
# ---------------------------------------------
#
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('model', RandomForestClassifier(
    random_state=42,
    class_weight="balanced_subsample",  # good for class imbalance
    n_jobs=64))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "model__n_estimators": np.linspace(300, 1200, 10, dtype=int),
    "model__max_depth": [None, 6, 8, 10, 12, 16],
    "model__min_samples_split": [2, 5, 10],
    "model__min_samples_leaf": [1, 2, 4],
    "model__max_features": ["sqrt", "log2", None],
    "model__bootstrap": [True, False]
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["RFC_PCA"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 10. Nested cross validation analysis of Random Forest Classifier with UMAP before to reduce dimensionality
# ---------------------------------------------
#
pipe = Pipeline([
    ('pre', column_transformer),
    ('umap',umap.UMAP(random_state=42)),
    ('model', RandomForestClassifier(
    random_state=42,
    class_weight="balanced_subsample",  # good for class imbalance
    n_jobs=64))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "umap__n_components": np.linspace(1,40,39,dtype=int),
    "model__n_estimators": np.linspace(300, 1200, 10, dtype=int),
    "model__max_depth": [None, 6, 8, 10, 12, 16],
    "model__min_samples_split": [2, 5, 10],
    "model__min_samples_leaf": [1, 2, 4],
    "model__max_features": ["sqrt", "log2", None],
    "model__bootstrap": [True, False]
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["RFC_UMAP"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 11. Nested cross validation analysis of Random Forest Classifier with tSNE before to reduce dimensionality
# ---------------------------------------------
#
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('tsne',TSNE(init="pca", learning_rate="auto", metric='cosine', random_state=42)),
    ('model', RandomForestClassifier(
    random_state=42,
    class_weight="balanced_subsample",  # good for class imbalance
    n_jobs=64))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__perplexity": np.linspace(1,100,99,dtype=int),
    "model__n_estimators": np.linspace(300, 1200, 10, dtype=int),
    "model__max_depth": [None, 6, 8, 10, 12, 16],
    "model__min_samples_split": [2, 5, 10],
    "model__min_samples_leaf": [1, 2, 4],
    "model__max_features": ["sqrt", "log2", None],
    "model__bootstrap": [True, False]
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["RFC_tSNE"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 12. Nested cross validation analysis of XGB Classifier with PCA before to reduce dimensionality
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('model', XGBClassifier(
    eval_metric="auc",
    tree_method="hist",
    enable_categorical=True,
    random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components":        np.linspace(1,40,39,dtype=int),
    "model__n_estimators":      np.linspace(200, 1200, 6, dtype=int),
    "model__learning_rate":     np.logspace(np.log10(0.005), np.log10(0.2), 8),
    "model__max_depth":         [3, 4, 5, 6, 7, 8],
    "model__subsample":         np.linspace(0.1,0.8,8),
    }

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["XGB_PCA"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 12. Nested cross validation analysis of XGB Classifier with UMAP before to reduce dimensionality
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('umap',umap.UMAP(random_state=42)),
    ('model', XGBClassifier(
    eval_metric="auc",
    tree_method="hist",
    enable_categorical=True,
    random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "umap__n_components":        np.linspace(1,40,39,dtype=int),
    "model__n_estimators":      np.linspace(200, 1200, 6, dtype=int),
    "model__learning_rate":     np.logspace(np.log10(0.005), np.log10(0.2), 8),
    "model__max_depth":         [3, 4, 5, 6, 7, 8],
    "model__subsample":         np.linspace(0.1,0.8,8),
    }

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["XGB_UMAP"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 12. Nested cross validation analysis of Random Forest Classifier with tSNE before to reduce dimensionality
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('tsne',TSNE(init="pca", learning_rate="auto", metric='cosine', random_state=42)),
    ('model', XGBClassifier(
    eval_metric="auc",
    tree_method="hist",
    enable_categorical=True,
    random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__perplexity": np.linspace(1,100,99,dtype=int),
    "model__n_estimators":      np.linspace(200, 1200, 6, dtype=int),
    "model__learning_rate":     np.logspace(np.log10(0.005), np.log10(0.2), 8),
    "model__max_depth":         [3, 4, 5, 6, 7, 8],
    "model__subsample":         np.linspace(0.1,0.8,8),
    }

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["XGB_tSNE"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 13. Nested Cross Validation of Logistic Regression Classifier with PCA
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('model', LogisticRegression(penalty="elasticnet",solver="saga",class_weight="balanced",max_iter=5000,random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "model__C": np.logspace(-3, 1, 10),        # inverse regularization strength (smaller = stronger)
    "model__l1_ratio": np.linspace(0, 1, 6)    # 0=L2, 1=L1, values in-between = elastic net
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["LGR_PCA"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 14. Nested Cross Validation of Logistic Regression Classifier with UMAP
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('umap',umap.UMAP(random_state=42)),
    ('model', LogisticRegression(penalty="elasticnet",solver="saga",class_weight="balanced",max_iter=5000,random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "umap__n_components": np.linspace(1,40,39,dtype=int),
    "model__C": np.logspace(-3, 1, 10),        # inverse regularization strength (smaller = stronger)
    "model__l1_ratio": np.linspace(0, 1, 6)    # 0=L2, 1=L1, values in-between = elastic net
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["LGR_UMAP"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 14. Nested Cross Validation of Logistic Regression Classifier with UMAP
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('tsne',TSNE(init="pca", learning_rate="auto", metric='cosine', random_state=42)),
    ('model', LogisticRegression(penalty="elasticnet",solver="saga",class_weight="balanced",max_iter=5000,random_state=42))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__perplexity": np.linspace(1,100,99,dtype=int),
    "model__C": np.logspace(-3, 1, 10),        # inverse regularization strength (smaller = stronger)
    "model__l1_ratio": np.linspace(0, 1, 6)    # 0=L2, 1=L1, values in-between = elastic net
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["LGR_tSNE"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]


# ---------------------------------------------
# 7. Nested cross validation analysis of Support Vector Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('model', SVC(probability=True))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    'pca__n_components': np.linspace(1,40,39,dtype=int),
    'model__C': np.logspace(-3, 3, 20),
    'model__gamma': np.logspace(-4, 1, 20),
    'model__kernel': ['rbf', 'poly', 'sigmoid','linear']
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["SVC_PCA"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 7. Nested cross validation analysis of Support Vector Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('umap',umap.UMAP(random_state=42)),
    ('model', SVC(probability=True))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    'umap__n_components': np.linspace(1,40,39,dtype=int),
    'model__C': np.logspace(-3, 3, 20),
    'model__gamma': np.logspace(-4, 1, 20),
    'model__kernel': ['rbf', 'poly', 'sigmoid','linear']
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["SVC_UMAP"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

# ---------------------------------------------
# 7. Nested cross validation analysis of Support Vector Classifier
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('pre', column_transformer),
    ('pca',PCA(random_state=42)),
    ('tsne',TSNE(init="pca", learning_rate="auto", metric='cosine', random_state=42)),
    ('model', SVC(probability=True))
])

# ---- Search space (prefix with 'model__' since RF is under 'model' in the Pipeline)
param_dist = {
    "pca__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__n_components": np.linspace(1,40,39,dtype=int),
    "tsne__perplexity": np.linspace(1,100,99,dtype=int),
    'model__C': np.logspace(-3, 3, 20),
    'model__gamma': np.logspace(-4, 1, 20),
    'model__kernel': ['rbf', 'poly', 'sigmoid','linear']
}

scores = nestedCV(pipe,param_dist,inner,outer,X_grades,y)
outcomes.loc["SVC_UMAP"] = [scores['test_ba'].mean(), 
                       scores['test_ba'].std(),
                       scores['test_auc'].mean(),
                       scores['test_auc'].std(),
                       scores['test_mcc'].mean(),
                       scores['test_mcc'].std()]

outcomes.to_csv("/VOSNE/Fischbein/bhislop/metabolomics/model_outcomes.csv",index=False)