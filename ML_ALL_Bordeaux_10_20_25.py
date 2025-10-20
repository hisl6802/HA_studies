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
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC,LinearSVC
import joblib
from sklearn.decomposition import PCA
import umap
from sklearn.manifold import TSNE


def nestedCV(pipe,param_dist,cv_inner, cv_outer,X,y):

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
X = X.loc[
    (X["prePost"]=='pre') & (X["Surgery"]=="Aneurysm"),
]
X = X.drop("prePost",axis=1)

X_grades = X.loc[
    X["Grades"].isin(["Mild","Mild-to-moderate","Moderate","Moderate-to-severe"]),
]

X_grades["Grades"].replace({
    "Moderate-to-severe": "Moderate",
    "Mild-to-moderate": "Mild"
}, inplace=True)
y = X_grades["Grades"] #extract the labels for training the model

X_grades = X_grades.iloc[:,list(range(1,3))+list(range(6,7))+list(range(4,5))+list(range(8,9))+list(range(9,X_grades.shape[1]))]
categorical_cols = ["Sex","Syndrome","Location","Size"]
X_grades[categorical_cols] = X_grades[categorical_cols].astype("category")
scaler = StandardScaler()
X_grades["Age"] = scaler.fit_transform(X_grades[["Age"]])
X_grades = pd.get_dummies(X_grades,columns=categorical_cols,drop_first=True)

#Setting the strategy for the K-folds making sure they are stratified to ensure good fits and that I have an inner and outer CV set.
inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=2)

#Label Encoder for the Mild and Moderate grades
le = LabelEncoder()
y  = le.fit_transform(y)

# ---------------------------------------------
# 3. Set up an outcomes Pandas DataFrame from which we can save all the model results.
# ---------------------------------------------
#
outcomes = pd.DataFrame(np.zeros((10,6)),columns=["BA_mean","BA_std","AUC_mean","AUC_std","MCC_mean","MCC_std"],index=['RFC','XGB','LGR','SVC','RFC_PCA','RFC_UMAP','RFC_tSNE','XGB_PCA','XGB_UMAP','XGB_tSNE'])


# ---------------------------------------------
# 4. Train Random Forest Classifier (no pre-variable selection)
# ---------------------------------------------
#
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
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