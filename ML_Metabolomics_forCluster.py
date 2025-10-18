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

def nestedCV(pipe,param_dist,cv_inner, cv_outer,X,y):

    # Randomized search that will run INSIDE each outer fold
    search = GridSearchCV(
        estimator=pipe,                 # your Pipeline
        param_grid=param_dist,
        scoring="roc_auc",    #outer scores computed separately below
        cv=cv_inner,                       # inner CV for tuning
        n_jobs=-1,
        verbose=1,
        refit=True
    )

    scoring = {
        "ba": make_scorer(balanced_accuracy_score),
        "auc": "roc_auc",
        "mcc": make_scorer(matthews_corrcoef)
    }

    scores = cross_validate(search, X, y, cv=cv_outer, scoring=scoring, n_jobs=-1)
    
    return scores

# Example:
X = pd.read_csv("/Users/bradyhislop/Library/CloudStorage/Box-Box/Brady Hislop's Files/HA_studies/Metabolomics/Analysis_spreadsheets/merged/Merged_path_added.csv")#read in the data

#filter out the post data
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
y = X_grades["Grades"]

X_grades = X_grades.iloc[:,list(range(1,3))+list(range(6,7))+list(range(4,5))+list(range(8,9))+list(range(9,X_grades.shape[1]))]
categorical_cols = ["Sex","Syndrome","Location","Size"]
X_grades[categorical_cols] = X_grades[categorical_cols].astype("category")
scaler = StandardScaler()
X_grades["Age"] = scaler.fit_transform(X_grades[["Age"]])
X_grades = pd.get_dummies(X_grades,columns=categorical_cols,drop_first=True)
le = LabelEncoder()
y=le.fit_transform(y)


inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=2)
# --- Pipeline to improve model performance evaluation
pipe = Pipeline([
    ('model', RandomForestClassifier(
    random_state=42,
    class_weight="balanced_subsample",  # good for class imbalance
    n_jobs=-1))
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

print(f"Nested BA : {scores['test_ba'].mean():.3f} ± {scores['test_ba'].std():.3f}")
print(f"Nested AUC : {scores['test_auc'].mean():.3f} ± {scores['test_auc'].std():.3f}")
print(f"Nested MCC: {scores['test_mcc'].mean():.3f} ± {scores['test_mcc'].std():.3f}")