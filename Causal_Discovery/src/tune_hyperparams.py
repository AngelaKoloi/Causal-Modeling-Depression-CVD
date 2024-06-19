import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import f1_score
from xgboost import XGBClassifier
import optuna
import joblib

data = pd.read_csv('../data/preprocessed/final_data_wave1.csv')

X = data.drop(['acidep09', 'amet_syn2'], axis=1)
y = data['acidep09']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

def objective(trial):
    param = {
        'n_estimators': trial.suggest_int('n_estimators', 50, 300),
        'max_depth': trial.suggest_int('max_depth', 3, 7),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.2),
        'subsample': trial.suggest_float('subsample', 0.5, 1.0),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
        'missing': np.nan,
        'use_label_encoder': False,
        'eval_metric': 'logloss'
    }

    model = XGBClassifier(**param)

    kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_f1 = cross_val_score(model, X_train, y_train, cv=kfold, scoring='f1_weighted')

    return cv_f1.mean()

study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=1000)
best_params = study.best_params
print(f'Best hyperparameters: {best_params}')