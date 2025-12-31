from xgboost import XGBClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import numpy as np

def calculate_metrics(true_labels, pred_labels):
    """
    Calculate sensitivity and specificity.
    :param true_labels: 0/1
    :param pred_labels: 0/1
    :return: (sensitivity, specificity)
    """
    if len(true_labels) != len(pred_labels):
        raise ValueError("the length of the labels were not matched")

    tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(true_labels, pred_labels) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(true_labels, pred_labels) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == 1 and p == 0)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

    return sensitivity, specificity

def summary(all_pred_probs, all_pred_labels, all_true_labels):
    # Convert to numpy arrays for scikit-learn
    all_pred_probs = np.array(all_pred_probs)
    all_true_labels = np.array(all_true_labels)

    fpr, tpr, _ = roc_curve(all_true_labels, all_pred_probs)
    roc_auc = auc(fpr, tpr)  # Calculate AUC score
    sensitivity, specificity = calculate_metrics(all_true_labels, all_pred_labels)
    return fpr, tpr, roc_auc, sensitivity, specificity

if __name__=="__main__":
    tpm_df = stad_tpm_df[stad_pheno_df["run_accession"]].T.reset_index(drop=True)
    pheno_df = stad_pheno_df

    y_true_ls = []
    y_pred_ls = []
    y_prob_ls = []
    df_ls = []
    kf = KFold(n_splits=10, random_state=None, shuffle=False)
    for i, (train_index, test_index) in enumerate(kf.split(pheno_df)):
        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss')
        train_tpm_df = tpm_df.loc[train_index]
        test_tpm_df = tpm_df.loc[test_index]
        train_df = pheno_df.loc[train_index]
        test_df = pheno_df.loc[test_index]
        model.fit(train_tpm_df, train_df["Label"])
        ## calculate feature importance
        importance_scores = model.feature_importances_
        feature_names = train_tpm_df.columns
        feature_importance_df = pd.DataFrame({
            'Feature': feature_names,
            f'Importance_{i}': importance_scores
        })
        feature_importance_df.set_index('Feature', inplace=True)
        df_ls.append(feature_importance_df)
        y_pred = model.predict(test_tpm_df)
        y_prob = model.predict_proba(test_tpm_df)[:,1]
        y_pred_ls.extend(y_pred)
        y_true_ls.extend(test_df["Label"].tolist())
        y_prob_ls.extend(y_prob)

    fpr, tpr, roc_auc, sensitivity, specificity = summary(y_prob_ls, y_pred_ls, y_true_ls)
    feature_mean = pd.DataFrame(pd.concat(df_ls, axis=1).mean(axis=1), columns=["importance"])
