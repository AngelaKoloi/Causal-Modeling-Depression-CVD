import pandas as pd
import numpy as np


def get_characteristic(df: pd.DataFrame) -> pd.DataFrame:
    df = df[[col for col in df.columns if col.endswith('_CO')]]
    percentage = 100

    female_df = df[df["sex_2007_CO"] == 1]
    male_df = df[df["sex_2007_CO"] == 2]

    female_row_count = female_df.shape[0]
    male_row_count = male_df.shape[0]

    return (pd.DataFrame({
        "N": [female_row_count, male_row_count],
        "Age (mean) in 2007": [female_df["age_2007_CO"].mean(), male_df["age_2007_CO"].mean()],
        "Age (sd) in 2007": [female_df["age_2007_CO"].std(), male_df["age_2007_CO"].std()],
        "Smoking (NO) in 2007": [
            (100 - (female_df["smoking_2007_CO"].sum() /
             female_row_count) * percentage),
            (100 - (male_df["smoking_2007_CO"].sum() /
             male_row_count) * percentage)
        ],
        "Smoking (NO) in 2011": [
            (100 - (female_df["smoking_2011_CO"].sum() /
             female_row_count) * percentage),
            (100 - (male_df["smoking_2011_CO"].sum() /
             male_row_count) * percentage)
        ],
        "Smoking (NO) in 2012": [
            (100 - (female_df["smoking_2012_CO"].sum() /
             female_row_count) * percentage),
            (100 - (male_df["smoking_2012_CO"].sum() /
             male_row_count) * percentage)
        ],
        "BMI (mean) in 2007": [female_df["bmi_2007_CO"].mean(), male_df["bmi_2007_CO"].mean()],
        "BMI (sd) in 2007": [female_df["bmi_2007_CO"].std(), male_df["bmi_2007_CO"].std()],
        "BMI (mean) in 2011": [female_df["bmi_2011_CO"].mean(), male_df["bmi_2011_CO"].mean()],
        "BMI (sd) in 2011": [female_df["bmi_2011_CO"].std(), male_df["bmi_2011_CO"].std()],
        "BMI (mean) in 2012": [female_df["bmi_2012_CO"].mean(), male_df["bmi_2012_CO"].mean()],
        "BMI (sd) in 2012": [female_df["bmi_2012_CO"].std(), male_df["bmi_2012_CO"].std()],
    },
        index=['Female', 'Male']
    ).round(2).T)


def descriptive_cvd(df: pd.DataFrame, round: int = 2) -> pd.DataFrame:
    df = df[[col for col in df.columns if col.endswith('_CVD')]]

    descriptives = {}
    for col in df.columns:
        col_data = df[col]
        Min = col_data.min()
        Max = col_data.max()
        Mean = col_data.mean()
        SD = col_data.std()

        descriptives[col] = {
            "Min.": Min,
            "Max.": Max,
            "Mean": Mean,
            "SD": SD
        }
    return (pd.DataFrame(descriptives).round(round).T)


def descriptive_depression(df: pd.DataFrame) -> pd.DataFrame:
    percentage = 100
    df = df[[col for col in df.columns if col.endswith('_D')]]
    df = df.round()
    df_row_count = df.shape[0]

    descriptives = {}
    for col in df.columns:
        col_data = df[col]
        value_counts = col_data.value_counts().reset_index()
        try:
            zero_count = value_counts[value_counts[col] == 0]['count'].iloc[0]
        except (IndexError) as e:
            zero_count = 0

        try:
            one_count = value_counts[value_counts[col] == 1]['count'].iloc[0]
        except (IndexError) as e:
            one_count = 0

        try:
            two_count = value_counts[value_counts[col] == 2]['count'].iloc[0]
        except (IndexError) as e:
            two_count = 0

        try:
            three_count = value_counts[value_counts[col] == 3]['count'].iloc[0]
        except (IndexError) as e:
            three_count = 0

        descriptives[col] = {
            "Max value percentage": (zero_count / df_row_count) * percentage,
            "Min value percentage": (three_count / df_row_count) * percentage,
            "Other": ((df_row_count - zero_count - three_count) / df_row_count) * percentage,

            # "Percentage 1": (one_count / df_row_count) * percentage,
            # "Percentage 2": (two_count / df_row_count) * percentage,
        }
    return (pd.DataFrame(descriptives).round(2).T)
