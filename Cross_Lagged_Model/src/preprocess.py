import pandas as pd
import statsmodels.api as sm
import re

from functools import reduce
from column_names import new_column_names, sorted_column_names


def recode_values(value: int) -> int:
    if value == 0:
        return 0
    elif value in [2, 3, 11, 12]:
        return 1
    elif value in [4, 5, 21, 22]:
        return 2
    elif value in [6, 7, 31, 32]:
        return 3
    else:
        return value


def process_depressive_data(depressive_data: pd.DataFrame) -> pd.DataFrame:
    depressive_data.rename(
        columns={depressive_data.columns[0]: 'patientID'}, inplace=True)
    depressive_data['patientID'] = depressive_data['patientID'].astype(int)

    columns_to_recode = ['b16', 'b18', 'b16_12',
                         'b18_12', 'deprb1611', 'deprb1811']
    for col in columns_to_recode:
        depressive_data[col] = depressive_data[col].apply(recode_values)

    columns = ['deprb0111', 'deprb0211', 'deprb0311', 'deprb0411', 'deprb0511', 'deprb0611', 'deprb0711',
               'deprb0811', 'deprb0911', 'deprb1011', 'deprb1111', 'deprb1211', 'deprb1311', 'deprb1411',
               'deprb1511', 'deprb1711', 'deprb1911', 'deprb2011', 'deprb2111']
    depressive_data[columns] = depressive_data[columns].apply(lambda x: x - 1)

    return depressive_data


def process_cvd_data(cvd_path: str) -> pd.DataFrame:
    cvd_excel = pd.ExcelFile(cvd_path)
    cvd_sheets = []

    for sheet in cvd_excel.sheet_names:
        sheet = pd.read_excel(cvd_path, sheet_name=sheet)
        cvd_sheets.append(sheet)

    for dataframe in cvd_sheets:
        dataframe.rename(
            columns={dataframe.columns[0]: 'patientID'}, inplace=True)

    cvd_data = reduce(lambda left, right: pd.merge(left, right, on=['patientID'],
                                                   how='outer'), cvd_sheets)
    return cvd_data


def preprocess_data(depressive_data: pd.DataFrame, cvd_data: pd.DataFrame) -> pd.DataFrame:
    depressive_data = process_depressive_data(depressive_data)
    merged_data = pd.merge(depressive_data, cvd_data, on='patientID')
    merged_data = merged_data.rename(columns=new_column_names)

    merged_data = merged_data.interpolate()
    return merged_data


def calculate_residu(data, year):
    print(f"year: {year}")

    covariance_var = [
        col for col in data if col.endswith("CO")] + ["patiendID"]
    print(f"len covariance list: {len(covariance_var)}")
    outcome_vars = list(set(data.columns) - set(covariance_var))

    print(f"len outcome list: {len(outcome_vars)}\n")
    X = data[[f'sex_{year}_CO', f'age_{year}_CO']]
    X = sm.add_constant(X)

    outcome_list = []
    for outcome in outcome_vars:
        y = data[outcome].values
        X_values = X.values
        model = sm.OLS(y, X_values)

        results = model.fit()
        residuals = results.resid
        outcome_list.append(pd.Series(residuals, name=outcome))

    outcome_list = pd.concat(outcome_list, axis=1)
    return outcome_list


def residu_data(data: pd.DataFrame) -> pd.DataFrame:
    years = list({re.search(r'_(\d{4})_', col).group(1) for col in data.columns
                  if re.search(r'_(\d{4})_', col)})
    years = sorted(years, key=lambda x: int(x))
    print(f"years in df: {years}\n")

    residu_data_list = [data["patientID"]]

    for year in years:
        columns_with_year = [
            col for col in data.columns if re.search(fr'_{year}_', col)]
        temp_data = data[columns_with_year]

        residu_data_list.append(calculate_residu(temp_data, year))

    residu_data = pd.concat(residu_data_list, axis=1)
    residu_data.head()

    sort_list = ["patientID"] + list(sorted_column_names.values())
    residu_data = residu_data[sort_list]
    return residu_data
