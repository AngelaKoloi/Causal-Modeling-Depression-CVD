import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from sklearn.impute import KNNImputer
from sklearn.feature_selection import mutual_info_regression

class DataAnalysis:
    @staticmethod
    def plot_missing_values(df):
        """Plots the number of missing values per feature in a DataFrame as a bar chart."""
        sum_missing_values = df.isnull().sum()
        total_missing_values = sum_missing_values[sum_missing_values > 0]

        if len(total_missing_values) > 0:
            plt.figure(figsize=(12, 6))
            total_missing_values.plot(kind='bar')
            plt.title('Columns with Missing Data')
            plt.xlabel('Features')
            plt.ylabel('Number of Missing Values')
            plt.xticks(rotation=45, ha='right')
            plt.show()
        else:
            print('No missing values')

    @staticmethod
    def plot_missing_heatmap(df):
        """Plots an interactive heatmap of missing values in a DataFrame using Plotly."""
        missing_data = df.isnull()

        fig = px.imshow(missing_data, 
                        labels=dict(x="Features", y="Number of Missing Values", color="Missing Data"),
                        x=df.columns, 
                        color_continuous_scale='Viridis')
        
        fig.update_layout(
            title='Heatmap of Missing Data',
            xaxis_nticks=50, 
            yaxis_nticks=10 
        )
        fig.update_xaxes(tickangle=90)
        fig.show()

    @staticmethod
    def get_list_of_features(df):
        features = df.columns
        return features.tolist()
    
    @staticmethod
    def describe_missingness(df):
        missing_data = df.columns[df.isnull().any()]
        return df[missing_data].describe()
        
    @staticmethod
    def plot_missing_data_distribution(df):
        """Plots the distribution for each column with missing data."""
        any_missing_data = df.columns[df.isnull().any()]

        num_cols = 3  
        num_rows = (len(any_missing_data) + num_cols - 1) // num_cols  # Rounds up if not a multiple of num_cols

        fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(15, num_rows * 4))
        axes = axes.flatten()  

        for i, col in enumerate(any_missing_data):
            ax = axes[i]
            df[col].dropna().plot(kind='hist', bins=30, ax=ax, color='blue', alpha=0.7)
            ax.set_title(f'Distribution of {col}')
            ax.set_xlabel(col)
            ax.set_ylabel('Frequency')

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_correlation_matrix(df, columns, title, label_dict=None):
        """
        Plots a correlation matrix for specified columns in a DataFrame.

        Parameters:
        df (pd.DataFrame): The DataFrame to analyze.
        columns (list): List of columns to include in the correlation matrix.
        title (string): Title of the plot.
        :param label_dict: dictionary with original column names as keys and new titles as values.
        """
        df_filtered = df[columns]
        corr_matrix = df_filtered.corr()

        plt.figure(figsize=(20, 20))
        sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", square=True,
                annot_kws={"size": 14}, cbar_kws={"shrink": 0.8})
        
        if label_dict:
            ticks = [label_dict.get(column, column) for column in columns]
        else:
            ticks = columns

        plt.title(f"Correlation Matrix of {title}", fontsize=20, pad=10)
        plt.xticks(ticks=range(len(ticks)), labels=ticks, rotation=45, ha='right', fontsize=14)
        plt.yticks(ticks=range(len(ticks)), labels=ticks, rotation=0, fontsize=14)
        # plt.xticks(rotation=45, ha='right', fontsize=14)
        # plt.yticks(rotation=0, fontsize=14)
        plt.subplots_adjust(top=0.9)
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_correlation_matrix_two_df(csv_file1, csv_file2, column_mapping, figsize=(12, 10), cmap='coolwarm'):
        """
        Plots a correlation heatmap between two datasets with shared variables but different column names.

        Parameters:
        - csv_file1 (str): Path to the first CSV file.
        - csv_file2 (str): Path to the second CSV file.
        - column_mapping (dict): A dictionary mapping column names from the second dataset to those in the first dataset.
        - figsize (tuple): Dimensions of the figure (width, height).
        - cmap (str): Colormap used for the heatmap visualization.
        """
        data_wave1 = pd.read_csv(csv_file1).set_index('pident')
        data_wave2 = pd.read_csv(csv_file2).set_index('pident')

        # Rename columns in the second dataset to match the first
        data_wave2_renamed = data_wave2.rename(columns=column_mapping)

        # Filter columns based on the mapping provided, ensuring columns exist in both datasets
        columns_wave1 = [column_mapping[key] for key in column_mapping if column_mapping[key] in data_wave1.columns]
        data_wave1_filtered = data_wave1[columns_wave1]
        data_wave2_filtered = data_wave2_renamed[columns_wave1]

        # Get original column names from Wave 2 that are present in the mapping
        original_wave2_labels = [key for key in column_mapping if key in data_wave2.columns]

        combined_data = pd.concat([data_wave1_filtered, data_wave2_filtered], axis=0)
        correlation_matrix = combined_data.corr()

        plt.figure(figsize=figsize)
        sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt=".2f", xticklabels=True, yticklabels=True)
        plt.xticks(ticks=np.arange(len(original_wave2_labels)) + 0.5, labels=original_wave2_labels, rotation=45, ha='right')
        plt.title('Correlation Matrix of the Shared Metabolites')
        plt.show()
    
    @staticmethod
    def plot_mutual_information_heatmap(file_path_wave1, file_path_wave2, column_mapping, n_neighbors=10):
        """
        Plots a heatmap of the mutual information matrix between two datasets.

        Parameters:
        - file_path_wave1 (str): Path to the first CSV file containing dataset.
        - file_path_wave2 (str): Path to the second CSV file containing dataset.
        - column_mapping (dict): A dictionary mapping column names from the second dataset to those in the first dataset.
        - n_neighbors (int): Number of neighbors to use for KNN imputation.
        """
        data_wave1 = pd.read_csv(file_path_wave1).set_index('pident')
        data_wave2 = pd.read_csv(file_path_wave2).set_index('pident')

        # only get the specified columns in column mapping from both datasets
        common_columns_wave1 = data_wave1.columns.intersection(list(column_mapping.values()))
        wave2_columns_mapped_to_wave1 = [column_mapping[col] for col in data_wave2.columns if col in column_mapping]
        df_wave1_filtered = data_wave1[common_columns_wave1]
        df_wave2_filtered = data_wave2[[col for col in column_mapping.keys() if col in data_wave2.columns]]
        df_wave2_filtered.columns = wave2_columns_mapped_to_wave1

        # impute missing values with the same method as what you would use during training
        df_combined = pd.concat([df_wave1_filtered, df_wave2_filtered])
        knn_imputer = KNNImputer(n_neighbors=n_neighbors)
        df_imputed = pd.DataFrame(knn_imputer.fit_transform(df_combined), columns=df_combined.columns)

        # calculate mutual information
        mi_matrix = DataAnalysis.compute_mutual_information(df_imputed)

        DataAnalysis.plot_heatmap(mi_matrix, column_mapping)

    @staticmethod
    def compute_mutual_information(df):
        mi = pd.DataFrame(index=df.columns, columns=df.columns)
        for col1 in df.columns:
            for col2 in df.columns:
                mi.loc[col1, col2] = mutual_info_regression(df[[col1]], df[[col2]].values.ravel())[0]
        return mi.astype(float)

    @staticmethod
    def plot_heatmap(mi_matrix, column_mapping):
        plt.figure(figsize=(12, 10))
        sns.heatmap(mi_matrix, annot=True, cmap='viridis', fmt=".2f")

        # reverse the mapping from altered names back to the original column names
        reverse_mapping = {v: k for k, v in column_mapping.items()}
        labels = [reverse_mapping.get(name, name) for name in mi_matrix.columns if name in reverse_mapping]
        plt.xticks(ticks=np.arange(len(labels)) + 0.5, labels=labels, rotation=45, ha='right')
        plt.title('Mutual Information Heatmap')
        plt.show()
    
    @staticmethod
    def generate_continuous_statistics(df, column_names, output_file):
        # Replace negative values (-1, -3) with NaN
        df[column_names] = df[column_names].replace([-1, -3], np.nan)
        
        missing_columns = [column for column in column_names if column not in df.columns]
        if missing_columns:
            raise ValueError(f"Columns not found in the dataset: {', '.join(missing_columns)}")
        
        with open(output_file, 'w') as file:
            for column in column_names:
                data = df[column].dropna()  # Drop NaN values
                statistics = {
                    'N': data.count(),
                    'Mean': round(data.mean(), 2),
                    'Standard Deviation': round(data.std(), 2)
                }
                
                file.write(f"\nDescriptive Statistics for '{column}':\n")
                for stat, value in statistics.items():
                    file.write(f"{stat}: {value}\n")

    @staticmethod
    def generate_categorical_statistics(df, column_names, output_file):
        missing_columns = [column for column in column_names if column not in df.columns]
        if missing_columns:
            raise ValueError(f"Columns not found in the dataset: {', '.join(missing_columns)}")
        
        with open(output_file, 'w') as file:
            for column in column_names:
                data = df[column]
                counts = data.value_counts()
                total = data.count()
                frequencies = counts / total
                
                file.write(f"\nFrequency Distribution for '{column}':\n")
                file.write(f"Total: {total}\n")
                for category, count in counts.items():
                    frequency = frequencies[category]
                    file.write(f"{category}: Count = {count}, Frequency = {frequency:.2f}\n")
    
    @staticmethod
    def plot_continuous_distribution(df, column_names, column_labels, output_folder):
        df[column_names] = df[column_names].replace([-1, -3], np.nan)
        
        missing_columns = [column for column in column_names if column not in df.columns]
        if missing_columns:
            raise ValueError(f"Columns not found in the dataset: {', '.join(missing_columns)}")
        
        num_cols = 3
        num_rows = (len(column_names) + num_cols - 1) // num_cols
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 5))
        
        for ax, column, label in zip(axes.flatten(), column_names, column_labels):
            sns.kdeplot(df[column].dropna(), ax=ax, bw_adjust=0.5)
            ax.set_xlabel(label)
            ax.set_ylabel('Density')
        
        for i in range(len(column_names), num_rows * num_cols):
            fig.delaxes(axes.flatten()[i])
        
        plt.tight_layout()
        output_file = f"{output_folder}/density_plots.png"
        plt.savefig(output_file)
        plt.close()
        