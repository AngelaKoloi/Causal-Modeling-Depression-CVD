## Installation
Clone the repository

Go to the project directory:
`cd Causal_Discovery`

Install the required Python packages:
`pip install -r requirements.txt`

## Usage
### Data Preparation
Place the raw datasets under the `data` directory. Ensure that the datasets are properly named and formatted as expected by the scripts.

### Exploratory Data Analysis (EDA)
To generate plots and statistics for exploratory data analysis, open and run all cells in the `eda.ipynb` notebook. This will provide insights into the raw data.

### Feature Elimination
To eliminate the metabolites from the raw dataset:

1. Open and run all cells in the `eliminate_metabolites.ipynb` notebook.

2. This notebook can only process one dataset at a time. If you have multiple datasets (e.g., from different cohorts), you need to:

- Update the file path to point to the next dataset.

- Rename the trained model for the second dataset to avoid overwriting.

3. Running this notebook will generate two trained XGBoost classifier models: one for depression and one for metabolic syndrome.

  
### Shared Metabolites Analysis
To analyze metabolites found in the trained models and find the overlap with metabolites suggested by expert opinions. Follow the instructions provided in the respective notebook (if available).

### Prepare Variables for Causal Graph Construction
1. Run all cells in the `prepare_wave1.ipynb` notebook to process and clean the baseline dataset.

2. Run all cells in the `prepare_wave2.ipynb` notebook to process and clean the follow up dataset.

These steps will produce cleaned and discretized datasets for each wave.

### Generate the Network
1. Navigate to the `src` directory:

```
cd src
```
2. Run the `causal_discovery.R` script to generate the network. This will produce a .RData file, which will be used in network.R. 

### Visualize the Network
Open demo.html in your browser