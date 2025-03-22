import numpy as np
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from geneticalgorithm import geneticalgorithm as ga
import pandas as pd

# Load and clean the data
df = pd.read_csv("C:/Users/Rober/Downloads/parameter_rmse_data.csv")
df_clean = df[df['RMSE_Total_Young'] < 100].copy()

# Define features and target
features = df_clean.drop(columns=['KTotal', 'small_wr', 'RMSE_Total_Young'])
target = df_clean['RMSE_Total_Young']

# Train the XGBoost model
X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
model = XGBRegressor(objective='reg:squarederror', random_state=42)
model.fit(X_train, y_train)

# Define the GA fitness function
def fitness_function(x):
    # Model expects 2D array
    x_reshaped = np.array(x).reshape(1, -1)
    prediction = model.predict(x_reshaped)[0]
    return prediction  # Since we want to minimize RMSE_Total_Young

# Define variable bounds (based on data min/max or logical limits)
var_bounds = []
for col in features.columns:
    col_min = features[col].min()
    col_max = features[col].max()
    var_bounds.append((col_min, col_max))

# GA parameters
model_g = ga(
    function=fitness_function,
    dimension=len(features.columns),
    variable_type='real',
    variable_boundaries=np.array(var_bounds),
    algorithm_parameters={
        'max_num_iteration': 200,
        'population_size': 50,
        'mutation_probability': 0.1,
        'elit_ratio': 0.01,
        'crossover_probability': 0.5,
        'parents_portion': 0.3,
        'crossover_type': 'uniform',
        'max_iteration_without_improv': 20
    }
)

# Run the GA
model_g.run()

# Show the best result
best_inputs = model_g.output_dict['variable']
best_prediction = model_g.output_dict['function']
print("\nBest Inputs:", dict(zip(features.columns, best_inputs)))
print("Predicted RMSE_Total_Young:", best_prediction)
