import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from kcat.txt and mw.txt
with open('kcat.txt', 'r') as file:
    kcat_data = np.array([float(line.strip()) for line in file.readlines() if float(line.strip()) <= 100])

with open('mw.txt', 'r') as file:
    mw_data = np.array([float(line.strip()) for line in file.readlines() if float(line.strip()) <= 100])

# Create a DataFrame for the data
df_kcat = pd.DataFrame({'Kcat': kcat_data})
df_mw = pd.DataFrame({'MW': mw_data})

# Plotting function with formatting
def plot_density(data, variable, color, title, xlabel, ylabel):
    plt.figure(figsize=(8, 5))
    sns.set_style("white")
    sns.kdeplot(data=data, x=variable, color=color, fill=True, linewidth=2)
    plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.xticks(fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.title(title)
    plt.tight_layout()
    plt.savefig('Kcat_model.png')
    plt.show()

# Plot original ridge plot for Kcat
plot_density(df_kcat, 'Kcat', 'skyblue', '', 'Kcat (1/s)', 'Density')

# Plot original ridge plot for MW
plot_density(df_mw, 'MW', 'salmon', '', 'MW (kDa)', 'Density')
# Perform Monte Carlo simulation for Kcat
num_samples_kcat = 373
num_simulations_kcat = 100
monte_carlo_results_kcat = []

for _ in range(num_simulations_kcat):
    kcat_samples = np.random.choice(kcat_data, num_samples_kcat, replace=True)
    monte_carlo_results_kcat.append(kcat_samples)

# Perform Monte Carlo simulation for MW
num_samples_mw = 12
num_simulations_mw = 100
monte_carlo_results_mw = []

for _ in range(num_simulations_mw):
    mw_samples = np.random.choice(mw_data, num_samples_mw, replace=True)
    monte_carlo_results_mw.append(mw_samples)

# Write Kcat Monte Carlo simulation results to a text file
with open('Kcat_MC.txt', 'w') as file:
    for i, result in enumerate(monte_carlo_results_kcat):
        file.write(f"Simulation {i+1}:\n")
        for kcat in result:
            file.write(f"{kcat}\n")
        file.write("\n")

# Write MW Monte Carlo simulation results to a text file
with open('MW_MC.txt', 'w') as file:
    for i, result in enumerate(monte_carlo_results_mw):
        file.write(f"Simulation {i+1}:\n")
        for mw in result:
            file.write(f"{mw}\n")
        file.write("\n")

# Plot the similarity of calculated Kcat distribution with original Kcat distribution
plt.figure(figsize=(8, 5))
sns.set_style("white")
sns.kdeplot(data=kcat_data, color='blue', fill=True, linewidth=2, label='Original Kcat Distribution')
for i in range(num_simulations_kcat):
    sns.kdeplot(data=monte_carlo_results_kcat[i], color='orange', alpha=0.3, linewidth=1)
plt.xlabel('Kcat (1/s)', fontsize=12, fontweight='bold')
plt.ylabel('Density', fontsize=12, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.legend(['Original Kcat Distribution', 'Monte Carlo Simulation'], loc='upper right')
plt.tight_layout()
plt.savefig('Kcat_distribution.png')
plt.show()

# Similar approach for the MW distribution
plt.figure(figsize=(8, 5))
sns.set_style("white")
sns.kdeplot(data=mw_data, color='blue', fill=True, linewidth=2, label='Original MW Distribution')
for i in range(num_simulations_mw):
    sns.kdeplot(data=monte_carlo_results_mw[i], color='orange', alpha=0.3, linewidth=1)
plt.xlabel('MW (kDa)', fontsize=12, fontweight='bold')
plt.ylabel('Density', fontsize=12, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.legend(['Original MW Distribution', 'Monte Carlo Simulation'], loc='upper right')
plt.tight_layout()
plt.savefig('MW_distribution.png')
plt.show()