import pandas as pd
import matplotlib.pyplot as plt

# Create a sample dataframe
data = {'A': [1, 2, 3, 4, 5],
        'B': [10, 20, 30, 40, 50],
        'C': [100, 200, 300, 400, 500]}
df = pd.DataFrame(data)

# Use the describe function to get the statistics
df_stats = df.describe()

# Visualize the statistics using plots and save the figures
# Box plots for showing min, 25th percentile, median, 75th percentile, and max
box_plot = df_stats.plot(kind='box', vert=False)
plt.title('Box Plot of Descriptive Statistics')
plt.xlabel('Value')
plt.ylabel('Statistics')
plt.savefig('box_plot.png', bbox_inches='tight')
plt.close()

# Bar plots for showing count, mean, standard deviation, and quartile ranges
bar_plot = df_stats.loc[['count', 'mean', 'std', '25%', '50%', '75%']].plot(kind='bar', legend=False)
plt.title('Bar Plot of Descriptive Statistics')
plt.xlabel('Statistics')
plt.ylabel('Value')
plt.savefig('bar_plot.png', bbox_inches='tight')
plt.close()

# Histograms for showing distribution
histograms = df.hist()
plt.suptitle('Histograms of Descriptive Statistics')
plt.tight_layout()
plt.savefig('histograms.png', bbox_inches='tight')
plt.close()