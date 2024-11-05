import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# read in filtered data
input_data = pd.read_csv('filtered_data.csv')

# pivot table to count occurrences of each cluster (event) in each group
try:    
    pivot_table = input_data.pivot_table(index='cluster', columns='genotype', aggfunc='size', fill_value=0)
    print('Successfully pivoted the filtered dataset')

except Exception as e:
    print(f"Error making heatmap: {e}")
    raise

# make heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(pivot_table, annot=True, linewidth=.5, cbar_kws={'label': 'Occurrences'})
#sns.heatmap(pivot_table, annot=True, cmap="YlGnBu", cbar_kws={'label': 'Occurrences'})
plt.title("Cluster Occurrences in Rescue Induced Groups")
plt.xlabel("Groups")
plt.ylabel("Clusters")

# save heatmap and clear memory
plt.savefig("rescued_heatmap.png", format="png", dpi=300)
plt.clf()

print("Successfully created and saved heatmap")