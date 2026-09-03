#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
import numpy as np

# Read the cleaned TSV data
data = pd.read_csv("doc/tables/ethnic_background_fixed.tsv", sep='\t')

# Add 1000G superpopulation assignments based on country/ethnicity
def assign_superpop(country):
    if country in ['Germany', 'Austria', 'France', 'Spain', 'Croatia', 'Poland', 'Russia']:
        return 'EUR'
    elif country in ['Taiwan', 'Vietnam']:
        return 'EAS'
    elif country in ['India']:
        return 'SAS'
    elif country in ['Mexico', 'Chile']:
        return 'AMR'
    else:
        return 'EUR'  # Default to European

data['SUPERPOP'] = data['Country'].apply(assign_superpop)

# Add jitter to coordinates to show overlapping points
np.random.seed(42)
data['Longitude_jitter'] = data['Longitude'] + np.random.normal(0, 0.5, len(data))
data['Latitude_jitter'] = data['Latitude'] + np.random.normal(0, 0.5, len(data))

# Color palette for superpopulations (same as PCA plot)
superpop_colors = {
    'EUR': '#1976D2',  # Blue
    'AFR': '#FF6F00',  # Orange
    'EAS': '#388E3C',  # Green
    'SAS': '#D32F2F',  # Red
    'AMR': '#7B1FA2'   # Purple
}

# Create the figure with Robinson projection (oval world map)
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.Robinson())

# Add map features with coarser resolution
ax.add_feature(cfeature.COASTLINE, linewidth=0.8, color='black')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, color='black', alpha=0.7)
ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.1)
ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)

# Set global extent
ax.set_global()

# Plot sample points by superpopulation
for superpop in data['SUPERPOP'].unique():
    subset = data[data['SUPERPOP'] == superpop]
    ax.scatter(subset['Longitude_jitter'], subset['Latitude_jitter'],
              c=superpop_colors[superpop], s=100, alpha=0.9,
              edgecolors='black', linewidth=1.0, 
              label=superpop, transform=ccrs.PlateCarree(),
              zorder=5)

# Customize the plot
plt.title('T2T Sample Geographic Distribution', fontsize=18, fontweight='bold', pad=30)

# Add legend
legend = plt.legend(title='1000G Superpopulation', 
                   title_fontsize=14, fontsize=12,
                   loc='lower center', bbox_to_anchor=(0.5, -0.1),
                   ncol=5, frameon=True, fancybox=True, shadow=True)
legend.get_title().set_fontweight('bold')

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('doc/figures/ethnicity_world_map.png', dpi=300, bbox_inches='tight', 
           facecolor='white', edgecolor='none')
plt.savefig('doc/figures/ethnicity_world_map.pdf', bbox_inches='tight', 
           facecolor='white', edgecolor='none')

# Print summary
print("Map saved as 'doc/figures/ethnicity_world_map.png' and '.pdf'")
print("Sample distribution by superpopulation:")
print(data['SUPERPOP'].value_counts())
print("\nSample distribution by country:")
print(data['Country'].value_counts())

plt.show()