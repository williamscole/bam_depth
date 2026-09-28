#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
from scipy import stats

# Check arguments
if len(sys.argv) != 2:
    print("Error: BUILD argument is required")
    print("Usage: python analyze_depth.py [b37|b38]")
    sys.exit(1)

BUILD = sys.argv[1]

if BUILD not in ["b37", "b38"]:
    print("Error: BUILD must be either 'b37' or 'b38'")
    print("Usage: python analyze_depth.py [b37|b38]")
    sys.exit(1)

# Define paths
PROJ_DIR = Path("/zivlab/data3/colew/Proj_GLIMPSE")
outputdir = PROJ_DIR / "_Output" / "depth" / BUILD / "mosdepth"
bamfile_list = PROJ_DIR / "_Data" / f"{BUILD}_geri_bamfiles.txt"
summary_file = outputdir / "depth_summary.tsv"
concorddir = PROJ_DIR / "_Output" / BUILD / "glimpse_concord"

# Load summary data
df_summary = pd.read_csv(summary_file, sep='\t')

# Load BAM file list and extract sample_id and batch
bam_data = []
suffix = f"_{BUILD}.DeDup.bam"

with open(bamfile_list, 'r') as f:
    for line in f:
        bam_path = line.strip()
        # Extract sample_id
        sample_id = Path(bam_path).name.replace(suffix, "")
        
        # Extract batch (e.g., "apr2024" from /zivlab/data3/geri_wes/apr2024/)
        parts = bam_path.split('/')
        batch = None
        for i, part in enumerate(parts):
            if part == "geri_wes" and i + 1 < len(parts):
                batch = parts[i + 1]
                break
        
        bam_data.append({'sample_id': sample_id, 'batch': batch})

df_bam = pd.DataFrame(bam_data)

# Merge with summary
df = df_summary.merge(df_bam, on='sample_id', how='left')

# Define month order for chronological sorting
month_order = ['jan2024', 'feb2024', 'mar2024', 'apr2024', 'may2024', 'jun2024',
               'july2024', 'aug2024', 'sep2024', 'oct2024', 'nov2024', 'dec2024',
               'jan2025', 'feb2025', 'mar2025', 'apr2025', 'may2025', 'jun2025',
               'july2025', 'aug2025', 'sep2025', 'oct2025', 'nov2025', 'dec2025']

# Filter to only batches present in data and maintain chronological order
present_batches = [b for b in month_order if b in df['batch'].values]
df['batch'] = pd.Categorical(df['batch'], categories=present_batches, ordered=True)

# Sort by batch for plotting
df = df.sort_values('batch')

print(f"Loaded {len(df)} samples")
print(f"Batches: {df['batch'].value_counts().sort_index()}")

# ============================================
# Load concordance data and compute R²
# ============================================

def compute_r2_from_sums(sum_x, sum_y, sum_x2, sum_y2, sum_xy, n):
    """Compute R² from sum statistics"""
    if n == 0:
        return np.nan
    
    # Pearson correlation coefficient
    numerator = n * sum_xy - sum_x * sum_y
    denominator = np.sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
    
    if denominator == 0:
        return np.nan
    
    r = numerator / denominator
    r2 = r ** 2
    
    return r2

concordance_data = []

for sample_id in df['sample_id']:
    concord_file = concorddir / f"{sample_id}_corr.txt"
    
    if not concord_file.exists():
        continue
    
    # Load concordance file
    df_concord = pd.read_csv(concord_file, sep='\t')
    
    # Process each region type
    for region_type in ['all', 'exome', 'nonexome']:
        # Subset by label prefix
        subset = df_concord[df_concord['label'].str.startswith(f"{region_type}/")]
        
        if len(subset) == 0:
            continue
        
        # Aggregate sum statistics
        total_n = subset['n'].sum()
        total_sum_x = subset['sum_x'].sum()
        total_sum_y = subset['sum_y'].sum()
        total_sum_x2 = subset['sum_x2'].sum()
        total_sum_y2 = subset['sum_y2'].sum()
        total_sum_xy = subset['sum_xy'].sum()
        
        # Compute R²
        r2 = compute_r2_from_sums(total_sum_x, total_sum_y, total_sum_x2, 
                                   total_sum_y2, total_sum_xy, total_n)
        
        concordance_data.append({
            'sample_id': sample_id,
            'region_type': region_type,
            'r2': r2
        })

df_concordance = pd.DataFrame(concordance_data)

# Pivot to wide format
df_concordance_wide = df_concordance.pivot(index='sample_id', 
                                           columns='region_type', 
                                           values='r2').reset_index()
df_concordance_wide.columns.name = None
df_concordance_wide = df_concordance_wide.rename(columns={
    'all': 'r2_all',
    'exome': 'r2_exome',
    'nonexome': 'r2_nonexome'
})

# Merge concordance with main dataframe
df = df.merge(df_concordance_wide, on='sample_id', how='left')

print(f"Samples with concordance data: {df_concordance_wide.shape[0]}")

# ============================================
# Plot 1: Boxplot of on-target depth by batch
# ============================================
fig, ax = plt.subplots(figsize=(12, 6))

# Create boxplot
bp = sns.boxplot(data=df, x='batch', y='mean_on', ax=ax, hue='batch')

# Calculate means and counts for each batch and add annotations
batch_stats = df.groupby('batch')['mean_on'].agg(['mean', 'count'])
for i, batch in enumerate(present_batches):
    if batch in batch_stats.index:
        mean_val = batch_stats.loc[batch, 'mean']
        count_val = int(batch_stats.loc[batch, 'count'])
        # Position text above the box
        y_pos = df[df['batch'] == batch]['mean_on'].max()
        ax.text(i, y_pos, f'mean={mean_val:.1f}\nn={count_val}', 
                ha='center', va='bottom', fontweight='bold', fontsize=9)

ax.set_xlabel('Batch', fontsize=12)
ax.set_ylabel('Mean On-Target Depth', fontsize=12)
ax.set_title(f'On-Target Depth Distribution by Batch ({BUILD})', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save plot
plot1_path = outputdir / "mosdepth_boxplot.png"
plt.savefig(plot1_path, dpi=500, bbox_inches='tight')
print(f"Saved boxplot to {plot1_path}")
plt.close()

# ============================================
# Plot 2: Scatter plot of on vs off target depth
# ============================================
fig, ax = plt.subplots(figsize=(10, 8))

# Create scatter plot colored by batch
for batch in present_batches:
    batch_data = df[df['batch'] == batch]
    ax.scatter(batch_data['mean_on'], batch_data['mean_off'], 
               label=batch, alpha=0.6, s=50)

# Add regression line
slope, intercept, r_value, p_value, std_err = stats.linregress(df['mean_on'], df['mean_off'])
x_range = np.array([df['mean_on'].min(), df['mean_on'].max()])
y_pred = slope * x_range + intercept
ax.plot(x_range, y_pred, 'k--', linewidth=2, 
        label=f'Regression (R²={r_value**2:.3f})')

ax.set_xlabel('Mean On-Target Depth', fontsize=12)
ax.set_ylabel('Mean Off-Target Depth', fontsize=12)
ax.set_title(f'On-Target vs Off-Target Depth ({BUILD})', fontsize=14)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()

# Save plot
plot2_path = outputdir / "mosdepth_scatter.png"
plt.savefig(plot2_path, dpi=500, bbox_inches='tight')
print(f"Saved scatter plot to {plot2_path}")
plt.close()

# ============================================
# Plot 3: R² vs Depth (3 panels)
# ============================================

# Filter to only samples with concordance data
df_with_concord = df.dropna(subset=['r2_all', 'r2_exome', 'r2_nonexome'])

print(f"Plotting R² vs depth for {len(df_with_concord)} samples with concordance data")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel 1: All regions vs genome depth
ax = axes[0]
for batch in present_batches:
    batch_data = df_with_concord[df_with_concord['batch'] == batch]
    ax.scatter(batch_data['mean_genome'], batch_data['r2_all'], 
               label=batch, alpha=0.6, s=50)
ax.set_xlabel('Mean Genome Depth', fontsize=12)
ax.set_ylabel('R² (All)', fontsize=12)
ax.set_title('All Sites', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()

# Panel 2: Exome regions vs on-target depth
ax = axes[1]
for batch in present_batches:
    batch_data = df_with_concord[df_with_concord['batch'] == batch]
    ax.scatter(batch_data['mean_on'], batch_data['r2_exome'], 
               label=batch, alpha=0.6, s=50)
ax.set_xlabel('Mean On-Target Depth', fontsize=12)
ax.set_ylabel('R² (Exome)', fontsize=12)
ax.set_title('Exome Sites', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()

# Panel 3: Non-exome regions vs off-target depth
ax = axes[2]
for batch in present_batches:
    batch_data = df_with_concord[df_with_concord['batch'] == batch]
    ax.scatter(batch_data['mean_off'], batch_data['r2_nonexome'], 
               label=batch, alpha=0.6, s=50)
ax.set_xlabel('Mean Off-Target Depth', fontsize=12)
ax.set_ylabel('R² (Non-Exome)', fontsize=12)
ax.set_title('Non-Exome Sites', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()

# Save plot
plot3_path = outputdir / "mosdepth_r2_by_depth.png"
plt.savefig(plot3_path, dpi=500, bbox_inches='tight')
print(f"Saved R² vs depth plot to {plot3_path}")
plt.close()

print("Analysis complete!")