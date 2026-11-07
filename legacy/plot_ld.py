import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Indlæs simuleringsdata
try:
    df_ld = pd.read_csv('ld_decay.csv')
    df_map = pd.read_csv('calculated_linkage_map.csv')
except FileNotFoundError as e:
    print(f"Fejl: Kunne ikke finde datafilerne. Kør C++ simulatoren først.\nDetaljer: {e}")
    exit()

# 2. Opsæt side-om-side layout (1 række, 2 kolonner)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
chromosomes = df_ld['Chromosome'].unique()

# =========================================================================
# GRAPH 1: LINKAGE DISEQUILIBRIUM (R²) DECAY (VENSTRE SIDE)
# =========================================================================
for idx, chr_id in enumerate(chromosomes):
    chr_ld = df_ld[df_ld['Chromosome'] == chr_id].sort_values('Pairwise_Distance_bp')
    
    if chr_ld.empty:
        continue
        
    x_ld = chr_ld['Pairwise_Distance_bp']
    y_ld = chr_ld['R2_Linkage_Disequilibrium']
    
    # Tegn rå scatter punkter for markør-par
    ax1.scatter(x_ld, y_ld, alpha=0.12, s=6, color=colors[idx % len(colors)],
                label=f'Markør-par (Kromosom {chr_id})')
    
    # Beregn og tegn rullende gennemsnit (Trendlinje)
    if len(chr_ld) > 10:
        window_size = max(5, len(chr_ld) // 20)
        rolling_x = x_ld.rolling(window=window_size, center=True).mean()
        rolling_y = y_ld.rolling(window=window_size, center=True).mean()
        ax1.plot(rolling_x, rolling_y, color=colors[idx % len(colors)], 
                 linewidth=2.5, label=f'Gns. LD Trend (Kromosom {chr_id})')

ax1.set_title('Linkage Disequilibrium ($R^2$) Decay Profile', fontsize=12, fontweight='bold')
ax1.set_xlabel('Fysisk afstand mellem markører (bp)', fontsize=10)
ax1.set_ylabel('LD Koefficient ($R^2$)', fontsize=10)
ax1.set_ylim(-0.02, 1.05)
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=9)

# =========================================================================
# GRAPH 2: CALCULATED GENETIC LINKAGE MAP (HØJRE SIDE)
# =========================================================================
for idx, chr_id in enumerate(chromosomes):
    chr_map = df_map[df_map['Chromosome'] == chr_id].sort_values('Physical_Position_bp')
    
    if chr_map.empty:
        continue
        
    x_map = chr_map['Physical_Position_bp']
    y_map = chr_map['Estimated_Genetic_Position_cM']
    
    # Tegn det estimerede genetiske kort
    ax2.plot(x_map, y_map, color=colors[idx % len(colors)], linewidth=3,
             label=f'Genetisk kort (Kromosom {chr_id})')

ax2.set_title('Calculated Genetic Linkage Map (cM vs bp)', fontsize=12, fontweight='bold')
ax2.set_xlabel('Fysisk position på kromosom (bp)', fontsize=10)
ax2.set_ylabel('Genetisk position (centiMorgan, cM)', fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='upper left', frameon=True, facecolor='white', edgecolor='none', fontsize=9)

# Tilføj en overordnet signatur og gem figuren
plt.suptitle('Side-by-Side Recombination Lookup & Validation', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('ld_and_linkage_map_validation.png')
plt.show()

print("📈 Side-om-side visualisering genereret med succes!")
print(" ↳ Gemt som filen: 'ld_and_linkage_map_validation.png'")
