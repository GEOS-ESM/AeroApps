import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def kde_plot(df_baseline, df_imp, df_oxh, df_oxm, spcname, outdir):
    """
    Create a 1x3 2D KDE (density) plot:
    - Left: Baseline vs IMPROVE observations
    - Middle: OxH vs IMPROVE observations
    - Right: OxM vs IMPROVE observations
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Common plot settings
    plot_configs = [
        (df_baseline, 'Baseline', axes[0]),
        (df_oxh, 'OxH', axes[1]),
        (df_oxm, 'OxM', axes[2])
    ]
    
    for df_model, label, ax in plot_configs:
        # Merge model and obs data
        merged = df_model[['time', 'station', 'PM25']].merge(
            df_imp[['time', 'station', 'PM25']],
            on=['time', 'station'],
            suffixes=('_model', '_obs')
        ).dropna()
        
        if len(merged) < 10:  # Need enough points for KDE
            ax.text(0.5, 0.5, 'Insufficient data for KDE', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"{label} vs IMPROVE")
            continue
        
        # Create 2D KDE plot
        sns.kdeplot(
            data=merged,
            x='PM25_obs',
            y='PM25_model',
            cmap='viridis',
            fill=True,
            levels=10,
            thresh=0.05,
            ax=ax
        )
        
        
        # Calculate statistics
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            merged['PM25_obs'], merged['PM25_model']
        )
        rmse = np.sqrt(((merged['PM25_model'] - merged['PM25_obs'])**2).mean())
        bias = (merged['PM25_model'] - merged['PM25_obs']).mean()
        n = len(merged)
        
        # Add 1:1 line
#        max_val = max(merged['PM25_obs'].max(), merged['PM25_model'].max())
#        min_val = min(merged['PM25_obs'].min(), merged['PM25_model'].min())
        max_val = 10
        min_val = 0
        ax.plot([min_val, max_val], [min_val, max_val], 
               'k--', linewidth=2, alpha=0.8, label='1:1 line')
        
        # Add regression line
        x_line = np.array([min_val, max_val])
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'r-', linewidth=2, 
               label=f'y={slope:.2f}x+{intercept:.2f}')
        
        # Add statistics text
        stats_text = f'R² = {r_value**2:.3f}\n'
        stats_text += f'RMSE = {rmse:.3f}\n'
        stats_text += f'Bias = {bias:.3f}\n'
        stats_text += f'N = {n}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='white', alpha=0.8, edgecolor='black'), 
               fontsize=9)
        
        # Labels and formatting
        ax.set_xlabel(r"IMPROVE ($\mu g/m^3$)")
        ax.set_ylabel(r"Model ($\mu g/m^3$)")
        ax.set_title(f"{label} vs IMPROVE")
        ax.legend(loc='lower right', fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Set equal aspect ratio
        ax.set_aspect('equal', adjustable='box')
    
    plt.suptitle(f"Surface {spcname} PM2.5 Density Plots", fontsize=14, y=1.02)
    plt.tight_layout()
    outf = f'{outdir}/{spcname}_kde.png'
    plt.savefig(outf, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
