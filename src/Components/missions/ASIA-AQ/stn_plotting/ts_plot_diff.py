import seaborn as sns
import matplotlib.pyplot as plt

def ts_plot_diff(plot_df, spcname, outdir):
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Top subplot: GEOS vs IMPROVE time series
    sns.lineplot(
        data=plot_df,
        x="time",
        y="PM25",
        hue="Source",
        estimator="median",
        errorbar=("pi", 50),
        alpha=0.8,
        ax=axes[0]
    )
    axes[0].set_title(f"Surface {spcname} PM2.5: GEOS vs. IMPROVE (Median & IQR)")
    axes[0].set_ylabel(r"Concentration ($\mu g/m^3$)")
    axes[0].set_xlabel("")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(title="Source")
    
    # Bottom subplot: Difference (GEOS - IMPROVE)
    # First, pivot the data to calculate differences
    df_geos = plot_df[plot_df['Source'] == 'GEOS'].copy()
    df_imp = plot_df[plot_df['Source'] == 'IMPROVE'].copy()
    
    # Merge to align time and station
    df_diff = df_geos[['time', 'station', 'PM25']].merge(
        df_imp[['time', 'station', 'PM25']],
        on=['time', 'station'],
        suffixes=('_geos', '_imp')
    )
    df_diff['PM25'] = df_diff['PM25_geos'] - df_diff['PM25_imp']
    df_diff['Source'] = 'GEOS - IMPROVE'
    
    sns.lineplot(
        data=df_diff,
        x="time",
        y="PM25",
        estimator="median",
        errorbar=("pi", 50),
        alpha=0.8,
        color='purple',
        ax=axes[1]
    )
    
    # Add zero reference line
    axes[1].axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    axes[1].set_title(f"Difference (GEOS - IMPROVE)")
    axes[1].set_ylabel(r"Difference ($\mu g/m^3$)")
    axes[1].set_xlabel("Date")
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    outf = f'{outdir}/{spcname}.png'
    plt.savefig(outf, dpi=150, bbox_inches='tight')
#    plt.show()
    plt.close()
