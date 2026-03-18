import seaborn as sns
import matplotlib.pyplot as plt

def ts_plot_delta(plot_df_top, plot_df_bottom, spcname, outdir):

    plt.figure(figsize=(12, 6))

    """
    Create a 2-panel time series plot:
    - Top: Baseline vs IMPROVE observations
    - Bottom: Experiment differences from baseline
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # Top subplot: Baseline vs Observations
    sns.lineplot(
        data=plot_df_top,
        x="time",
        y="PM25",
        hue="Source",
        estimator="median",
        errorbar=("pi", 50),
        alpha=0.8,
        ax=axes[0]
    )
    axes[0].set_title(f"Surface {spcname} PM2.5: Baseline vs. IMPROVE (Median & IQR)")
    axes[0].set_ylabel(r"Concentration ($\mu g/m^3$)")
    axes[0].set_xlabel("")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(title="Source")

    # Bottom subplot: Experiment Differences
    sns.lineplot(
        data=plot_df_bottom,
        x="time",
        y="PM25",
        hue="Source",
        estimator="median",
        errorbar=("pi", 50),
        alpha=0.8,
        ax=axes[1]
    )
    axes[1].axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    axes[1].set_title(f"Experiment Differences from Baseline")
    axes[1].set_ylabel(r"Difference ($\mu g/m^3$)")
    axes[1].set_xlabel("Date")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(title="Experiment")

    plt.tight_layout()
    outf = f'{outdir}/{spcname}_delta.png'
    plt.savefig(outf, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
