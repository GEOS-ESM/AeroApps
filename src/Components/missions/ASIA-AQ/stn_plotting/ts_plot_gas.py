import seaborn as sns
import matplotlib.pyplot as plt

def ts_plot(plot_df,spcname,units,outdir):

    plt.figure(figsize=(12, 6))

    # Plot both lines with 25th-75th percentile shading
    # ("pi", 50) creates an interval covering 50% of the data (25th to 75th)
    sns.lineplot(
        data=plot_df,
        x="time",
        y=spcname,
        hue="Source",
        estimator="median",
        errorbar=("pi", 50),
        alpha=0.8
    )

    plt.title(f"Comparison of Surface {spcname}: GEOS (Median & 25-75th Percentile)")
    plt.ylabel(f"Concentration ({units})")
    plt.xlabel("Date")
    plt.grid(True, alpha=0.3)
    outf = f'{outdir}/{spcname}.png'
    plt.savefig(outf)
    plt.show()
