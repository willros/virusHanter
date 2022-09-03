import pandas as pd
import altair as alt

alt.data_transformers.disable_max_rows()


def bar_chart_kaiju_raw(
    file: str, cutoff: float = 0.01, number: int = 10
) -> alt.vegalite.v4.api.Chart:
    """
    Returns a bar chart of the taxnomies from the kaiju raw csv in the cleaned_files folder

    :param str file: Path to the kaiju raw csv in the cleaned_files folder
    :param float cutoff: Cutoff of percent the taxonomy level is present in. Default = 0.01
    :param int number: The number bars to plot. Default = 10
    :return: altair.vegalite.v4.api.Chart
    """

    df = (
        pd.read_csv(file)
        .groupby(["taxon_id", "percent", "taxon_name"], as_index=False)
        .agg(taxonomy=("taxonomy", list))
        .sort_values("percent", ascending=False)
        .loc[lambda x: x.percent > cutoff]
    )

    return (
        alt.Chart(df.head(number))
        .mark_bar()
        .encode(
            alt.X("percent:Q", axis=alt.Axis(format=".1%"), title="Percent of reads"),
            alt.Y("taxon_name:N", sort="-x", title=None),
            alt.Color("taxon_name:N", title=None),
            tooltip=[alt.Tooltip("taxonomy:O"), alt.Tooltip("percent:Q", format=".1%")],
        )
        .properties(width=500, height=500, title="Kaiju classification raw")
    )
