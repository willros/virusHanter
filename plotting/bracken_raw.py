import pandas as pd
import altair as alt

alt.data_transformers.disable_max_rows()


def df_bracken_species_raw(
    file: str, level: str = "species", cutoff: float = 0.001, virus_only: bool = True
) -> pd.DataFrame:
    """
    Returns a df for the taxonomy found in the cleaned raw bracken report.
    Used to generate bar plots of the different taxonomies.

    :param str file: Path to the bracken report.
    :param str level: Level of taxonomy. [domain, phylum, class, order, family, genus, species]. Default = 'species'
    :param float cutoff: Cutoff of percent the taxonomy level is present in. Default = 0.05
    :param bool virus_only: Only include Viruses. Default = True
    :return: pd.DataFrame
    """

    taxonomy = {
        "domain": "D",
        "phylum": "P",
        "class": "K",
        "order": "O",
        "family": "F",
        "genus": "G",
        "species": "S",
    }

    df = (
        pd.read_csv(file)
        .loc[lambda x: x.level == taxonomy[level]]
        .loc[lambda x: x.percent > cutoff]
        .sort_values("percent", ascending=False)
    )

    if virus_only:
        return df.loc[lambda x: x.domain == "Virus"]
    return df


def bar_chart_bracken_raw(
    file: str,
    level: str = "species",
    cutoff: float = 0.001,
    number: int = 10,
    virus_only=True,
) -> alt.vegalite.v4.api.Chart:
    """
    Returns a bar chart of the taxnomies from the bracken species file in the cleaned_files folder.

    :param str file: Path to the cleaned bracken report in the cleaned_files folder.
    :param str level: Level of taxonomy. [domain, phylum, class, order, family, genus, species]. Default = 'species'
    :param float cutoff: Cutoff of percent the taxonomy level is present in. Default = 0.05
    :param int number: The number bars to plot. Default = 10
    :return: altair.vegalite.v4.api.Chart
    """

    df = df_bracken_species_raw(file, level, cutoff, virus_only).head(number)

    return (
        alt.Chart(df, title="Bracken classification raw")
        .mark_bar()
        .encode(
            alt.X("percent:Q", axis=alt.Axis(format=".1%"), title="Percent of reads"),
            alt.Y("name:N", sort="-x", title=None),
            alt.Color("name:N", title=None),
            tooltip=[alt.Tooltip("domain:N"), alt.Tooltip("percent:Q", format=".1%")],
        )
        .properties(width=500, height=500)
    )


def pie_chart_bracken_raw(file: str) -> alt.vegalite.v4.api.Chart:
    """
    Returns a pie chart of the kingdoms from the bracken species file in the cleaned_files folder.

    :param str file: Path to the cleaned bracken report in the cleaned_files folder.
    :return: altair.vegalite.v4.api.Chart
    """

    df = df_bracken_species_raw(file, virus_only=False, level="domain")

    pie = (
        alt.Chart(df, title="Abundance of Kingdoms")
        .mark_arc(innerRadius=50)
        .encode(
            theta=alt.Theta(field="percent", type="quantitative"),
            color=alt.Color(field="name", type="nominal", legend=None),
        )
        .properties(width=500)
    )

    text = pie.mark_text(radius=174, size=10).encode(text="name:N")

    return pie + text
