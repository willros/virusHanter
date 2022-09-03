import pandas as pd
import altair as alt
import numpy as np

alt.data_transformers.disable_max_rows()


def bar_chart_kaiju_megahit(file: str) -> alt.vegalite.v4.api.Chart:
    """
    Plots taxonomy abundancy of the contigs from the megahit assembly.
    Needs the "kaiju_out" file.
    :param str file: Path to the csv file for the kaiju out file on the contigs.
    :return: Altair bar chart
    """

    kaiju_raw = pd.read_csv(
        file,
        sep="\t",
        header=None,
        usecols=[1, 2, 7],
        names=["name", "taxon_id", "taxonomy"],
    )

    kaiju = (
        kaiju_raw.dropna()
        .assign(taxonomy=lambda x: x.taxonomy.str.split(";").str[:-1])
        .assign(last_level=lambda x: x.taxonomy.str[-1])
        .assign(second_level=lambda x: x.taxonomy.str[-2])
        .assign(third_level=lambda x: x.taxonomy.str[-3])
        .assign(
            kingdom=lambda x: np.select(
                [x.taxonomy.str[0] != "cellular organisms"],
                [x.taxonomy.str[0]],
                default=x.taxonomy.str[1],
            )
        )
        .drop(columns="taxonomy")
    )

    return (
        alt.Chart(kaiju, title="Kaiju classification on MEGAHIT contigs")
        .mark_bar()
        .encode(
            alt.X(
                "count(last_level):Q",
                title="Number of occurancies",
                axis=alt.Axis(values=np.arange(0, 200, 1), format=".0f"),
            ),
            alt.Y(
                "last_level:N",
                sort=alt.EncodingSortField(
                    field="last_level:N", op="count", order="descending"
                ),
                title=None,
            ),
            alt.Color("last_level:N", title=None),
            tooltip=[alt.Tooltip("second_level"), alt.Tooltip("third_level")],
        )
        .properties(width=500, height=500)
    )
