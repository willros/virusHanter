import pandas as pd
import altair as alt
import numpy as np

alt.data_transformers.disable_max_rows()


def bar_chart_cat_megahit(file: str) -> alt.vegalite.v4.api.Chart:
    """
    Plots taxonomy abundancy predicted by CAT of the contigs from the megahit assembly.
    Needs the CAT file.
    :param str file: Path to the CAT file made on megahit contigs.
    :return: Altair bar chart
    """
    cat_raw = pd.read_csv(file, sep="\t")

    cat = (
        cat_raw.rename(
            columns={
                "# contig": "name",
                "species": "last_level_cat",
                "genus": "second_level_cat",
                "family": "third_level_cat",
            }
        )
        .loc[lambda x: x.classification != "no taxid assigned"]
        .loc[lambda x: x["superkingdom"] != "no support"]
        .loc[lambda x: x["phylum"] != "no support"]
        .drop(columns=["classification", "lineage", "lineage scores"])
        .assign(kingdom_cat=lambda x: x["superkingdom"].str[:-6])
        .drop(columns="superkingdom")
        .loc[lambda x: x.last_level_cat != "no support"]
    )

    return (
        alt.Chart(cat, title="CAT classification on MEGAHIT contigs")
        .mark_bar()
        .encode(
            alt.X(
                "count(last_level_cat):Q",
                title="Number of occurancies",
                axis=alt.Axis(values=np.arange(0, 200, 1), format=".0f"),
            ),
            alt.Y(
                "last_level_cat:N",
                sort=alt.EncodingSortField(
                    field="last_level_cat:N", op="count", order="descending"
                ),
                title=None,
            ),
            alt.Color("last_level_cat:N", title=None),
            tooltip=[alt.Tooltip("second_level_cat"), alt.Tooltip("third_level_cat")],
        )
        .properties(width=500, height=500)
    )
