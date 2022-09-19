import pandas as pd
import altair as alt
import numpy as np

alt.data_transformers.disable_max_rows()


def megahit_contig_histogram(file: str) -> alt.vegalite.v4.api.Chart:
    """
    Plots histogram of the contigs from the megahit assembled contigs.
    :param str file: Path to the csv file for the megahit contigs.
    :return: Altair histogram
    """
    contigs = pd.read_csv(file)

    return (
        alt.Chart(contigs, title="Megahit contigs size")
        .mark_bar()
        .encode(
            alt.X("length:Q", bin=alt.Bin(step=500), title="Length (nt)"),
            alt.Y("count(length):Q", title="Number of contigs"),
        )
        .properties(width=500, height=500)
    )


def megahit_contig_boxplot(file: str) -> alt.vegalite.v4.api.Chart:
    """
    Returns boxplot of the contigs from the megahit assembled contigs.
    :param str file: Path to the csv file for the megahit contigs.
    :return: Altair boxplot
    """
    contigs = pd.read_csv(file).assign(group="group1")

    return (
        alt.Chart(contigs, title="Boxplot of contigs in Megahit")
        .mark_boxplot(extent="min-max")
        .encode(
            alt.X("group:N", axis=alt.Axis(labels=False, title="Contigs", ticks=False)),
            alt.Y("length:Q", title="Length (nt)"),
        )
        .properties(width=500, height=500)
    )


def megahit_contig_coverage_facet(file: str) -> list[alt.vegalite.v4.api.Chart]:
    """
    Returns a list of plots for every contig in the csv file generated from samtools mpileup and the
    wrangle_contig_info.py script.
    :params str file: Path to the samtools mpileup file in the cleaned_files folder.
    :returns: List of altair facet wrap with plots in the categories: short, medium and long.
    """
    plots = []

    contigs_coverage = (
        pd.read_csv(file)
        .assign(
            category=lambda x: pd.cut(
                x.length,
                bins=[1, 500, 2500, np.inf],
                labels=["short", "medium", "long"],
            )
        )
        .assign(
            category=lambda x: pd.Categorical(x.category, ["short", "medium", "long"])
        )
        .sort_values("category")
    )

    step_size = {"short": 50, "medium": 150, "long": 2000}
    colors = {"short": "blue", "medium": "green", "long": "purple"}

    for contig in ["short", "medium", "long"]:

        base = (
            alt.Chart(contigs_coverage.loc[lambda x, contig=contig: x.category == contig])
            .mark_line(color=colors[contig])
            .encode(
                alt.X(
                    "position:N",
                    axis=alt.Axis(values=np.arange(0, 20000, step_size[contig])),
                ),
                alt.Y("coverage:Q"),
                alt.Facet("contig:N"),
            )
            .properties(width=300, height=300, title=f"Coverage of {contig} contigs")
            .resolve_scale(y="independent", x="independent")
        )

        plots.append(base)

    return plots


# plot with the mean value as red line. Does not work in Flask.

# def megahit_contig_coverage_facet(file: str) -> list[alt.vegalite.v4.api.Chart]:
#    """
#    Returns a list of plots for every contig in the csv file generated from samtools mpileup and the
#    wrangle_contig_info.py script.
#    :params str file: Path to the samtools mpileup file in the cleaned_files folder.
#    :returns: List of altair facet wrap with plots in the categories: short, medium and long.
#    """
#    plots = []
#
#    contigs_coverage = (
#        pd.read_csv(file)
#        .assign(
#            category=lambda x: pd.cut(
#                x.length,
#                bins=[1, 500, 2500, np.inf],
#                labels=["short", "medium", "long"],
#            )
#        )
#        .assign(
#            category=lambda x: pd.Categorical(x.category, ["short", "medium", "long"])
#        )
#        .sort_values("category")
#    )
#
#    step_size = {"short": 50, "medium": 150, "long": 2000}
#
#    for contig in contigs_coverage.category.unique():
#
#        base = (
#            alt.Chart()
#            .mark_line()
#            .encode(
#                alt.X(
#                    "position:N",
#                    axis=alt.Axis(values=np.arange(0, 20000, step_size[contig])),
#                ),
#                alt.Y("coverage:Q"),
#            )
#            .properties(width=125, height=125)
#        )
#
#        rule = (
#            alt.Chart()
#            .mark_rule(color="red")
#            .encode(alt.Y("mean(coverage):Q"))
#            .properties(width=125, height=125)
#        )
#
#        layered = (
#            alt.layer(
#                base, rule, data=contigs_coverage.loc[lambda x: x.category == contig]
#            )
#            .encode(alt.Y(title="Coverage"), alt.X(title="Position"))
#            .facet(
#                "contig:N", columns=6, title=f"Coverage of contigs with {contig} length"
#            )
#            .resolve_scale(y="independent", x="independent")
#        )
#
#        plots.append(layered)
#
#    return plots
