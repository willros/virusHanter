from flask import Flask, render_template, url_for, request, flash, redirect, session
from datetime import timedelta
import altair as alt
import pathlib
import pandas as pd

from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit

# init app and objects
app = Flask(__name__)

# config
app.secret_key = "70c6a968ed1ada341dbcbf252b3ea3cf"
app.config.SAMPLES = [
    str(x) for x in pathlib.Path("static/data").iterdir() if x.is_dir()
]

# ldap-login
@app.before_request
def make_session_expire():
    """
    Makes the session expire after 5 minutes, if no GET REQUEST is made.
    """
    session.permanent = True
    app.permanent_session_lifetime = timedelta(minutes=5)


# routes
@app.route("/")
@app.route("/home", methods=["GET"])
def home():
    """
    Route to the home page.
    """
    if not "SAMPLE" in session:
        session["SAMPLE"] = app.config.SAMPLES[0]

    return render_template("home.html")


@app.route("/choose_sample", methods=["GET", "POST"])
def choose_sample():
    """
    Stores the sample to analyze in the session variable.
    """
    session["SAMPLE"] = request.form.get("sample_name")

    return render_template("choose_sample.html", samples=app.config.SAMPLES)


@app.route("/raw_data", methods=["GET", "POST"])
def raw_data():
    """
    Shows information about the Kaiju and Kraken runs on the raw file (not the contigs)
    """

    # files
    cleaned_bracken_report = f"{session['SAMPLE']}/results/cleaned_files/Bat-Guano-15_S6_L001_R_bracken_raw.csv"
    clenaed_kaiju_report = f"{session['SAMPLE']}/results/cleaned_files/Bat-Guano-15_S6_L001_R_kaiju_raw.csv"

    # plots
    # the pie chart doesnt work...why?
    bracken_bar_plot = bracken_raw.bar_chart_bracken_raw(cleaned_bracken_report)
    bracken_domain_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, level="domain", virus_only=False
    )
    kaiju_bar_plot = kaiju_raw.bar_chart_kaiju_raw(file=clenaed_kaiju_report).to_json()

    # bracken_pie_chart = bracken_raw.pie_chart_bracken_raw(cleaned_bracken_report).to_json()
    # bar_and_pie = alt.vconcat(bracken_bar_plot, bracken_pie_chart).to_json()
    species_and_domain_bracken = (
        alt.hconcat(bracken_bar_plot, bracken_domain_bar_plot)
        .resolve_scale(color="independent")
        .to_json()
    )

    return render_template(
        "raw_data.html",
        species_and_domain_bracken=species_and_domain_bracken,
        kaiju_bar_plot=kaiju_bar_plot,
        current_sample=session["SAMPLE"],
    )


@app.route("/contig_inspection", methods=["GET"])
def contig_inspection():
    """
    Shows information and plots about the quality and coverage of the contigs made by MEGAHIT.
    """
    # files
    megahit_csv = f"{session['SAMPLE']}/results/megahit/Bat-Guano-15_S6_L001_R.csv"
    short = f"{session['SAMPLE']}/results/plots/short.pdf"
    medium = f"{session['SAMPLE']}/results/plots/medium.pdf"
    long = f"{session['SAMPLE']}/results/plots/long.pdf"

    # plots
    histogram = contig_quality.megahit_contig_histogram(file=megahit_csv)
    boxplot = contig_quality.megahit_contig_boxplot(file=megahit_csv)
    histo_and_box = alt.hconcat(histogram, boxplot).to_json()

    # checkv dataframe
    checkv_tsv = f"{session['SAMPLE']}//results/checkv/Bat-Guano-15_S6_L001_R_quality_summary.tsv"
    checkv_df = (
        pd.read_csv(checkv_tsv, sep="\t")
        .sort_values("contig_length", ascending=False)
        .drop(
            columns=["warnings", "kmer_freq", "completeness_method", "miuvig_quality"]
        )
    )

    return render_template(
        "contig_inspection.html",
        histo_and_box=histo_and_box,
        short=short,
        medium=medium,
        long=long,
        checkv_df=checkv_df.to_html(),
    )


@app.route("/contig_data", methods=["GET"])
def contig_data():
    """
    Shows analysis of the contig data.
    """

    # files
    kaiju_megahit_report = f"{session['SAMPLE']}/results/kaiju/megahit/Bat-Guano-15_S6_L001_R_names_megahit.out"
    cat_megahit_out = (
        f"{session['SAMPLE']}/results/cat/CAT_Bat-Guano-15_S6_L001_R_contigs_names.txt"
    )

    # plots
    kaiju_bar_plot = kaiju_megahit.bar_chart_kaiju_megahit(file=kaiju_megahit_report)
    cat_bar_plot = cat_megahit.bar_chart_cat_megahit(file=cat_megahit_out)
    kaiju_and_cat = (
        alt.hconcat(kaiju_bar_plot, cat_bar_plot)
        .resolve_scale(color="independent")
        .to_json()
    )

    # test dataframe
    csv = f"{session['SAMPLE']}/results/cleaned_files/Bat-Guano-15_S6_L001_R_cat_kaiju_merged.csv"
    df = pd.read_csv(csv)[
        ["name", "taxon_id", "length", "last_level_kaiju", "last_level_cat"]
    ]

    return render_template(
        "contig_data.html",
        kaiju_and_cat=kaiju_and_cat,
        current_sample=session["SAMPLE"],
        df=df.to_html(),
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)
