from flask import Flask, render_template, url_for, request, flash, redirect, session
from datetime import timedelta
import altair as alt
from pathlib import Path
import pandas as pd

from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit
from helper_functions import move_files

# init app and objects
app = Flask(__name__)

# config
app.secret_key = "70c6a968ed1ada341dbcbf252b3ea3cf"
app.config.SAMPLES = [str(x) for x in Path("static/data").iterdir() if x.is_dir() and not x.stem.startswith(".")]

# ldap-login
@app.before_request
def make_session_expire():
    """
    Makes the session expire after 5 minutes, if no GET REQUEST is made.
    """
    session.permanent = True
    app.permanent_session_lifetime = timedelta(minutes=5)

    # choose the last sample if sample is not chosen.
    if not "SAMPLE" in session:
        session["SAMPLE"] = app.config.SAMPLES[-1]


# routes
@app.route("/")
@app.route("/home", methods=["GET"])
def home():
    """
    Route to the home page.
    """

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
    # number of species to show in the bar plots
    number = 10  # default
    if request.method == "POST":
        number = int(request.form.get("number"))

    # files
    cleaned_bracken_report = list(Path(session["SAMPLE"]).rglob("*bracken_raw.csv"))[0]
    cleaned_kaiju_report = list(Path(session["SAMPLE"]).rglob("*kaiju_raw.csv"))[0]

    # plots
    bracken_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, number=number
    )
    bracken_domain_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, level="domain", virus_only=False
    )
    kaiju_bar_plot = kaiju_raw.bar_chart_kaiju_raw(file=cleaned_kaiju_report).to_json()

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
    megahit_csv = list(
        Path(session["SAMPLE"]).joinpath("megahit").rglob("*.csv")
    )[0]
    
    # pdf placeholder for the samples that lack a certain figure of contig size
    short = Path("static/no-contigs-placeholder.pdf")
    medium = Path("static/no-contigs-placeholder.pdf")
    long = Path("static/no-contigs-placeholder.pdf")
    
    try:
        short = list(Path(session["SAMPLE"]).rglob("*short*.pdf"))[0]
    except:
        pass

    try:
        medium = list(Path(session["SAMPLE"]).rglob("*medium*.pdf"))[0]
    except:
        pass
    
    try:
        long = list(Path(session["SAMPLE"]).rglob("*long*.pdf"))[0]
    except:
        pass
    
    checkv_tsv = list(Path(session["SAMPLE"]).rglob("*quality_summary.tsv"))[0]

    # plots
    histogram = contig_quality.megahit_contig_histogram(file=megahit_csv)
    boxplot = contig_quality.megahit_contig_boxplot(file=megahit_csv)
    histo_and_box = alt.hconcat(histogram, boxplot).to_json()

    # checkv dataframe
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


@app.route("/contig_data", methods=["GET", "POST"])
def contig_data():
    """
    Shows analysis of the contig data.
    """

    # files
    kaiju_megahit_report = list(Path(session["SAMPLE"]).rglob("*megahit.out"))[0]
    cat_megahit_out = list(Path(session["SAMPLE"]).rglob("*contigs_names.txt"))[0]
    cat_kaiju_csv = list(Path(session["SAMPLE"]).rglob("*cat_kaiju_merged.csv"))[0]

    # plots
    kaiju_bar_plot = kaiju_megahit.bar_chart_kaiju_megahit(file=kaiju_megahit_report)
    cat_bar_plot = cat_megahit.bar_chart_cat_megahit(file=cat_megahit_out)
    kaiju_and_cat = (
        alt.hconcat(kaiju_bar_plot, cat_bar_plot)
        .resolve_scale(color="independent")
        .to_json()
    )

    # cat and kaiju dataframe
    df = pd.read_csv(cat_kaiju_csv)

    # choose fasta sequence to look at
    contigs = df.name.to_list()
    if request.method == "POST":
        contig = request.form.get("contig_name")
        sequence = df.loc[lambda x: x.name == contig].sequence.squeeze()
        return render_template("fasta_sequence.html", contig=contig, sequence=sequence)

    return render_template(
        "contig_data.html",
        kaiju_and_cat=kaiju_and_cat,
        current_sample=session["SAMPLE"],
        df=df[
            ["name", "taxon_id", "length", "last_level_kaiju", "last_level_cat"]
        ].to_html(),
        contigs=contigs,
    )


@app.route("/read_quality", methods=["GET"])
def read_quality():
    """
    Shows the content of the report generated by fastp.
    """
    fastp_report = list(Path(session["SAMPLE"]).rglob("*fastp.html"))[0]
    template_folder = "templates/"

    move_files.copy_reports(fastp_report, template_folder, name="fastp_report.html")

    return render_template("read_quality.html")


@app.route("/krona_raw", methods=["GET"])
def krona_raw():
    """
    Displays the krona report for the raw reads.
    """
    # copy the krona chart to the template folder
    krona_raw = list(Path(session["SAMPLE"]).rglob("*krona_raw.html"))[0]
    template_folder = "templates/"
    move_files.copy_reports(krona_raw, template_folder, name="krona_raw.html")

    return render_template("krona_raw.html")


@app.route("/krona_contig", methods=["GET"])
def krona_contig():
    """
    Displays the krona report for the contigs.
    """
    # copy the krona chart to the template folder
    krona_contig = list(Path(session["SAMPLE"]).rglob("*krona_megahit.html"))[0]
    template_folder = "templates/"
    move_files.copy_reports(krona_contig, template_folder, name="krona_contig.html")

    return render_template("krona_contig.html")


# route for testing
@app.route("/test_template", methods=["GET", "POST"])
def test_template():
    """
    For testing different templates
    """
    cat_kaiju_csv = list(Path(session["SAMPLE"]).rglob("*cat_kaiju_merged.csv"))[0]

    df = pd.read_csv(cat_kaiju_csv)

    contigs = df.name.to_list()

    if request.method == "POST":
        contig = request.form.get("contig_name")
        sequence = df.loc[lambda x: x.name == contig].sequence.squeeze()
        return render_template("fasta_sequence.html", contig=contig, sequence=sequence)

    return render_template("test_template.html", contigs=contigs)


if __name__ == "__main__":
    app.run()
