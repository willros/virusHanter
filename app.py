from flask import Flask, render_template, url_for, request, flash, redirect, session
from datetime import timedelta
import altair as alt
import pathlib

from plotting import bracken_raw, contig_quality
from plotting import kaiju_raw

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

    return render_template("home.html")


@app.route("/raw_data", methods=["GET", "POST"])
def raw_data():
    """
    Shows information about the Kaiju and Kraken runs on the raw file (not the contigs)
    """
    # make this into a session variable which the user can switch between in the right sidebar
    current_sample = "sample1"

    # files
    cleaned_bracken_report = f"static/data/{current_sample}/results/cleaned_files/Bat-Guano-15_S6_L001_R_bracken_raw.csv"
    clenaed_kaiju_report = f"static/data/{current_sample}/results/cleaned_files/Bat-Guano-15_S6_L001_R_kaiju_raw.csv"

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
        current_sample=current_sample,
    )


@app.route("/contig_inspection", methods=["GET"])
def contig_inspection():
    """
    Shows information and plots about the quality and coverage of the contigs made by MEGAHIT.
    """
    # make this into a session variable which the user can switch between in the right sidebar
    current_sample = "sample1"

    # files (need to make this automated later)
    megahit_csv = (
        f"static/data/{current_sample}/results/megahit/Bat-Guano-15_S6_L001_R.csv"
    )
    short = f"static/data/{current_sample}/results/plots/short.pdf"
    medium = f"static/data/{current_sample}/results/plots/medium.pdf"
    long = f"static/data/{current_sample}/results/plots/long.pdf"

    # plots
    histogram = contig_quality.megahit_contig_histogram(file=megahit_csv)
    boxplot = contig_quality.megahit_contig_boxplot(file=megahit_csv)
    histo_and_box = alt.hconcat(histogram, boxplot).to_json()

    return render_template(
        "contig_inspection.html",
        histo_and_box=histo_and_box,
        short=short,
        medium=medium,
        long=long,
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)
