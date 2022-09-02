from flask import Flask, render_template, url_for, request, flash, redirect, session
from datetime import timedelta
import altair as alt

from plotting import bracken_raw, contig_quality

# init app and objects
app = Flask(__name__)

# config
app.secret_key = "70c6a968ed1ada341dbcbf252b3ea3cf"


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

    return render_template("raw_data.html")


@app.route("/contig_inspection", methods=["GET"])
def contig_inspection():
    """
    Shows information and plots about the quality and coverage of the contigs made by MEGAHIT.
    """
    # files (need to make this automated later)
    megahit_csv = "data/sample1/results/megahit/Bat-Guano-15_S6_L001_R.csv"
    # mpileup = "data/sample1/results/cleaned_files/Bat-Guano-15_S6_L001_R_contigs_coverage_mpileup.csv"

    # plots
    histogram = contig_quality.megahit_contig_histogram(file=megahit_csv).to_json()
    # facet_wraps = contig_quality.megahit_contig_coverage_facet(file=mpileup)
    # short = facet_wraps[0].to_json()

    return render_template("contig_inspection.html", histogram=histogram)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)
