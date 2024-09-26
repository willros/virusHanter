from pathlib import Path
import re
import subprocess
import os
import pandas as pd
import numpy as np
import pyfastx
import janitor
import seaborn as sns
import matplotlib.pyplot as plt
import panel as pn
import altair as alt
from bs4 import BeautifulSoup
import warnings
from io import StringIO
warnings.filterwarnings("ignore")


def read_file_as_blob(file_path):
    with open(file_path, 'rb') as file:
        blob_data = file.read().hex()
    return blob_data

def common_suffix(folder: str) -> str:
    samples = sorted(
        [file.name for file in Path(folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", file.name)]
    )

    suffix = ""
    test = samples[0]
    
    for i in range(1, len(test) + 1):
        index = -i
        if any(sample[index] != test[index] for sample in samples):
            break
        suffix += test[index]
        
    return suffix[::-1]

def paired_reads(folder: str) -> list[str]:
    def common_name(str1: str, str2: str) -> str:
        name = ""
        for a, b in zip(str1, str2):
            if a != b:
                break
            name += a
        return name
    
    samples = sorted(
        [x.stem for x in Path(folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", x.name)]
    )
    
    prefixes = []
    for i in range(0, len(samples), 2):
        read1, read2 = samples[i], samples[i + 1]
        common = common_name(read1, read2)
        prefixes.append(common)
        
    return prefixes

def kaiju_db_files(kaiju_db: str) -> tuple[str, str, str]:
    files = [x for x in Path(kaiju_db).iterdir() if x.is_file()]
    fmi = [x.resolve for x in files if x.suffix == ".fmi"][0]
    names = [x.resolve for x in files if x.name == "names.dmp"][0]
    nodes = [x.resolve for x in files if x.name == "nodes.dmp"][0]
    return fmi, names, nodes

def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1]] for x in fastx]))

    df = (
        pd.DataFrame({"name": reads[0], "sequence": reads[1]})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )
    
    return df

def wrangle_kraken(kraken_file: str) -> pd.DataFrame:
    kraken = (
        pd.read_csv(
            kraken_file, sep="\t", header=None,
            names=["percent", "count_clades", "count", "tax_lvl", "taxonomy_id", "name"]
        )
        .assign(name=lambda x: x.name.str.strip())
        .assign(
            domain=lambda x: np.select(
                [x.tax_lvl.isin(["D", "U", "R"])],
                [x.name],
                default=pd.NA
            )
        )
        .fillna(method="ffill")
    )
    
    return kraken

def run_blastn(contigs: str, db: str, temp_file: str, threads: int) -> pd.DataFrame:
    # In order for blastn to find the files
    os.environ["BLASTDB"] = config["BLASTDB_ENVIRON_VARIABLE"]
    df = pd.read_csv(contigs)
    if df.shape[0] == 0:
        return df
    
    matches = []
    
    for contig in df.itertuples():
        with open(temp_file, "w+") as f:
            print(f">{contig.name}", file=f)
            print(contig.sequence, file=f)
        command = [
            "blastn", "-num_threads", str(threads), "-task", "megablast", "-query", temp_file, "-db", db, "-max_target_seqs", "1",
            "-outfmt", "6 stitle sacc pident slen"
        ]
        
        match = subprocess.check_output(command, universal_newlines=True).strip()
        matches.append(match)
        
    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.shape[0] == 0:
        return df
    
    df[["match_name", "accession", "percent_identity", "sequence_len"]] = (
        df.matches.str.split("\t", expand=True).loc[:, :3]
    )
    df = df.assign(sequence_len=lambda x: [y[0] for y in x.sequence_len.str.split("\n")])
    
    return df

            
def parse_bwa_flagstat(bwa_flagstat: str) -> tuple[int, float]:
    pattern_total = r"(\d*) \+ \d* paired in sequencing"
    pattern_mapped = r"(\d*) \+ \d* with itself and mate mapped"
    #pattern_mapped = r"\d+ \+ \d+ mapped \((\d+\.\d+)%"
    
    with open(bwa_flagstat) as f:
        flagstat = f.read()
        
    total_reads = int(re.findall(pattern_total, flagstat)[0])
    total_mapped = int(re.findall(pattern_mapped, flagstat)[0])
    
    # give the fraction
    total_mapped = (total_mapped / total_reads) * 100
    
    return total_reads, total_mapped


def parse_fastp(fastp_report: str) -> pd.DataFrame:
    summary_information = {}

    with open(fastp_report, "r") as f:
        soup = BeautifulSoup(f, "html.parser")

    for table in soup.find_all("table", class_="summary_table")[:4]:
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            key = cells[0].text
            value = cells[1].text
            summary_information[key] = value
            
    df = pd.DataFrame.from_dict(summary_information, orient="index", columns=["value"])
    df = df.rename_axis("description").reset_index()
    
    return df


def plot_flagstat(flagstat: str):
    total_reads, percent_aligned = parse_bwa_flagstat(flagstat)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned

    # Create dataframe
    df = pd.DataFrame(
        {"amount": [number_unaligned, number_aligned], "type": ["unaligned", "aligned"]}
    )

    # Generate plot
    plot = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                "sum(amount)",
                stack="normalize",
                axis=alt.Axis(format="%"),
                title=None,
            ),
            color=alt.Color("type:N", scale=alt.Scale(scheme='dark2')),
            tooltip=[
                alt.Tooltip("amount:Q", title="Number of reads aligned"),
                alt.Tooltip("type:N"),
            ],
        )
        .properties(
            title="Reads aligned to host",
            width="container",
            height=100,
        )
    )
        
    return plot


def plot_kaiju(
    kaiju: str, cutoff: float = 0.01, max_entries: int = 10
):
    df = (
        pd.read_csv(kaiju, sep="\t")
        .assign(percent=lambda x: x.percent / 100)
    )
    
    unclassified = df.loc[lambda x: x.taxon_name == "unclassified"].percent.squeeze()
    
    df = (
        df
        .drop(columns=["file"])
        .sort_values("percent", ascending=False)
        .loc[lambda x: x.taxon_name != "unclassified"]
        .loc[lambda x: x.percent > cutoff]
        .head(max_entries)
    )
    
    plot= (
        alt.Chart(df) 
        .mark_bar()
        .encode(
            alt.X("percent:Q", 
                  axis=alt.Axis(format="%"), 
                  title=f"Percent of reads ({unclassified * 100:.1f}% not classified)",
                  scale=alt.Scale(zero=True) 
            ),
            alt.Y("taxon_name:N", sort="-x", title=None),
            alt.Color("taxon_name:N", title=None, legend=None),
            tooltip=[alt.Tooltip("taxon_name:N"), alt.Tooltip("reads:Q", title="Number of reads")],
        )
        .properties(
            width="container",
            title=f"""
            Kaiju classification
            """
        )
    ) 
    return plot


def kraken_df(
    kraken: str, 
    level: str = "species", 
    cutoff: float = 0.01, 
    max_entries: int = 10,
    virus_only: bool = True,
) -> pd.DataFrame:
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
        pd.read_csv(kraken)
        .assign(percent=lambda x: x.percent / 100)
    )
    
    unclassified = df.loc[lambda x: x.domain == "unclassified"]
    if unclassified.shape[0] == 0:
        unclassified = 0
    else:
        unclassified = unclassified.percent.squeeze()
    
    df = (
        df
        .loc[lambda x: x.tax_lvl == taxonomy[level]]
        .loc[lambda x: x.percent > cutoff]
        .sort_values("percent", ascending=False)
    )

    if virus_only:
        return df.loc[lambda x: x.domain == "Viruses"].head(max_entries), unclassified
    return df.head(max_entries), unclassified


def plot_kraken(
    kraken: str,
    level: str = "species",
    cutoff: float = 0.001,
    max_entries: int = 10,
    virus_only=True,
):
    
    df, unclassified = kraken_df(
        kraken,
        level, 
        cutoff,
        max_entries,
        virus_only
    )

    return (
        alt.Chart(df, title="Kraken classification")
        .mark_bar()
        .encode(
            alt.X(
                "percent:Q",
                axis=alt.Axis(format="%"),
                title=f"Percent of reads ({unclassified * 100:.1f}% not classified)"
            ),
            alt.Y("name:N", sort="-x", title=None),
            alt.Color("name:N", title=None, legend=None),
            tooltip=[alt.Tooltip("domain:N"), alt.Tooltip("count_clades:Q", title="Number of reads")],
        )
        .properties(width="container")
    )

def plot_blastn(blastn):
    df = pd.read_csv(blastn)
    
    if df.shape[0] == 0:
        return alt.Chart(df, title="No classified contigs").mark_bar()
    
    plot = (
        alt.Chart(df, title="BLASTN of Contigs")
        .mark_bar()
        .encode(
            alt.X(
                "count(match_name):Q",
                title="Number of contigs",
                axis=alt.Axis(format="d", tickMinStep=1),
            ),
            alt.Y("match_name:N", sort="-x", title=None),
            alt.Color("match_name:N", title=None, legend=None),
        )
        .properties(
            width="container",
        )
    )
    
    return plot


def alignment_stats(flagstat_log, species: str) -> tuple:
    
    total_reads, percent_aligned = parse_bwa_flagstat(flagstat_log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned
    
    # Markdown with above text
    alignment_stats = pn.pane.Markdown(
        f"""
        ### Total Number of Reads: 
        {total_reads:,}
        ### Reads aligned to {species} Genome: 
        {number_aligned:,} ({percent_aligned:.2f}%)
        ### Reads NOT aligned to {species} Genome:
        {number_unaligned:,} ({100 - percent_aligned:.2f}%)
        """,
        name=f"{species} Alignment Stats"
    )
    
    # human bwa flagstat plot
    flagstat_plot = plot_flagstat(flagstat_log).interactive()
    flagstat_pane = pn.pane.Vega(
        flagstat_plot, sizing_mode="stretch_both", name=f"{species} Alignment Plot"
    )
    
    return alignment_stats, flagstat_pane

def panel_report(
    result_folder: str,
    blastn_file: str,
    kraken_file: str,
    kaiju_table: str,
    secondary_host: str = None
) -> None:
    result_folder = Path(result_folder)
    sample_name = result_folder.parts[-1]
    
    pn.extension("tabulator")
    pn.extension("vega", sizing_mode="stretch_width", template="fast")
    pn.widgets.Tabulator.theme = 'modern'
    
    def header(
        text: str, 
        bg_color: str = "#04c273",
        height: int = 150,
        fontsize: str = "px20",
        textalign: str = "center"
    ):
        return (
            pn.pane.Markdown(
                f"""
                {text}
                """,
                styles={
                    "color": "white", 
                    "padding": "10px",
                    "text-align": f"{textalign}",
                    "font-size": f"{fontsize}",
                    "background": f"{bg_color}",
                    "margin": f"10",
                    "height": f"{height}",
                }
            )
       )
    
    # --- Alignment and Read Statistics --- #
    
    # human
    flagstat_human_log = list(result_folder.rglob("human_contamination_flagstat.txt"))[0]
    
    human_stats, human_pane = alignment_stats(flagstat_human_log, species="Human")
    
    # If secondary alignment
    flagstat_secondary_log = [
        x for x in result_folder.rglob("*") if x.stem == "secondary_contamination_flagstat"
    ]
    
    if flagstat_secondary_log:
        secondary_stats, secondary_pane = alignment_stats(flagstat_secondary_log[0], species=secondary_host)
    else:
        secondary_stats = None
        
    # fastp report
    fastp_report = list(result_folder.rglob("*FASTP/*.html"))[0]
    
    fastp_df = parse_fastp(fastp_report)
    
    fastp_table = pn.widgets.Tabulator(
        fastp_df, 
        layout='fit_columns',
        show_index=False,
        name="Read Summary from FASTP"
    )
    
    # Header for this section
    alignment_subheader = header(
        text=
        f"""
        ## Alignment and Read statistics
        Reads were aligned to Human and optionally other host species with bwa
        """, 
        bg_color="#04c273",
        height=80, 
        textalign="left"
    )
    
    flagstat_and_stats = pn.Column(
        human_stats,
        pn.layout.Divider(),
        human_pane,
        name="Alignment"
    )
    
    if secondary_stats:
        flagstat_and_stats.extend(
            [
                secondary_stats,
                pn.layout.Divider(),
                secondary_pane,
            ]
        )
    
    # SECTION
    alignment_tab = pn.Tabs(flagstat_and_stats, fastp_table)
    alignment_section = pn.Column(alignment_subheader, alignment_tab)
    
    # --- Raw Classification --- #
    
    kraken_plot = plot_kraken(kraken_file, virus_only=True).interactive()

    kraken_domain_plot = plot_kraken(
        kraken_file, level="domain", virus_only=False
    ).interactive()

    kaiju_plot = plot_kaiju(kaiju_table).interactive()
    
    # Vega panes
    kraken_pane = pn.pane.Vega(
        kraken_plot, sizing_mode="stretch_both", name="Kraken Virus Only"
    )
    kraken_domain_pane = pn.pane.Vega(
        kraken_domain_plot, sizing_mode="stretch_both", name="Kraken All Domains"
    )
    kaiju_pane = pn.pane.Vega(kaiju_plot, sizing_mode="stretch_both", name="Kaiju")
    
    # Header for this section
    raw_header = header(
        text=
        f"""
        ## Classification of Raw Reads
        Reads from the sequences were classfied with kaiju and Kraken2
        """, 
        bg_color="#04c273",
        height=80, 
        textalign="left"
    )
    
    # Section
    raw_tab = pn.Tabs(kraken_pane, kraken_domain_pane, kaiju_pane)
    raw_section = pn.Column(raw_header, raw_tab)
    
    # --- Contig Classification --- #
    
    # plots
    blastn_plot = plot_blastn(blastn_file).interactive()

    # Vega panes
    blastn_pane = pn.pane.Vega(blastn_plot, sizing_mode="stretch_both", name="BLASTN")
    
    # BLASTN dataframe
    blastn_df = pd.read_csv(blastn_file)
    
    if blastn_df.shape[0] == 0:
        blastn_df = pd.DataFrame({"sequence": ["NO SEQUENCES GENERATED"]})
    else:
        blastn_df = blastn_df.drop(columns=["name", "matches"])
        
    blastn_table = pn.widgets.Tabulator(
        blastn_df, 
        editors={
            'sequence': {'type': 'editable', 'value': False}
        },
        layout='fit_columns',
        pagination='local', 
        page_size=15,
        show_index=False,
        name="Contig Table"
    )
    
    # Header for this section
    contig_header = header(
        text=
        f"""
        ## Classification of Contigs
        #### MEGAHIT was used to generate contigs which were classified with BLASTN 
        #### To get the sequence: Copy (CTRL-C) the column containing the sequence
        """, 
        bg_color="#04c273",
        height=120, 
        textalign="left"
    )
    
    # Section
    contig_tab = pn.Tabs(blastn_pane, blastn_table)
    contig_section = pn.Column(contig_header, contig_tab)
    
    
    # --- Coverage plots --- #
    
    coverage_plot_path = result_folder / "COVERAGE_PLOTS"
    coverage_plots = [x for x in coverage_plot_path.iterdir() if x.suffix == ".svg"]
    
    coverage_tab = pn.Tabs()
    if coverage_plots:
        for plot in coverage_plots:
            name = plot.stem[:20]
            SVG_pane = pn.pane.SVG(plot, sizing_mode='stretch_width', name=name)
            coverage_tab.append(SVG_pane)
    else:
        no_plots = pn.pane.Markdown(
            f"""
            ## No Coverage plots Available
            """, 
            name="No Coverage Plots"
        )
        coverage_tab.append(no_plots)
    
    # Header for this section
    coverage_header = header(
        text=
        f"""
        ## Alignment Coverage
        """, 
        bg_color="#04c273",
        height=80, 
        textalign="left"
    )
    
    # Section
    coverage_section = pn.Column(coverage_header, coverage_tab)
    
    # --- Create the report --- #
    
    # header
    head = header(
        text=
        f"""
        # Virushanter report
        ## Report of {sample_name}
        """,
        fontsize="20px",
        bg_color="#011a01",
        height=185
    )
    
    all_tabs = pn.Tabs(
        ("Alignment Stats", alignment_section),
        ("Classification of Raw Reads", raw_section),
        ("Classification of Contigs", contig_section),
        ("Alignment Coverage", coverage_section),
        tabs_location="left",
    )

    rep = pn.Column(
        head,
        pn.layout.Divider(),
        all_tabs,
    )
    
    return rep

# --- CONFIG AND SETUP
configfile: "/scratch/virusHanter/config.yaml"

SAMPLES = paired_reads(config["SAMPLES"])
SUFFIX = common_suffix(config["SAMPLES"])

SAMPLES_FOLDER = config["SAMPLES"]
RESULT_FOLDER = f"{config['RESULTS_FOLDER']}/{Path(SAMPLES_FOLDER).parts[-1]}" 

# --- RUN EVERYTHING
clean_list = [f"{RESULT_FOLDER}/analysis_done.txt"] if config["CLEAN"] == "TRUE" else [] 
rule all:
    input:
        expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),

        # aggregation
        f"{RESULT_FOLDER}/run_information_{Path(SAMPLES_FOLDER).parts[-1]}.csv",
        # clean folder after completion
        clean_list,

rule fastp:
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/fastp.log"
    input:
        r1=f"{SAMPLES_FOLDER}/{{sample}}1{SUFFIX}",
        r2=f"{SAMPLES_FOLDER}/{{sample}}2{SUFFIX}",
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r1_trimmed.fq",
        r2=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r2_trimmed.fq",
        html_report=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.html",
        json_report=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.json",
    threads: 
        config["THREADS"]
    shell:
        """
        fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1} \
        --out2 {output.r2} \
        --report_title {wildcards.sample} \
        --thread {threads} \
        --html {output.html_report} \
        --json {output.json_report} \
        > {log} 2>&1
        """

rule bwa_human:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    params:
        index=config["HUMAN_INDEX"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_human_contamination.log"
    threads: 
        config["THREADS"]
    output:
        mapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human.bam",
    shell:
        """
        # NB: added -k 26
        bwa mem -t {threads} -k 26 {params.index} {input.r1} {input.r2} > {output} 2> {log}
        """
        # bowtie2 --threads {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.mapped} 2> {log} 
        
rule remove_host:
    input:
        mapped_reads=rules.bwa_human.output.mapped
    threads: 
        config["THREADS"]
    output:
        unmapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_removed.bam",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/human_contamination_flagstat.txt",
        stats=f"{RESULT_FOLDER}/{{sample}}/logs/human_contamination_stats.txt",
    shell:
        """
        samtools flagstat {input.mapped_reads} > {output.flagstat}
        samtools stats {input.mapped_reads} > {output.stats}
        
        samtools view -f 12 -O bam {input.mapped_reads} > {output.unmapped}
        """
        
rule bam_to_fastq_human:
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bamtofastq_human.log"
    input:
        unmapped=rules.remove_host.output.unmapped,
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_trimmed_human_1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_trimmed_human_2.fastq",
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} {input.unmapped} > {log} 2>&1
        """
        
SECONDARY_HOST = config["SECONDARY_HOST_INDEX"] if "SECONDARY_HOST_INDEX" in config.keys() else ""
rule bwa_secondary_host:
    input:
        r1=rules.bam_to_fastq_human.output.r1,
        r2=rules.bam_to_fastq_human.output.r2,
    params:
        index=SECONDARY_HOST
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_secondary_host_contamination.log"
    threads: 
        config["THREADS"]
    output:
        mapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_secondary.bam",
    shell:
        """
        # bowtie2 --threads {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.mapped} 2> {log} 
        bwa mem -t {threads} {params.index} {input.r1} {input.r2} > {output} 2> {log}
        """
        
rule remove_secondary_host:
    input:
        mapped_reads=rules.bwa_secondary_host.output.mapped
    threads: 
        config["THREADS"]
    output:
        unmapped=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_secondary_removed.bam",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/secondary_contamination_flagstat.txt",
        stats=f"{RESULT_FOLDER}/{{sample}}/logs/secondary_contamination_stats.txt"
    shell:
        """
        samtools flagstat {input.mapped_reads} > {output.flagstat}
        samtools stats {input.mapped_reads} > {output.stats}
        
        samtools view -f 12 -O bam {input.mapped_reads} > {output.unmapped}
        """

rule bam_to_fastq_secondary:
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bamtofastq_secondary.log"
    input:
        unmapped=rules.remove_secondary_host.output.unmapped,
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_trimmed_secondary_1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_trimmed_secondary_2.fastq",
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} {input.unmapped} > {log} 2>&1
        """

        
SECONDARY_HOST_OR_NOT = True if "SECONDARY_HOST_INDEX" in config.keys() else False
FASTQ_1 = (
    rules.bam_to_fastq_human.output.r1 if not SECONDARY_HOST_OR_NOT
    else rules.bam_to_fastq_secondary.output.r1
)
FASTQ_2 = (
    rules.bam_to_fastq_human.output.r2 if not SECONDARY_HOST_OR_NOT
    else rules.bam_to_fastq_secondary.output.r2
)


fmi, names, nodes = kaiju_db_files(config["KAIJU_DB"])
rule kaiju:
    input:
        r1=FASTQ_1,
        r2=FASTQ_2,
        fmi=fmi,
        nodes=nodes,
    threads: 
        config["THREADS"]
    output:
        f"{RESULT_FOLDER}/{{sample}}/KAIJU/{{sample}}.kaiju.out"
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/kaiju.log"
    shell:
        """
        kaiju \
        -t {input.nodes} \
        -f {input.fmi} \
        -i {input.r1} \
        -j {input.r2} \
        -v \
        -z {threads} \
        -o {output} \
        2> {log}
        """

rule kaiju_to_table:
    input:
        names=names,
        nodes=nodes,
        kaiju=rules.kaiju.output,
    output:
        table=f"{RESULT_FOLDER}/{{sample}}/KAIJU/{{sample}}.kaiju.table.tsv",
    shell:
        """
        kaiju2table \
        -t {input.nodes} \
        -n {input.names} \
        -r genus \
        -e \
        -o {output.table} \
        {input.kaiju}
        """
        
rule kraken:
    input:
        r1=FASTQ_1,
        r2=FASTQ_2,
        db=config["KRAKEN_DB"],
    threads: 
        config["THREADS"]
    params:
        output_folder=f"{RESULT_FOLDER}/{{sample}}/KRAKEN"
    output:
        tsv=f"{RESULT_FOLDER}/{{sample}}/KRAKEN/{{sample}}.kraken.tsv"
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/kraken.log"
    shell:
        """
        kraken2 \
        --db {input.db} \
        --threads {threads} \
        --report {output} \
        --use-names \
        --paired \
        {input.r1} \
        {input.r2} \
        2> {log}
        """

rule wrangle_kraken:
    input:
        kraken=rules.kraken.output.tsv
    output:
        csv=f"{RESULT_FOLDER}/{{sample}}/KRAKEN/{{sample}}.kraken.csv"
    run:
        df = wrangle_kraken(input.kraken)
        df.to_csv(output.csv, index=False)


rule megahit:
    input:
        r1=FASTQ_1,
        r2=FASTQ_2,
    output:
        # make temp
        contigs=f"{RESULT_FOLDER}/{{sample}}/MEGAHIT/{{sample}}.contigs.fa",
    params:
        f"{RESULT_FOLDER}/{{sample}}/MEGAHIT"
    threads: 
        config["THREADS"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/megahit.log"
    run:
        # megahit cannot continue if the dir exists
        shell("rm -rf {params}")
        shell("megahit -1 {input.r1} -2 {input.r2} -o {params} --out-prefix {wildcards.sample} 2> {log}")
        # remove everything except contigs file
        shell("ls -d -1 {params}/* | grep -v .fa | xargs rm -rf")
        
        # if contigs.fa is empty
        from pathlib import Path
        if Path(output.contigs).read_text() == "":
            with open(output.contigs, "w+") as f:
                print(">DUMMY_CONTIG", file=f)
                print("TTAACCTTGG" * 20, file=f)

rule pilon:
    input:
        contigs=rules.megahit.output.contigs,
        r1=FASTQ_1,
        r2=FASTQ_2,
    output:
        contigs_bam=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_contigs.bam",
        improved_contigs=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_improved_contigs.fasta"
    params:
        index_folder=f"{RESULT_FOLDER}/{{sample}}/PILON/bwa",
        pilon_folder=f"{RESULT_FOLDER}/{{sample}}/PILON",
    threads: 
        config["THREADS"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_contigs.log"
    run:
        # bwa index
        # need to create dir
        shell("rm -rf {params.index_folder}")
        shell("mkdir {params.index_folder}")
        index_folder = f"{params.index_folder}/{wildcards.sample}"
        shell("bwa index -p {index_folder} {input.contigs}")

        # map reads to index of contigs
        shell("bwa mem -t {threads} {index_folder} {input.r1} {input.r2} 2> {log} | samtools view -h -O bam | samtools sort -o {output.contigs_bam}")
        shell("samtools index {output.contigs_bam}")
        # pilon
        shell("pilon -Xmx50G --threads {threads} --genome {input.contigs} --frags {output.contigs_bam} --outdir {params.pilon_folder} --output {wildcards.sample}_improved_contigs")

        # remove the index
        shell("rm -rf {params.index_folder}")


rule wrangle_pilon:
    input:
        contigs=rules.pilon.output.improved_contigs
    output:
        csv=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}.contigs.csv"
    params:
        min_len=config["CONTIG_LENGTH"],
    run:
        df = fastx_file_to_df(input.contigs)
        df = df.assign(sample_id=wildcards.sample)
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.csv, index=False)


rule blastn:
    input:
        contigs=rules.wrangle_pilon.output.csv
    output:
        blast=f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.csv"
    params:
        temp_file=f"{RESULT_FOLDER}/{{sample}}/BLASTN/temp.blastn",
        db=config["BLASTN_DB"],
        threads=config["THREADS"]
    run:
        df = run_blastn(contigs=input.contigs, db=params.db, temp_file=params.temp_file, threads=params.threads)
        df.to_csv(output.blast, index=False)
        
        if Path(params.temp_file).exists():
            os.remove(params.temp_file)

rule checkv:
    input:
        contigs=rules.pilon.output.improved_contigs
    params:
        db=config["CHECKV_DB"],
        folder=f"{RESULT_FOLDER}/{{sample}}/CHECKV",
    output:
        checkv=f"{RESULT_FOLDER}/{{sample}}/CHECKV/{{sample}}.contamination.tsv"
    threads: 
        config["THREADS"]
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/checkv.log"
    shell:
        """
        checkv contamination \
        -d {params.db} \
        {input.contigs} \
        {params.folder} \
        -t {threads} \
        2> {log}

        mv {params.folder}/contamination.tsv {output.checkv}
        # remove everything except contamination file
        ls -d -1 {params.folder}/* | grep -v .tsv | xargs rm -rf 
        """   
        
rule merge_checkv_blastn:
    input:
        checkv=rules.checkv.output.checkv,
        blastn=rules.blastn.output.blast,
    output:
        merged=f"{RESULT_FOLDER}/{{sample}}/CHECKV/{{sample}}.merged.csv"
    run:
        blastn = pd.read_csv(input.blastn)
        checkv = (
            pd.read_csv(
                input.checkv,
                sep="\t"
            )
            .rename(columns={"contig_id": "name"})
            [["name", "total_genes", "viral_genes", "host_genes", "provirus"]]
        )

        merged = pd.merge(blastn, checkv)
        merged.to_csv(output.merged, index=False)
    
        

rule bwa_align_to_kraken_hits:
    input:
        kraken=rules.wrangle_kraken.output.csv,
	r1=FASTQ_1,
	r2=FASTQ_2,
    output:
        fasta=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/kraken.fasta",
        bam=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/kraken.bam",
    params:
        virus_db=config["VIRUS_PARQUET"],
        index_folder=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/bwa",
    run:
        kraken_list = (
            pd.read_csv(input.kraken)
            .loc[lambda x: x.domain == "Viruses"]
            .sort_values("percent", ascending=False)
            .head(20)
            ["taxonomy_id"]
            .to_list()
        )
        virus_from_kraken = (
            pd.read_parquet(params.virus_db)
            .loc[lambda x: x.tax_id.isin(kraken_list)]
        )
        
        with open(output.fasta, "a+") as f:
            for x in virus_from_kraken.itertuples():
                print(f">{x.name.strip()}", file=f)
                print(x.sequence, file=f)
                
        # create bwa index
        shell("rm -rf {params.index_folder}")
        shell("mkdir {params.index_folder}")
        index_folder = f"{params.index_folder}/{wildcards.sample}"
        shell("bwa index -p {index_folder} {output.fasta}")
        
        # map reads to index of contigs
        shell("bwa mem -t {threads} {index_folder} {input.r1} {input.r2} | samtools sort -o {output.bam} -")
        shell("samtools index {output.bam}")
        
        # remove the index
        shell("rm -rf {params.index_folder}")
        
        
rule bam2plot:
    input:
        bam=rules.bwa_align_to_kraken_hits.output.bam,
    params:
        threshold=config["PLOT_THRESHOLD"],
        n_refs=config["NUMBER_OF_PLOTS"],
        out_folder=f"{RESULT_FOLDER}/{{sample}}/COVERAGE_PLOTS",
    output:
        out_folder=f"{RESULT_FOLDER}/{{sample}}/COVERAGE_PLOTS",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bam2plot.log",
    shell:
        """
        bam2plot from_bam -b {input.bam} -o {params.out_folder} -t {params.threshold} -p svg -n {params.n_refs} > {log}
        """
        
SECONDARY_HOST_NAME = config["SECONDARY_HOST_NAME"] if "SECONDARY_HOST_NAME" in config.keys() else ""
rule generate_report:
    input:
        flagstat=rules.remove_host.output.flagstat,
        plot=rules.bam2plot.params.out_folder,
        blastn_file=rules.merge_checkv_blastn.output.merged,
        kraken_file=rules.wrangle_kraken.output.csv,
        kaiju_table=rules.kaiju_to_table.output.table,
    params:
        result_folder=f"{RESULT_FOLDER}/{{sample}}"
    output:
        report=f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html"
    run:
        panel_report(
            result_folder=params.result_folder,
            blastn_file=input.blastn_file,
            kraken_file=input.kraken_file,
            kaiju_table=input.kaiju_table,
            secondary_host=SECONDARY_HOST_NAME,
        ).save(output.report, title=f"Report of {wildcards.sample}")


rule aggregate_run_information:
    input:
        reports=expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
    params:
        results_folder=f"{RESULT_FOLDER}"
    output:
        run_information=f"{RESULT_FOLDER}/run_information_{Path(SAMPLES_FOLDER).parts[-1]}.csv",
    run:
        def aggregate_information_sample(folder):
            folder = Path(folder).resolve()

            run_name = folder.parts[-2]
            date = run_name.split("_")[0]
            sample_name = folder.parts[-1]
            html_report = read_file_as_blob(list(folder.rglob("REPORT/*.html"))[0])
            kaiju_report = read_file_as_blob(list(folder.rglob("KAIJU/*.tsv"))[0])
            blastn_report = read_file_as_blob(list(folder.rglob("BLASTN/*.csv"))[0])

            df = parse_fastp(
                list(folder.rglob("FASTP/*html"))[0]
            )

            sequencing_length = re.findall(r"\d*", df.iloc[2].value)[0]
            number_reads = df.iloc[6].value

            human_contamination = parse_bwa_flagstat(
                list(folder.rglob("*/human_contamination_flagstat.txt"))[0]
            )

            classified_by_kraken = (
                pd.read_csv(list(folder.rglob("KRAKEN/*.csv"))[0])
                .loc[lambda x: x.name == "Viruses"]
                .percent.squeeze()
            )

            kaiju_df = (
                pd.read_csv(list(folder.rglob("KAIJU/*.tsv"))[0], sep="\t")
                .assign(name_and_count=lambda x: [f"{y.taxon_name} ({y.reads})" for y in x.itertuples()])
            )

            classified_by_kaiju = (
                kaiju_df
                .dropna()
                .percent.sum()
            )
            
            top_virus_kaiju = "||".join(kaiju_df.name_and_count.head(10).to_list())

            blastn_df = (
                pd.read_csv(list(folder.rglob("BLASTN/*blastn.csv"))[0])
                .assign(name_and_len=lambda x: [f"{y.match_name} ({y.read_len})" for y in x.itertuples()])
            )

            number_contigs = (
                blastn_df 
                .shape[0]
            )

            top_contigs_blastn = "||".join(blastn_df.name_and_len.head(5).to_list())

            df = (
                pd.DataFrame()
                .assign(run_name=[run_name])
                .assign(sample_name=[sample_name])
                .assign(date=[date])
                .assign(read_len=[sequencing_length])
                .assign(number_reads=[number_reads])
                .assign(mapped_to_human_percent=[human_contamination[1]])
                .assign(kraken_virus_percent=[classified_by_kraken])
                .assign(kaiju_virus_percent=[classified_by_kaiju])
                .assign(number_of_contigs=[number_contigs])
                .assign(top_contigs_blastn=[top_contigs_blastn])
                .assign(top_virus_kaiju=[top_virus_kaiju])
                .assign(html_report=[html_report])
                .assign(kaiju_report=[kaiju_report])
                .assign(blastn_report=[blastn_report])
            )
            
            return df

        def aggregate_information_whole_run(folder):
            folder = Path(folder)
            
            sample_folders = [x for x in folder.iterdir() if x.is_dir()]
            dfs = []
            for sample_folder in sample_folders:
                df = aggregate_information_sample(sample_folder)
                dfs.append(df)
                
            return pd.concat(dfs, ignore_index=True)

        df = aggregate_information_whole_run(params.results_folder)
        df.to_csv(output.run_information, index=False)


rule clean_everything:
    input:
        rules.aggregate_run_information.output.run_information,
    params:
        results_folder=f"{RESULT_FOLDER}"
    output:
        done=f"{RESULT_FOLDER}/analysis_done.txt"
    shell:
        """
        find {params.results_folder} -type f | grep -v "logs\|.html\|.csv\|.table.tsv\|.tar.gzip\|flagstat" | xargs rm -rf

        echo "Analysis done: $(date +%y%m%d)" > {output.done} 
        """
