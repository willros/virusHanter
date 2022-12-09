import re

def parse_alignments(log_file: str) -> tuple[int, float]:
    """Parses out the total number of reads and the percentage of reads that were aligned to the reference genome from a bowtie2 log file.
    
    Args:
        log_file: The path to the log file.
    
    Returns:
        A tuple containing the total number of reads as an int and the percentage of reads that were aligned to the
        reference genome as a float, or None if the percentage could not be parsed.
    """
    total_reads = 0
    percent_aligned = None
    pattern_total = r"^(\d+) reads; of these:$"
    pattern_aligned = r"^[\d.]+\% overall alignment rate$"
    with open(log_file) as f:
        for line in f:
            match = re.search(pattern_total, line)
            if match:
                total_reads = int(match.group(1))
            match = re.search(pattern_aligned, line)
            if match:
                percent_aligned = float(match.group(0).split("%")[0])
    return total_reads, percent_aligned

