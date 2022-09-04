from pathlib import Path
import shutil
from typing import Union

def copy_fastp(fastp_file: Union[str, Path], dest: Union[str, Path]) -> None:
    """
    Copies the fastp report to the template folder, in order to include it in the template.
    :param str|Path fastp_file: Path to the fastp report in the result folder
    :param str|Path dest: Path to the destination folder.
    """
    fastp_file = Path(fastp_file)
    destination = Path(dest) / 'fastp_report.html'
    
    shutil.copy(fastp_file, destination)
    
    
