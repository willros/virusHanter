from pathlib import Path
import shutil
from typing import Union


def copy_reports(
    source_file: Union[str, Path], dest: Union[str, Path], name: str
) -> None:
    """
    Copies the fastp report to the template folder, in order to include it in the template.
    :param str|Path fastp_file: Path to the fastp report in the result folder
    :param str|Path dest: Path to the destination folder.
    :param str name: Name of the copied file.
    """
    source_file = Path(source_file)
    destination = Path(dest) / name

    shutil.copy(source_file, destination)
