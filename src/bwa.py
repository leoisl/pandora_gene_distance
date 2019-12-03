import subprocess
from typing import Tuple, List

import pysam
from pathlib import Path


class BWA:
    def __init__(self, threads=1):
        self.threads = threads
        self.reference = ""

    def index(self, reference: str):
        self.reference = reference

        completed_process = subprocess.run(
            ["bwa", "index", self.reference],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        completed_process.check_returncode()

    def align(self, query: str) -> Tuple[str, str]:
        options = self.get_options()
        bwa_mem = subprocess.run(
            ["bwa", "mem", *options, str(self.reference), "-"],
            input=query.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if bwa_mem.returncode != 0:
            if b"fail to locate the index" in bwa_mem.stderr:
                raise IndexError("Reference must be indexed by BWA before alignment.")
            else:
                bwa_mem.check_returncode()

        return bwa_mem.stdout.decode(), bwa_mem.stderr.decode()

    def get_options(self):
        options = []
        options.extend(["-t", str(self.threads), "-k", "8", "-T", "15"])

        return options

    @staticmethod
    def parse_sam_string(
        sam_string: str
    ) -> Tuple[pysam.VariantHeader, List[pysam.AlignedSegment]]:
        header = ""
        sam_lines = []
        for line in sam_string.split("\n"):
            if line.startswith("@"):
                header += line + "\n"
            else:
                sam_lines.append(line)

        header = pysam.AlignmentHeader.from_text(header)

        return (
            header,
            [pysam.AlignedSegment.fromstring(sam, header) for sam in sam_lines if sam],
        )

    @staticmethod
    def map_query_to_ref(
        query: Path, ref: Path, output: Path, threads: int = 1
    ) -> Tuple[pysam.VariantHeader, List[pysam.AlignedSegment]]:
        bwa = BWA(threads)
        bwa.reference = str(ref)
        stdout, stderr = bwa.align(query.read_text())

        # write sam to file if output path given
        if output.name:
            output.write_text(stdout)

        return bwa.parse_sam_string(stdout)
