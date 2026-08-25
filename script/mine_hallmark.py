#!/usr/bin/env python3
"""
mine_hallmark.py - Mine one NCLDV hallmark protein family from protein or genome FASTA files.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

HELP_TEXT = """\
mine_hallmark.py - Mine one NCLDV hallmark from protein or genome FASTA files.

Usage: mine_hallmark.py <protein/genome FASTA | input_directory> -db <one_marker.hmm> [OPTIONS]

Mandatory:
  input                 One protein/genome FASTA file, or a directory containing FASTA files
  -db, --db             One specific NCLDV hallmark HMM file, for example D5.hmm, RNR.hmm, mRNAc.hmm

Optionals:
  -o, --outdir          Output directory                                    [default: ./HallmarkMine_<YYYYMMDD>]
  -e, --evalue          Domain i-Evalue cutoff                              [default: 1e-10]
  -c, --hmm-coverage    Minimum HMM coverage, 0-1 scale                     [default: 0.80]
  -t, --threads         CPU threads for hmmsearch                           [default: 4]
  --input-type          Input sequence type: auto, protein, or genome       [default: auto]
  --marker              Marker name to write in output headers              [default: inferred from HMM NAME]
  --min-score           Optional minimum domain bit score                   [default: disabled]
  --keep-all            Keep all filtered hits instead of best per genome
  --keep-temp           Keep temporary normalized FASTA and domtblout files
  -h, --help            Show this help and exit

Main output:
  <outdir>/hallmark/<marker>.faa
"""

USAGE_TEXT = "Usage: mine_hallmark.py <protein/genome FASTA | input_directory> -db <one_marker.hmm> [OPTIONS]\n"

DEFAULTS = dict(
    evalue=1e-10,
    hmm_coverage=0.80,
    threads=4,
)

_LOG = logging.getLogger("mine_hallmark")
OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")

def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)

logging.Logger.output = _output


def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.setLevel(logging.DEBUG)
    _LOG.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    _LOG.addHandler(sh)
    if log_path is not None:
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        _LOG.addHandler(fh)

_NATKEY_RE = re.compile(r"(\d+)")

def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in _NATKEY_RE.split(str(s))]

@dataclass
class ProteinRecord:
    internal_id: str
    source_file: str
    original_header: str
    protein_id: str
    genome_id: str
    sequence: str

@dataclass
class DomainHit:
    marker: str
    internal_id: str
    genome_id: str
    protein_id: str
    source_file: str
    target_len: int
    hmm_len: int
    full_evalue: float
    full_score: float
    domain_ievalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    hmm_cov: float
    target_cov: float


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).replace("*", "")
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.replace("*", ""))
        if header is not None:
            yield header, "".join(chunks).replace("*", "")

def discover_fasta_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        return [input_path]
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input not found: {input_path}")

    suffixes = (
        ".faa", ".faa.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz",
        ".pep", ".pep.gz", ".aa", ".aa.gz", ".fna", ".fna.gz",
        ".ffn", ".ffn.gz", ".fsa", ".fsa.gz", ".fas", ".fas.gz"
    )
    files = [p for p in input_path.iterdir() if p.is_file() and str(p.name).endswith(suffixes)]
    return sorted(files, key=lambda p: _natural_key(p.name))


def fasta_stem(path: Path) -> str:
    name = path.name
    for suffix in (
        ".fasta.gz", ".faa.gz", ".fna.gz", ".ffn.gz", ".pep.gz", ".aa.gz",
        ".fsa.gz", ".fas.gz", ".fa.gz", ".fasta", ".faa", ".fna", ".ffn",
        ".pep", ".aa", ".fsa", ".fas", ".fa"
    ):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def detect_fasta_type(path: Path) -> str:
    """Detect whether a FASTA file contains protein or nucleotide sequences."""
    lower = path.name.lower()
    if lower.endswith((".faa", ".faa.gz", ".pep", ".pep.gz", ".aa", ".aa.gz")):
        return "protein"
    if lower.endswith((".fna", ".fna.gz", ".ffn", ".ffn.gz")):
        return "genome"

    total = 0
    nucleotide = 0
    dna_alphabet = set("ACGTUNRYKMSWBDHV")
    for _, seq in iter_fasta(path):
        clean = re.sub(r"[^A-Za-z]", "", seq.upper())
        if not clean:
            continue
        clean = clean[:10000]
        total += len(clean)
        nucleotide += sum(base in dna_alphabet for base in clean)
        if total >= 10000:
            break
    if total == 0:
        raise RuntimeError(f"No sequences found in FASTA file: {path}")
    return "genome" if nucleotide / total >= 0.90 else "protein"


def sanitize_nucleotide_sequence(seq: str) -> str:
    seq = re.sub(r"\s+", "", seq.upper())
    return re.sub(r"[^ACGT]", "N", seq)


def load_genomes_with_pyrodigal(fasta_files: List[Path]) -> Dict[str, ProteinRecord]:
    try:
        import pyrodigal
    except ImportError as exc:
        raise RuntimeError(
            "Genome input requires Pyrodigal. Install it, e.g. pip install pyrodigal "
            "or conda install -c bioconda pyrodigal"
        ) from exc

    records: Dict[str, ProteinRecord] = {}
    counter = 0
    n_contigs = 0
    finder = pyrodigal.GeneFinder(meta=True)

    for fp in fasta_files:
        genome_id = fasta_stem(fp)
        _LOG.info(f"Predicting proteins with Pyrodigal: {fp}")
        file_gene_count = 0
        for header, seq in iter_fasta(fp):
            n_contigs += 1
            contig_id = first_token(header) or f"contig{n_contigs}"
            nt = sanitize_nucleotide_sequence(seq)
            if not nt:
                continue
            genes = finder.find_genes(nt)
            for gene_no, gene in enumerate(genes, start=1):
                protein = sanitize_sequence(gene.translate(include_stop=False))
                if not protein:
                    continue
                counter += 1
                file_gene_count += 1
                internal_id = f"seq{counter:012d}"
                protein_id = f"{contig_id}_{gene_no}"
                strand = "+" if gene.strand == 1 else "-"
                original_header = (
                    f"{protein_id} # {gene.begin} # {gene.end} # {strand} "
                    f"# pyrodigal_meta;source={fp.name}"
                )
                records[internal_id] = ProteinRecord(
                    internal_id=internal_id,
                    source_file=fp.name,
                    original_header=original_header,
                    protein_id=protein_id,
                    genome_id=genome_id,
                    sequence=protein,
                )
        _LOG.info(f"Predicted {file_gene_count:,} protein(s) from genome file: {fp.name}")

    _LOG.info(
        f"Pyrodigal predicted {len(records):,} protein sequence(s) "
        f"from {n_contigs:,} nucleotide sequence(s) in {len(fasta_files):,} file(s)"
    )
    if not records:
        raise RuntimeError("Pyrodigal did not predict any protein sequences from the genome input.")
    return records


def infer_genome_id(header: str, source_file: Path) -> str:
    """Infer GVDB genome_id from a protein FASTA header.

    Typical GVDB headers:
      AbALV.fna|LC506465.1_1 # ...                  -> AbALV
      ERX552244.9.dc.fa|contig_8865_48 # ...        -> ERX552244.9.dc
      GVMAG-M-3300023174-56.fltr JOINED_PROTEIN     -> GVMAG-M-3300023174-56
    """
    first = header.split()[0]
    left = first.split("|", 1)[0]

    for suffix in (".fna", ".fa", ".faa", ".fltr"):
        if left.endswith(suffix):
            left = left[: -len(suffix)]

    # If the header is unexpectedly uninformative, fall back to filename.
    if not left or left in {"-", "unknown", "UNKNOWN"}:
        left = source_file.name
        for suffix in (".faa.gz", ".faa", ".fa.gz", ".fa", ".fasta.gz", ".fasta", ".pep.gz", ".pep"):
            if left.endswith(suffix):
                left = left[: -len(suffix)]
        if left.endswith(".fltr"):
            left = left[:-5]

    return left


def first_token(header: str) -> str:
    return header.split()[0]


def sanitize_sequence(seq: str) -> str:
    seq = seq.upper().replace("*", "")
    # HMMER protein alphabet accepts common AA plus X; convert rare/ambiguous residues conservatively to X.
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZUOJ]", "X", seq)


def load_proteins(fasta_files: List[Path]) -> Dict[str, ProteinRecord]:
    records: Dict[str, ProteinRecord] = {}
    n_empty = 0
    counter = 0

    for fp in fasta_files:
        _LOG.info(f"Reading proteins: {fp}")
        for header, seq in iter_fasta(fp):
            seq = sanitize_sequence(seq)
            if not seq:
                n_empty += 1
                continue
            counter += 1
            internal_id = f"seq{counter:012d}"
            protein_id = first_token(header)
            genome_id = infer_genome_id(header, fp)
            records[internal_id] = ProteinRecord(
                internal_id=internal_id,
                source_file=fp.name,
                original_header=header,
                protein_id=protein_id,
                genome_id=genome_id,
                sequence=seq,
            )

    _LOG.info(f"Loaded {len(records):,} protein sequence(s) from {len(fasta_files):,} file(s)")
    if n_empty:
        _LOG.warning(f"Skipped {n_empty:,} empty protein sequence(s)")
    if not records:
        raise RuntimeError("No protein sequences were loaded. Check input FASTA files.")
    return records

def write_search_fasta(records: Dict[str, ProteinRecord], path: Path) -> None:
    with open(path, "w") as fh:
        for rid in sorted(records.keys()):
            rec = records[rid]
            fh.write(f">{rec.internal_id}\n")
            seq = rec.sequence
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

def infer_hmm_names(hmm_path: Path) -> List[str]:
    names: List[str] = []
    with open(hmm_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split(maxsplit=1)
                if len(parts) == 2:
                    names.append(parts[1].strip())
    return names

def run_hmmsearch(hmm_path: Path, fasta_path: Path, domtblout: Path, stdout_path: Path, evalue: float, threads: int) -> None:
    if not shutil.which("hmmsearch"):
        raise RuntimeError(
            "hmmsearch was not found in PATH. Install HMMER, e.g. conda install -c bioconda hmmer"
        )

    cmd = [
        "hmmsearch",
        "--cpu", str(max(1, threads)),
        "-E", str(evalue),
        "--domE", str(evalue),
        "--domtblout", str(domtblout),
        "-o", str(stdout_path),
        str(hmm_path),
        str(fasta_path),
    ]
    _LOG.info("Running HMMER: " + " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if proc.stdout.strip():
        _LOG.debug(proc.stdout)
    if proc.stderr.strip():
        _LOG.debug(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"hmmsearch failed with exit code {proc.returncode}; see run.log for details")


def parse_domtblout(domtblout: Path, records: Dict[str, ProteinRecord], marker: str) -> List[DomainHit]:
    hits: List[DomainHit] = []
    with open(domtblout, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 23:
                continue

            internal_id = parts[0]
            rec = records.get(internal_id)
            if rec is None:
                continue

            try:
                target_len = int(parts[2])
                hmm_len = int(parts[5])
                full_evalue = float(parts[6])
                full_score = float(parts[7])
                domain_ievalue = float(parts[12])
                domain_score = float(parts[13])
                hmm_from = int(parts[15])
                hmm_to = int(parts[16])
                ali_from = int(parts[17])
                ali_to = int(parts[18])
                env_from = int(parts[19])
                env_to = int(parts[20])
            except ValueError:
                continue

            hmm_cov = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0.0
            target_cov = (ali_to - ali_from + 1) / target_len if target_len > 0 else 0.0

            hits.append(DomainHit(
                marker=marker,
                internal_id=internal_id,
                genome_id=rec.genome_id,
                protein_id=rec.protein_id,
                source_file=rec.source_file,
                target_len=target_len,
                hmm_len=hmm_len,
                full_evalue=full_evalue,
                full_score=full_score,
                domain_ievalue=domain_ievalue,
                domain_score=domain_score,
                hmm_from=hmm_from,
                hmm_to=hmm_to,
                ali_from=ali_from,
                ali_to=ali_to,
                env_from=env_from,
                env_to=env_to,
                hmm_cov=hmm_cov,
                target_cov=target_cov,
            ))
    return hits

def hit_sort_key(h: DomainHit):
    # Best hit means strongest domain score, then better coverage, then smaller i-Evalue.
    return (-h.domain_score, -h.hmm_cov, h.domain_ievalue, h.protein_id)

def select_hits(
    hits: List[DomainHit],
    evalue: float,
    min_hmm_cov: float,
    min_score: Optional[float],
    keep_all: bool,
) -> Tuple[List[DomainHit], List[DomainHit], List[DomainHit]]:
    # First keep the best domain per protein, because domtblout may contain multiple domain rows.
    best_by_protein: Dict[str, DomainHit] = {}
    for h in hits:
        current = best_by_protein.get(h.internal_id)
        if current is None or hit_sort_key(h) < hit_sort_key(current):
            best_by_protein[h.internal_id] = h

    protein_best = list(best_by_protein.values())
    filtered = []
    for h in protein_best:
        if h.domain_ievalue > evalue:
            continue
        if h.hmm_cov < min_hmm_cov:
            continue
        if min_score is not None and h.domain_score < min_score:
            continue
        filtered.append(h)

    filtered.sort(key=lambda h: (h.genome_id, -h.domain_score, -h.hmm_cov, h.domain_ievalue, h.protein_id))

    if keep_all:
        retained = filtered
    else:
        best_by_genome: Dict[str, DomainHit] = {}
        for h in filtered:
            current = best_by_genome.get(h.genome_id)
            if current is None or hit_sort_key(h) < hit_sort_key(current):
                best_by_genome[h.genome_id] = h
        retained = sorted(best_by_genome.values(), key=lambda h: _natural_key(h.genome_id))

    return protein_best, filtered, retained


def write_hits_table(hits: List[DomainHit], path: Path) -> None:
    header = [
        "marker", "genome_id", "protein_id", "source_file", "target_len", "hmm_len",
        "full_evalue", "full_score", "domain_ievalue", "domain_score",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "hmm_cov", "target_cov", "internal_id",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for h in hits:
            row = [
                h.marker, h.genome_id, h.protein_id, h.source_file,
                str(h.target_len), str(h.hmm_len),
                f"{h.full_evalue:.3g}", f"{h.full_score:.1f}",
                f"{h.domain_ievalue:.3g}", f"{h.domain_score:.1f}",
                str(h.hmm_from), str(h.hmm_to), str(h.ali_from), str(h.ali_to),
                str(h.env_from), str(h.env_to),
                f"{h.hmm_cov:.4f}", f"{h.target_cov:.4f}", h.internal_id,
            ]
            fh.write("\t".join(row) + "\n")

def format_evalue(x: float) -> str:
    if x == 0.0:
        return "0"
    return f"{x:.2e}".replace("e-0", "e-").replace("e+0", "e+")

def write_retained_fasta(hits: List[DomainHit], records: Dict[str, ProteinRecord], path: Path) -> None:
    with open(path, "w") as fh:
        for h in hits:
            rec = records[h.internal_id]
            header = (
                f">{h.genome_id}|{h.protein_id}|{h.marker} "
                f"score={h.domain_score:.1f};"
                f"ievalue={format_evalue(h.domain_ievalue)};"
                f"hmmcov={h.hmm_cov:.3f};"
                f"targetcov={h.target_cov:.3f}"
            )
            fh.write(header + "\n")
            seq = rec.sequence
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = _Parser(prog="mine_hallmark.py", add_help=False)
    p.add_argument("-h", "--help", action="help", help="Show this help and exit")
    p.add_argument("input", type=Path)
    p.add_argument("-db", "--db", type=Path, required=True)
    p.add_argument("-o", "--outdir", type=Path, default=Path(f"HallmarkMine_{datetime.now().strftime('%Y%m%d')}"))
    p.add_argument("-e", "--evalue", type=float, default=DEFAULTS["evalue"])
    p.add_argument("-c", "--hmm-coverage", dest="hmm_coverage", type=float, default=DEFAULTS["hmm_coverage"])
    p.add_argument("-t", "--threads", type=int, default=DEFAULTS["threads"])
    p.add_argument("--input-type", choices=("auto", "protein", "genome"), default="auto")
    p.add_argument("--marker", type=str, default=None)
    p.add_argument("--min-score", type=float, default=None)
    p.add_argument("--keep-all", action="store_true")
    p.add_argument("--keep-temp", action="store_true")
    return p.parse_args(argv)

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "tables").mkdir(exist_ok=True)
    (args.outdir / "hallmark").mkdir(exist_ok=True)
    (args.outdir / "logs").mkdir(exist_ok=True)

    setup_logging(args.outdir / "run.log")
    t0 = time.time()

    _LOG.info(
        f"mine_hallmark started | input={args.input} | db={args.db} | "
        f"evalue={args.evalue:g} | hmm_coverage={args.hmm_coverage:.2f} | "
        f"threads={args.threads} | keep_all={args.keep_all} | input_type={args.input_type}"
    )

    if not args.db.is_file():
        _LOG.error(f"HMM file not found: {args.db}")
        return 2

    if args.hmm_coverage < 0 or args.hmm_coverage > 1:
        _LOG.error("--hmm-coverage must be between 0 and 1, for example 0.80")
        return 2

    hmm_names = infer_hmm_names(args.db)
    if args.marker is not None:
        marker = args.marker
    elif len(hmm_names) == 1:
        marker = hmm_names[0]
    elif len(hmm_names) == 0:
        marker = args.db.stem
        _LOG.warning(f"Could not infer HMM NAME from {args.db}; using marker={marker}")
    else:
        _LOG.error(
            f"The HMM file contains multiple profiles: {', '.join(hmm_names[:10])}. "
            "Please provide a single-marker HMM, or set --marker explicitly after extracting one HMM."
        )
        return 2

    _LOG.info(f"Marker name: {marker}")

    try:
        fasta_files = discover_fasta_files(args.input)
        if not fasta_files:
            _LOG.error(f"No FASTA files found in {args.input}")
            return 2

        if args.input_type == "auto":
            by_type: Dict[str, List[Path]] = {"protein": [], "genome": []}
            for fp in fasta_files:
                detected = detect_fasta_type(fp)
                by_type[detected].append(fp)
                _LOG.info(f"Detected input type: {fp.name} -> {detected}")
        else:
            by_type = {"protein": [], "genome": []}
            by_type[args.input_type] = fasta_files

        records: Dict[str, ProteinRecord] = {}
        if by_type["protein"]:
            records.update(load_proteins(by_type["protein"]))
        if by_type["genome"]:
            genome_records = load_genomes_with_pyrodigal(by_type["genome"])
            # Re-number internal IDs if protein and genome inputs are mixed.
            if records:
                offset = len(records)
                for idx, rec in enumerate(genome_records.values(), start=1):
                    new_id = f"seq{offset + idx:012d}"
                    rec.internal_id = new_id
                    records[new_id] = rec
            else:
                records.update(genome_records)

        if not records:
            raise RuntimeError("No protein sequences were available for HMM search.")
    except Exception as exc:
        _LOG.error(str(exc))
        return 2

    work_parent = args.outdir / "tmp" if args.keep_temp else Path(tempfile.mkdtemp(prefix="mine_hallmark_"))
    work_parent.mkdir(parents=True, exist_ok=True)
    search_fasta = work_parent / "search_input.faa"
    domtblout = work_parent / f"{marker}.domtblout"
    hmm_stdout = args.outdir / "logs" / f"{marker}.hmmsearch.out"

    try:
        write_search_fasta(records, search_fasta)
        run_hmmsearch(args.db, search_fasta, domtblout, hmm_stdout, args.evalue, args.threads)

        all_domains = parse_domtblout(domtblout, records, marker)
        protein_best, filtered, retained = select_hits(
            all_domains,
            evalue=args.evalue,
            min_hmm_cov=args.hmm_coverage,
            min_score=args.min_score,
            keep_all=args.keep_all,
        )

        write_hits_table(all_domains, args.outdir / "tables" / f"{marker}.all_domains.tsv")
        write_hits_table(protein_best, args.outdir / "tables" / f"{marker}.best_domain_per_protein.tsv")
        write_hits_table(filtered, args.outdir / "tables" / f"{marker}.filtered.tsv")
        retained_table = args.outdir / "tables" / (f"{marker}.retained_all.tsv" if args.keep_all else f"{marker}.best_per_genome.tsv")
        write_hits_table(retained, retained_table)

        out_faa = args.outdir / "hallmark" / f"{marker}.faa"
        write_retained_fasta(retained, records, out_faa)

        n_genomes = len({r.genome_id for r in records.values()})
        n_hit_genomes = len({h.genome_id for h in filtered})
        _LOG.info("Result Summary")
        _LOG.info(f"  input proteins:             {len(records):,}")
        _LOG.info(f"  input genomes inferred:     {n_genomes:,}")
        _LOG.info(f"  raw domain hits:            {len(all_domains):,}")
        _LOG.info(f"  best domain per protein:    {len(protein_best):,}")
        _LOG.info(f"  filtered protein hits:      {len(filtered):,}")
        _LOG.info(f"  filtered hit genomes:       {n_hit_genomes:,}")
        _LOG.info(f"  retained FASTA sequences:   {len(retained):,}")
        _LOG.output(f"Retained hallmark FASTA -> {out_faa}")
        _LOG.output(f"Retained hit table       -> {retained_table}")

        if not retained:
            notice = args.outdir / f"No_{marker}_hallmark_was_found.txt"
            notice.write_text(
                f"No {marker} hallmark passed the filters.\n"
                f"Filters: domain_iEvalue <= {args.evalue:g}; HMM_coverage >= {args.hmm_coverage:.3f}\n"
                "Check tables/*.all_domains.tsv and logs/*.hmmsearch.out for weak or partial hits.\n"
            )
            _LOG.warning(f"No retained {marker} hits; notice written to {notice}")

    except Exception as exc:
        _LOG.error(str(exc))
        return 1
    finally:
        if not args.keep_temp:
            shutil.rmtree(work_parent, ignore_errors=True)

    _LOG.info(f"Finished in {time.time() - t0:.1f} sec")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
