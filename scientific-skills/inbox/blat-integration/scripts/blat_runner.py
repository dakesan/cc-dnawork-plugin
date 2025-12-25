#!/usr/bin/env python
"""
BLAT (BLAST-Like Alignment Tool) Integration Module

Provides both UCSC REST API and local command-line BLAT execution
with automatic mode selection based on query size.

Author: cc-dnawork-plugin
License: MIT
"""

import subprocess
import requests
import pandas as pd
import time
import os
from pathlib import Path
from typing import Union, List, Dict, Optional
from Bio import SeqIO


class BlatRunner:
    """
    BLAT sequence search runner with automatic mode selection.

    Supports:
    - UCSC REST API (small-scale, rate-limited)
    - Local BLAT command (large-scale, unlimited)
    - Automatic mode selection
    """

    UCSC_API_BASE = "https://api.genome.ucsc.edu/getData"
    UCSC_BLAT_URL = "https://genome.ucsc.edu/cgi-bin/hgBlat"

    # Rate limits
    API_RETRY_DELAY = 15  # seconds
    API_DAILY_LIMIT = 5000  # hits per day

    def __init__(
        self,
        mode: str = "auto",
        max_api_daily: int = API_DAILY_LIMIT,
        retry_on_limit: bool = True,
        retry_delay: int = API_RETRY_DELAY,
        verbose: bool = False
    ):
        """
        Initialize BlatRunner.

        Args:
            mode: "api", "local", or "auto"
            max_api_daily: Maximum API calls per day
            retry_on_limit: Retry on rate limit
            retry_delay: Delay between retries (seconds)
            verbose: Print debug messages
        """
        self.mode = mode
        self.max_api_daily = max_api_daily
        self.retry_on_limit = retry_on_limit
        self.retry_delay = retry_delay
        self.verbose = verbose
        self.api_call_count = 0

    def search(
        self,
        sequence: Optional[str] = None,
        fasta_file: Optional[str] = None,
        genome: str = "hg38",
        database: Optional[str] = None,
        mode: Optional[str] = None,
        min_identity: float = 0.0,
        max_gap_size: Optional[int] = None,
        min_tile_size: int = 11,
        output_format: str = "psl",
        output_file: Optional[str] = None,
        max_results: Optional[int] = None,
        sort_by: str = "score"
    ) -> Union[Dict, pd.DataFrame, List[Dict]]:
        """
        Search sequence(s) against genomic database.

        Args:
            sequence: Single DNA sequence string
            fasta_file: Path to FASTA file with sequences
            genome: Genome assembly (hg38, mm39, rn7, etc.)
            database: Path to local 2bit database file
            mode: Override self.mode if provided
            min_identity: Minimum % identity (0-100)
            max_gap_size: Maximum gap size in bp
            min_tile_size: Minimum perfect match size
            output_format: Output format (psl, json, bed, dataframe)
            output_file: Save results to file
            max_results: Maximum number of results
            sort_by: Sort results by (score, identity, coverage)

        Returns:
            Results in specified format
        """
        mode = mode or self.mode

        # Read sequences
        if sequence:
            sequences = [sequence]
        elif fasta_file:
            sequences = self._read_fasta(fasta_file)
        else:
            raise ValueError("Provide either 'sequence' or 'fasta_file'")

        # Auto-select mode if needed
        if mode == "auto":
            mode = self._select_mode(sequences)
            if self.verbose:
                print(f"Auto-selected mode: {mode}")

        # Execute search
        if mode == "api":
            results = self._search_api(sequences, genome, min_identity)
        elif mode == "local":
            if not database:
                raise ValueError("'database' required for local mode")
            results = self._search_local(
                sequences,
                database,
                min_tile_size,
                max_gap_size
            )
        else:
            raise ValueError(f"Invalid mode: {mode}")

        # Filter and sort
        if results:
            results = self._filter_results(results, min_identity, max_results)
            results = self._sort_results(results, sort_by)

        # Convert format
        output = self._format_output(results, output_format)

        # Save to file
        if output_file:
            self._save_results(output, output_format, output_file)

        return output

    def _read_fasta(self, fasta_file: str) -> List[str]:
        """Read sequences from FASTA file."""
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
        return sequences

    def _select_mode(self, sequences: List[str]) -> str:
        """Auto-select mode based on query size."""
        # Use local if >1000 sequences or total >100kb
        total_length = sum(len(s) for s in sequences)
        if len(sequences) > 1000 or total_length > 100000:
            return "local"
        return "api"

    def _search_api(self, sequences: List[str], genome: str, min_identity: float) -> List[Dict]:
        """Search using UCSC REST API."""
        results = []

        for i, seq in enumerate(sequences):
            if self.verbose:
                print(f"[API] Searching sequence {i+1}/{len(sequences)}")

            try:
                # Make API request
                params = {
                    "genome": genome,
                    "sequence": seq,
                    "type": "DNA"
                }

                response = requests.get(
                    f"{self.UCSC_API_BASE}/blat",
                    params=params,
                    timeout=30
                )

                if response.status_code == 429:  # Rate limit
                    if self.retry_on_limit:
                        if self.verbose:
                            print(f"Rate limited. Waiting {self.retry_delay}s...")
                        time.sleep(self.retry_delay)
                        response = requests.get(
                            f"{self.UCSC_API_BASE}/blat",
                            params=params,
                            timeout=30
                        )
                    else:
                        raise RuntimeError("API rate limit exceeded")

                if response.status_code == 200:
                    hits = response.json().get("blat", [])
                    results.extend(hits)
                    self.api_call_count += 1
                else:
                    print(f"Warning: API returned {response.status_code}")

            except Exception as e:
                print(f"Error searching sequence {i+1}: {e}")
                continue

        return results

    def _search_local(
        self,
        sequences: List[str],
        database: str,
        min_tile_size: int = 11,
        max_gap_size: Optional[int] = None
    ) -> List[Dict]:
        """Search using local BLAT command."""

        # Check BLAT is installed
        if not self._check_blat_installed():
            raise RuntimeError("BLAT not found. Install with: brew install blat")

        # Write sequences to temp FASTA
        temp_fasta = "/tmp/blat_query.fasta"
        with open(temp_fasta, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq_{i}\n{seq}\n")

        # Run BLAT
        output_psl = "/tmp/blat_output.psl"
        cmd = ["blat", database, temp_fasta, output_psl]

        if self.verbose:
            print(f"Running: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"BLAT failed: {e.stderr.decode()}")

        # Parse PSL results
        results = self._parse_psl(output_psl)

        # Cleanup
        os.remove(temp_fasta)

        return results

    def _parse_psl(self, psl_file: str) -> List[Dict]:
        """Parse PSL format output from BLAT."""
        results = []

        with open(psl_file) as f:
            # Skip header lines
            for _ in range(5):
                f.readline()

            for line in f:
                if not line.strip():
                    continue

                fields = line.strip().split()
                if len(fields) < 21:
                    continue

                # PSL format fields
                result = {
                    "match": int(fields[0]),
                    "mismatch": int(fields[1]),
                    "rep_match": int(fields[2]),
                    "n_count": int(fields[3]),
                    "q_num_insert": int(fields[4]),
                    "q_base_insert": int(fields[5]),
                    "t_num_insert": int(fields[6]),
                    "t_base_insert": int(fields[7]),
                    "strand": fields[8],
                    "q_name": fields[9],
                    "q_size": int(fields[10]),
                    "q_start": int(fields[11]),
                    "q_end": int(fields[12]),
                    "t_name": fields[13],
                    "t_size": int(fields[14]),
                    "t_start": int(fields[15]),
                    "t_end": int(fields[16]),
                }

                # Calculate metrics
                total_matches = result["match"] + result["mismatch"]
                if total_matches > 0:
                    result["identity"] = 100 * result["match"] / total_matches
                else:
                    result["identity"] = 0

                query_insert = result["q_base_insert"]
                target_insert = result["t_base_insert"]
                gap_penalty = query_insert + target_insert

                # Calculate score (UCSC metric)
                result["score"] = result["match"] - result["mismatch"] - gap_penalty

                results.append(result)

        return results

    def _filter_results(
        self,
        results: List[Dict],
        min_identity: float = 0.0,
        max_results: Optional[int] = None
    ) -> List[Dict]:
        """Filter results by criteria."""

        # Filter by identity
        filtered = [r for r in results if r.get("identity", 0) >= min_identity]

        # Limit results
        if max_results:
            filtered = filtered[:max_results]

        return filtered

    def _sort_results(self, results: List[Dict], sort_by: str = "score") -> List[Dict]:
        """Sort results by specified field."""

        valid_sort_keys = ["score", "identity", "match", "t_name"]
        if sort_by not in valid_sort_keys:
            raise ValueError(f"sort_by must be one of {valid_sort_keys}")

        return sorted(results, key=lambda x: x.get(sort_by, 0), reverse=True)

    def _format_output(
        self,
        results: List[Dict],
        output_format: str = "psl"
    ) -> Union[pd.DataFrame, Dict, List[Dict]]:
        """Convert results to requested format."""

        if output_format == "dataframe":
            return pd.DataFrame(results)
        elif output_format == "json":
            return {
                "count": len(results),
                "hits": results
            }
        elif output_format == "bed":
            bed_data = []
            for r in results:
                bed_data.append({
                    "chrom": r.get("t_name"),
                    "chromStart": r.get("t_start"),
                    "chromEnd": r.get("t_end"),
                    "name": r.get("q_name"),
                    "score": r.get("score"),
                    "strand": r.get("strand")
                })
            return pd.DataFrame(bed_data)
        else:  # psl
            return results

    def _save_results(
        self,
        results: Union[pd.DataFrame, List[Dict], Dict],
        output_format: str,
        output_file: str
    ):
        """Save results to file."""

        if isinstance(results, pd.DataFrame):
            if output_file.endswith(".csv"):
                results.to_csv(output_file, index=False)
            elif output_file.endswith(".xlsx"):
                results.to_excel(output_file, index=False)
            else:
                results.to_csv(output_file, index=False, sep="\t")
        else:
            import json
            with open(output_file, "w") as f:
                json.dump(results, f, indent=2)

        if self.verbose:
            print(f"Results saved to {output_file}")

    def _check_blat_installed(self) -> bool:
        """Check if BLAT is installed."""
        try:
            subprocess.run(["blat"], capture_output=True)
            return True
        except FileNotFoundError:
            return False

    def check_installation(self) -> Dict[str, str]:
        """Check local BLAT installation."""
        result = {}

        try:
            output = subprocess.run(
                ["blat"],
                capture_output=True,
                text=True
            )
            result["blat_installed"] = True
            result["blat_version"] = output.stderr.split("\n")[0]
        except FileNotFoundError:
            result["blat_installed"] = False
            result["blat_version"] = "Not found"

        # Try to find blat in PATH
        try:
            path = subprocess.run(
                ["which", "blat"],
                capture_output=True,
                text=True
            ).stdout.strip()
            result["blat_path"] = path
        except:
            result["blat_path"] = "Not found"

        return result

    def test_api_connection(self) -> Dict[str, str]:
        """Test UCSC API connectivity."""
        try:
            response = requests.get(
                f"{self.UCSC_API_BASE}/blat",
                params={"genome": "hg38", "sequence": "ATGC", "type": "DNA"},
                timeout=10
            )

            if response.status_code == 200:
                return {
                    "status": "OK",
                    "api_available": True
                }
            else:
                return {
                    "status": f"Error {response.status_code}",
                    "api_available": False
                }
        except Exception as e:
            return {
                "status": f"Connection error: {e}",
                "api_available": False
            }


if __name__ == "__main__":
    # Test example
    runner = BlatRunner(verbose=True)

    # Test single sequence (API mode)
    print("Testing API mode...")
    sequence = "ATGCGTACGATCGATCGATCGATCGATCG"
    results = runner.search(
        sequence=sequence,
        genome="hg38",
        mode="api",
        output_format="dataframe"
    )

    if isinstance(results, pd.DataFrame):
        print(f"Found {len(results)} hits")
        print(results.head())
    else:
        print("No results")

    # Check installation
    print("\nChecking installation...")
    print(runner.check_installation())
    print(runner.test_api_connection())
