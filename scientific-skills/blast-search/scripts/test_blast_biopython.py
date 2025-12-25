#!/usr/bin/env python3
"""
Simple test script for BioPython BLAST (qblast + NCBIXML parser)
"""
import json
from pathlib import Path
from Bio.Blast import NCBIWWW, NCBIXML


def main():
    # Test sequence (short E. coli 16S rRNA fragment)
    test_sequence = """
    AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCT
    """
    test_sequence = test_sequence.replace("\n", "").replace(" ", "").strip()

    print("=" * 60)
    print("BioPython BLAST Test")
    print("=" * 60)
    print(f"Query sequence: {test_sequence[:50]}...")
    print(f"Length: {len(test_sequence)} bp")
    print()

    # Run BLAST search
    print("Submitting BLAST search to NCBI...")
    print("Program: blastn")
    print("Database: nt (nucleotide collection)")
    print()

    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=test_sequence,
        expect=0.001,
        hitlist_size=10,
    )

    # Parse results
    print("Parsing results...")
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    # Display results
    print("=" * 60)
    print("BLAST Results")
    print("=" * 60)
    print(f"Query: {blast_record.query}")
    print(f"Query length: {blast_record.query_length}")
    print(f"Database: {blast_record.database}")
    print(f"Number of hits: {len(blast_record.alignments)}")
    print()

    if not blast_record.alignments:
        print("No hits found!")
        return

    # Convert results to JSON format
    results = {
        "query": blast_record.query,
        "query_length": blast_record.query_length,
        "database": blast_record.database,
        "num_hits": len(blast_record.alignments),
        "hits": []
    }

    # Display and collect top 10 hits
    print("Top 10 hits:")
    print("-" * 60)

    for i, alignment in enumerate(blast_record.alignments[:10], 1):
        hsp = alignment.hsps[0]  # Best HSP for this alignment
        percent_identity = (hsp.identities / hsp.align_length) * 100

        hit_data = {
            "rank": i,
            "accession": alignment.accession,
            "title": alignment.title,
            "length": alignment.length,
            "e_value": float(hsp.expect),
            "bit_score": float(hsp.bits),
            "score": float(hsp.score),
            "identities": hsp.identities,
            "align_length": hsp.align_length,
            "percent_identity": round(percent_identity, 2),
            "gaps": hsp.gaps,
            "query_start": hsp.query_start,
            "query_end": hsp.query_end,
            "subject_start": hsp.sbjct_start,
            "subject_end": hsp.sbjct_end,
        }
        results["hits"].append(hit_data)

        print(f"\n[{i}] {alignment.title[:80]}")
        print(f"    Accession: {alignment.accession}")
        print(f"    E-value: {hsp.expect:.2e}")
        print(f"    Bit score: {hsp.bits:.1f}")
        print(f"    Identity: {hsp.identities}/{hsp.align_length} ({percent_identity:.1f}%)")
        print(f"    Gaps: {hsp.gaps}")
        print(f"    Query range: {hsp.query_start}-{hsp.query_end}")
        print(f"    Subject range: {hsp.sbjct_start}-{hsp.sbjct_end}")

    # Save results as JSON
    output_file = Path("blast_test_results.json")
    with output_file.open("w") as f:
        json.dump(results, f, indent=2)

    print()
    print("=" * 60)
    print(f"Results saved to: {output_file}")
    print("=" * 60)
    print("Test completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
