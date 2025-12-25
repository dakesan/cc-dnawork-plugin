#!/usr/bin/env python3
"""
Test BLAST against human genome (hg38/CHM13)
"""
import json
from pathlib import Path
from Bio.Blast import NCBIWWW, NCBIXML


def test_blast_human(sequence, database, entrez_query=None, description=""):
    """Run BLAST against human genome database."""
    print("=" * 70)
    print(f"BLAST Test: {description}")
    print("=" * 70)
    print(f"Query sequence: {sequence[:60]}...")
    print(f"Length: {len(sequence)} bp")
    print(f"Database: {database}")
    if entrez_query:
        print(f"Entrez query: {entrez_query}")
    print()

    print("Submitting BLAST search to NCBI...")
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database=database,
        sequence=sequence,
        expect=0.001,
        hitlist_size=10,
        entrez_query=entrez_query,
    )

    print("Parsing results...")
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    # Display results
    print("-" * 70)
    print(f"Query: {blast_record.query}")
    print(f"Query length: {blast_record.query_length}")
    print(f"Database: {blast_record.database}")
    print(f"Number of hits: {len(blast_record.alignments)}")
    print()

    if not blast_record.alignments:
        print("No hits found!")
        return None

    # Convert results to JSON format
    results = {
        "description": description,
        "query_sequence": sequence,
        "query_length": blast_record.query_length,
        "database": blast_record.database,
        "entrez_query": entrez_query,
        "num_hits": len(blast_record.alignments),
        "hits": []
    }

    # Display and collect top 10 hits
    print(f"Top {min(10, len(blast_record.alignments))} hits:")
    print("-" * 70)

    for i, alignment in enumerate(blast_record.alignments[:10], 1):
        hsp = alignment.hsps[0]
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

        print(f"\n[{i}] {alignment.title[:100]}")
        print(f"    Accession: {alignment.accession}")
        print(f"    E-value: {hsp.expect:.2e}")
        print(f"    Bit score: {hsp.bits:.1f}")
        print(f"    Identity: {hsp.identities}/{hsp.align_length} ({percent_identity:.1f}%)")

    return results


def main():
    # Test sequence: Human BRCA1 gene fragment (first 100bp of coding sequence)
    # This is a well-known human gene sequence
    brca1_fragment = """
    ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAA
    """.replace("\n", "").replace(" ", "").strip()

    all_results = []

    # Test 1: Search in nt database with Homo sapiens filter
    print("\n" + "=" * 70)
    print("TEST 1: nt database with Homo sapiens filter")
    print("=" * 70 + "\n")

    results1 = test_blast_human(
        sequence=brca1_fragment,
        database="nt",
        entrez_query="Homo sapiens[Organism]",
        description="BRCA1 fragment vs nt database (Homo sapiens only)"
    )
    if results1:
        all_results.append(results1)

    print("\n" + "=" * 70)
    print("Saving combined results...")
    print("=" * 70)

    # Save all results
    if all_results:
        output_file = Path("blast_human_genome_results.json")
        with output_file.open("w") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nResults saved to: {output_file}")

    print("\n" + "=" * 70)
    print("Test completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
