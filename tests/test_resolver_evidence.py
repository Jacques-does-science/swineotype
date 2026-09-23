"""Distinct physical copies of cpsK, and what the tool does when they disagree.

Stage 2 used to keep only the highest-bit-score eligible HSP. If the assembly
carried two copies of the locus with different diagnostic codons, one was
silently chosen and reported with full confidence; which one depended on bit
score, and therefore on which reference happened to align best.

Two things have to be separated to fix that:
  * several references hitting ONE copy (the normal case -- the four resolver
    references are 98-99.5% identical, so two or three of them find every
    copy). Different query names are not evidence of multiple copies.
  * one reference hitting SEVERAL copies at different subject coordinates.
    That is a real disagreement, and the exact label must be withheld.
"""
from pathlib import Path
from unittest.mock import patch

import pytest

from helpers import REF_1_14, REF_2_12, REF_2_12_ALT
from swineotype.stages import (
    assess_coding_integrity,
    collapse_to_loci,
    resolver_status,
    stage2_resolver_call,
)

CFG = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}


def resolver_row(qseqid, pos, triplet, contig="c1", sstart=1, bitscore=2000,
                 qlen=1000, pident=100):
    qseq = "A" * (pos - 3) + "TGG" + "A" * (qlen - pos)
    sseq = "A" * (pos - 3) + triplet + "A" * (qlen - pos)
    send = sstart + qlen - 1
    return "\t".join([qseqid, contig, str(pident), str(qlen), str(qlen), "0",
                      str(bitscore), "1", str(qlen), str(sstart), str(send), qseq, sseq])


def run(rows, allowed_pair="2_vs_1_2"):
    with patch("swineotype.stages.ensure_tool"), \
         patch("swineotype.stages.make_db_if_needed", return_value="db"), \
         patch("swineotype.stages.run_blast", return_value="\n".join(rows)):
        return stage2_resolver_call("a.fa", "r.fa", 1, Path("run"), CFG, allowed_pair)


# --- collapsing redundant alignments -----------------------------------

def test_two_references_hitting_one_contig_are_one_locus():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="c1", bitscore=2000),
              resolver_row(REF_2_12_ALT, 483, "TGG", contig="c1", bitscore=1900)])
    assert ev["n_loci"] == 1
    assert ev["conflict"] is False
    assert sorted(ev["loci"][0]["ref_ids"]) == sorted([REF_2_12, REF_2_12_ALT])
    assert ev["loci"][0]["n_alignments"] == 2


def test_overlapping_subject_spans_on_one_contig_are_one_locus():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="c1", sstart=1),
              resolver_row(REF_2_12_ALT, 483, "TGG", contig="c1", sstart=400)])
    assert ev["n_loci"] == 1


def test_distinct_subject_spans_on_one_contig_are_distinct_loci():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="c1", sstart=1),
              resolver_row(REF_2_12, 483, "TGG", contig="c1", sstart=50_000)])
    assert ev["n_loci"] == 2
    assert ev["conflict"] is False, "they agree, so there is nothing to withhold"


def test_a_locus_keeps_its_best_alignment_as_the_representative():
    ev = run([resolver_row(REF_2_12, 483, "TGG", bitscore=1500),
              resolver_row(REF_2_12_ALT, 483, "TGG", bitscore=1900)])
    assert ev["ref_id"] == REF_2_12_ALT
    assert ev["bitscore"] == 1900


# --- conflict ----------------------------------------------------------

def test_two_copies_with_different_codons_conflict():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
              resolver_row(REF_2_12, 483, "TGT", contig="cB")])
    assert ev["conflict"] is True
    assert ev["conflicting_readings"] == ["1/2", "2"]
    assert ev["n_loci"] == 2
    assert resolver_status(ev) == "CONFLICTING_COPIES"


def test_conflict_is_not_resolved_by_bit_score():
    """A clear bit-score winner must not settle a real disagreement."""
    high = run([resolver_row(REF_2_12, 483, "TGG", contig="cA", bitscore=9999),
                resolver_row(REF_2_12, 483, "TGT", contig="cB", bitscore=100)])
    assert high["conflict"] is True
    assert high["conflicting_readings"] == ["1/2", "2"]


def test_conflict_is_not_resolved_by_fasta_order():
    forward = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
                   resolver_row(REF_2_12, 483, "TGT", contig="cB")])
    reverse = run([resolver_row(REF_2_12, 483, "TGT", contig="cB"),
                   resolver_row(REF_2_12, 483, "TGG", contig="cA")])
    assert forward["conflict"] is reverse["conflict"] is True
    assert forward["conflicting_readings"] == reverse["conflicting_readings"]


def test_both_pieces_of_evidence_are_preserved():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
              resolver_row(REF_2_12, 483, "TGT", contig="cB")])
    by_contig = {locus["contig"]: locus for locus in ev["loci"]}
    assert by_contig["cA"]["triplet"] == "TGG"
    assert by_contig["cB"]["triplet"] == "TGT"
    assert by_contig["cA"]["implied_serotype"] == "2"
    assert by_contig["cB"]["implied_serotype"] == "1/2"


def test_an_agreeing_second_copy_is_not_a_conflict():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
              resolver_row(REF_2_12, 483, "TGG", contig="cB")])
    assert ev["conflict"] is False
    assert resolver_status(ev) == "OK"


def test_one_valid_and_one_broken_copy_conflict():
    """A second copy whose codon is AGG is not agreement, even though it
    implies no serotype of its own."""
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
              resolver_row(REF_2_12, 483, "AGG", contig="cB")])
    assert ev["conflict"] is True
    assert ev["conflicting_readings"] == ["2", "AGG"], \
        "the undocumented codon is reported as itself, not dropped"


@pytest.mark.parametrize("codon_a, codon_b, serotype", [
    ("TGT", "TGC", "1/2"),
    ("TGC", "TGT", "1/2"),
])
def test_synonymous_codons_in_two_copies_agree(codon_a, codon_b, serotype):
    """Regression: codon SPELLINGS were compared, so a TGT copy and a TGC copy
    were reported as conflicting -- with `conflicting_serotypes: ['1/2']`, a
    single serotype listed as disagreeing with itself -- and the call was
    withheld. Both encode Cys161; they imply the same serotype."""
    from swineotype.stages import interpret_resolver
    ev = run([resolver_row(REF_2_12, 483, codon_a, contig="cA"),
              resolver_row(REF_2_12, 483, codon_b, contig="cB")])

    assert ev["n_loci"] == 2
    assert ev["conflict"] is False
    assert ev["conflicting_readings"] == []
    assert resolver_status(ev) == "OK"
    assert interpret_resolver(ev, {}) == serotype


def test_no_qualifying_alignment_returns_nothing():
    assert run([]) is None
    assert run([resolver_row(REF_2_12, 483, "TGG", pident=50)]) is None
    assert run([resolver_row(REF_1_14, 492, "TGG")], allowed_pair="2_vs_1_2") is None
    assert resolver_status(None) == "NO_HSP_OR_LOW_QUAL"


def test_collapse_to_loci_is_stable_for_an_empty_input():
    assert collapse_to_loci([]) == []


# --- coding integrity --------------------------------------------------

def test_full_length_clean_alignment_is_intact():
    q = "ATG" + "AAA" * 10 + "TAA"
    r = assess_coding_integrity(q, q, 1, len(q), len(q))
    assert r["coding_status"] == "INTACT"


def test_a_frameshift_is_detected():
    """One unpaired subject insertion puts everything downstream out of frame."""
    q = "ATG" + "AAA" * 10 + "-" + "TAA"
    s = "ATG" + "AAA" * 10 + "C" + "TAA"
    r = assess_coding_integrity(q, s, 1, 36, 36)
    assert r["coding_status"] == "DISRUPTED"
    assert "frameshift" in r["coding_detail"]


def test_a_balanced_three_base_indel_is_not_a_frameshift():
    q = "ATG" + "AAA" * 10 + "---" + "TAA"
    s = "ATG" + "AAA" * 10 + "CCC" + "TAA"
    r = assess_coding_integrity(q, s, 1, 36, 36)
    assert r["coding_status"] == "INTACT"


def test_a_premature_stop_is_detected():
    q = "ATG" + "AAA" * 10 + "TAA"
    s = "ATG" + "AAA" * 4 + "TAA" + "AAA" * 5 + "TAA"
    r = assess_coding_integrity(q, s[:len(q)], 1, len(q), len(q))
    assert r["coding_status"] == "DISRUPTED"
    assert "premature_stop" in r["coding_detail"]


def test_a_partial_alignment_is_unassessed_not_intact():
    """Recovering a coordinate across an indel is not proof of an intact gene.
    Saying so is the whole point of this state."""
    q = "AAA" * 10
    r = assess_coding_integrity(q, q, 100, 129, 1005)
    assert r["coding_status"] == "UNASSESSED"
    assert "alignment_covers_query_100-129_of_1005" in r["coding_detail"]


def test_a_partial_alignment_can_still_be_positively_disrupted():
    q = "AAA" * 5 + "-" + "AAA" * 4
    s = "AAA" * 5 + "C" + "AAA" * 4
    r = assess_coding_integrity(q, s, 100, 127, 1005)
    assert r["coding_status"] == "DISRUPTED"


@pytest.mark.parametrize("status, expected", [
    ("OK", "OK"),
])
def test_resolver_status_of_a_clean_locus(status, expected):
    ev = {"triplet_status": status, "coding_status": "INTACT", "conflict": False}
    assert resolver_status(ev) == expected


def test_resolver_status_precedence():
    """Conflict outranks everything: it is the state a user must not miss."""
    assert resolver_status({"conflict": True, "triplet_status": "OK",
                            "coding_status": "INTACT"}) == "CONFLICTING_COPIES"
    assert resolver_status({"conflict": False, "triplet_status": "AMBIGUOUS",
                            "coding_status": "INTACT"}) == "INVALID_TRIPLET"
    assert resolver_status({"conflict": False, "triplet_status": "OK",
                            "coding_status": "DISRUPTED"}) == "CODING_DISRUPTED"


def test_coding_disruption_withholds_the_call():
    from swineotype.stages import interpret_resolver
    ev = {"ref_id": REF_2_12, "triplet": "TGG", "triplet_status": "OK",
          "coding_status": "DISRUPTED", "conflict": False}
    assert interpret_resolver(ev, {}) is None


# --- a locus whose codon could not be read is not a disagreement --------

def truncated_row(qseqid, pos, contig, sstart=1, bitscore=800, qlen=1000, alen=400):
    """An HSP that reaches the diagnostic base but stops before its codon.

    This is what one physical gene looks like when the assembly breaks inside
    it and the two contigs overlap across the diagnostic position.
    """
    qstart = pos - 1                      # starts between pos-2 and pos
    qseq = "TG" + "A" * (alen - 2)
    sseq = qseq
    return "\t".join([qseqid, contig, "100", str(alen), str(qlen), "0", str(bitscore),
                      str(qstart), str(qstart + alen - 1), str(sstart),
                      str(sstart + alen - 1), qseq, sseq])


def test_a_truncated_second_locus_does_not_create_a_false_conflict():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA", bitscore=2000),
              truncated_row(REF_2_12, 483, contig="cB")])

    assert ev["n_loci"] == 2, "the second locus is still recorded"
    statuses = {locus["contig"]: locus["triplet_status"] for locus in ev["loci"]}
    assert statuses["cB"] == "INCOMPLETE"
    assert ev["conflict"] is False, \
        "a codon that could not be read says nothing about whether it agrees"
    assert resolver_status(ev) == "OK"


def test_the_primary_locus_is_still_the_best_scoring_one():
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA", bitscore=2000),
              truncated_row(REF_2_12, 483, contig="cB", bitscore=800)])
    assert ev["contig"] == "cA"
    assert ev["triplet"] == "TGG"


def test_two_fully_read_codons_that_differ_still_conflict():
    """The guard above must not weaken the case it exists for."""
    ev = run([resolver_row(REF_2_12, 483, "TGG", contig="cA"),
              resolver_row(REF_2_12, 483, "TGT", contig="cB")])
    assert ev["conflict"] is True


def test_the_readable_locus_becomes_primary_even_at_a_lower_bit_score():
    """Regression: the longest alignment on a fragmented assembly is often the
    one that stops just short of the diagnostic codon. Letting it become the
    primary turned a perfectly readable call into a no-call, and which locus
    won depended on the order BLAST emitted its rows in."""
    ev = run([truncated_row(REF_2_12, 483, contig="cLong", bitscore=5000),
              resolver_row(REF_2_12, 483, "TGG", contig="cShort", bitscore=900)])

    assert ev["contig"] == "cShort"
    assert ev["triplet"] == "TGG"
    assert resolver_status(ev) == "OK"


def test_locus_order_does_not_depend_on_blast_row_order():
    rows = [resolver_row(REF_2_12, 483, "TGG", contig="cA", bitscore=1000),
            resolver_row(REF_2_12, 483, "TGG", contig="cB", bitscore=1000)]
    forward = run(rows)
    reverse = run(list(reversed(rows)))
    assert [l["contig"] for l in forward["loci"]] == [l["contig"] for l in reverse["loci"]]
    assert forward["contig"] == reverse["contig"]
