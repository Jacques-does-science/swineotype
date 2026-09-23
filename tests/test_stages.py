import pytest

from helpers import (REF_1_14, REF_2_12, STAGE2_CFG, blast_row, resolver_row, run_stage1,
                     run_stage2, stage1_cfg)
from swineotype.stages import interpret_resolver, triplet_at_query_pos


def test_stage1_score():
    result = run_stage1({"q1": ("1", "wzx"), "q2": ("14", "wzy"), "q3": ("2", "wzx")},
                        [blast_row("q1", 1000, pid=90, contig="s1"),
                         blast_row("q2", 2000, pid=95, contig="s2"),
                         blast_row("q3", 200, pid=80, contig="s3")],
                        stage1_cfg(require_wzy=0))

    assert result["top"] == "14"
    assert result["second"] == "1"
    # 1 and 14 are ONE family, so this is a within-family tie, not competition.
    assert result["family_top"] == "1_vs_14"
    assert result["family_decisive"] is True


def test_stage1_duplicate_copies_do_not_inflate_score():
    """Regression: a duplicated gene copy used to count twice.

    Two HSPs covering the SAME query interval are two copies in the assembly,
    not two halves of a split gene, so only one may be counted. Left summing,
    the inflated score pushed the wrong serotype to the top of the ranking and
    Stage 2 was then pointed at the wrong resolver pair.
    """
    result = run_stage1({"dup": ("1", "wzx"), "single": ("2", "wzx")},
                        [blast_row("dup", 1000, contig="c1"),      # copy 1
                         blast_row("dup", 1000, contig="c2"),      # copy 2, same query span
                         blast_row("single", 1500, contig="c3")],
                        stage1_cfg(require_wzy=0))

    assert result["scores"]["1"] == 1000.0, "duplicate copy must not be counted twice"
    assert result["top"] == "2"


def test_stage1_split_gene_still_counts_both_parts():
    """Two HSPs covering DIFFERENT query intervals are a gene split across
    contigs and must both count, or fragmented assemblies get under-scored."""
    result = run_stage1({"split": ("9", "wzy")},
                        ["split\tc1\t100\t50\t100\t0\t400\t1\t50\t1\t50",
                         "split\tc2\t100\t50\t100\t0\t400\t51\t100\t1\t50"])

    assert result["scores"]["9"] == 800.0
    assert result["allele_evidence"]["split"]["split"] is True, \
        "a cross-contig chain must be flagged as such"
    assert result["allele_evidence"]["split"]["contigs"] == ["c1", "c2"]


def test_stage2_resolver_call_1_vs_14():
    result = run_stage2([resolver_row(REF_1_14, 492, "TGG")], "1_vs_14")

    assert result["ref_id"] == REF_1_14
    assert result["base"] == "G"
    assert result["triplet"] == "TGG"
    assert result["triplet_status"] == "OK"
    assert interpret_resolver(result, STAGE2_CFG) == "14"


def test_stage2_resolver_call_2_vs_1_2():
    result = run_stage2([resolver_row(REF_2_12, 483, "TGT")], "2_vs_1_2")

    assert result["ref_id"] == REF_2_12
    assert result["base"] == "T"
    assert result["triplet"] == "TGT"
    assert interpret_resolver(result, STAGE2_CFG) == "1/2"


def test_stage2_ignores_references_outside_the_allowed_pair():
    """The four resolver refs are 98-99.5% identical, so the pair filter is the
    only thing stopping a cross-pair call. Make sure it actually filters."""
    assert run_stage2([resolver_row(REF_1_14, 492, "TGG")], "2_vs_1_2") is None


# --- codon interpretation ----------------------------------------------

def ev(ref_id, triplet, status="OK", **over):
    e = {"ref_id": ref_id, "triplet": triplet, "triplet_status": status,
         "base": triplet[-1], "coding_status": "INTACT"}
    e.update(over)
    return e


@pytest.mark.parametrize("ref_id, triplet, expected_serotype", [
    (REF_1_14, "TGG", "14"),
    (REF_1_14, "TGC", "1"),
    (REF_1_14, "TGT", "1"),
    (REF_2_12, "TGG", "2"),
    (REF_2_12, "TGC", "1/2"),
    (REF_2_12, "TGT", "1/2"),
])
def test_interpret_resolver_accepts_the_three_documented_codons(ref_id, triplet, expected_serotype):
    assert interpret_resolver(ev(ref_id, triplet), {}) == expected_serotype


@pytest.mark.parametrize("triplet", ["AGG", "CGG", "GGG", "TAG", "TGA", "TAA", "TTG"])
def test_interpret_resolver_rejects_undocumented_codons(triplet):
    """Regression: only the wobble base was read, so ANY triplet ending in G
    was accepted as Trp. An AGG (Arg) subject codon was called serotype 2."""
    assert interpret_resolver(ev(REF_2_12, triplet, "UNEXPECTED_CODON"), {}) is None


@pytest.mark.parametrize("triplet, status", [
    ("TNG", "AMBIGUOUS"), ("NGG", "AMBIGUOUS"),
    ("T-G", "DELETED"), ("TG-", "DELETED"),
    ("?GG", "INCOMPLETE"), ("TG?", "INCOMPLETE"),
    ("TGG", "NON_CONTIGUOUS"),
])
def test_interpret_resolver_rejects_uninterpretable_triplets(triplet, status):
    assert interpret_resolver(ev(REF_2_12, triplet, status), {}) is None


def test_interpret_resolver_withholds_on_detected_coding_disruption():
    """A recovered coordinate is not proof the gene around it still reads."""
    assert interpret_resolver(ev(REF_2_12, "TGG", coding_status="DISRUPTED"), {}) is None


def test_interpret_resolver_allows_a_call_when_integrity_is_unassessed():
    """Most alignments are partial; UNASSESSED is reported, not penalised."""
    assert interpret_resolver(ev(REF_2_12, "TGG", coding_status="UNASSESSED"), {}) == "2"


def test_interpret_resolver_withholds_on_conflict():
    assert interpret_resolver(ev(REF_2_12, "TGG", conflict=True), {}) is None


def test_interpret_resolver_of_nothing():
    assert interpret_resolver(None, {}) is None


# --- gap-aware read-out -------------------------------------------------
#
# These pinned the retired single-base reader; they now run against
# triplet_at_query_pos(), which is what production calls.

def test_readout_ungapped():
    t = triplet_at_query_pos("ACGTACGTAC", "ACGTGCGTAC", 1, 101, 110, 5)
    assert (t["base"], t["contig_pos"], t["strand"]) == ("G", 105, "+")


def test_readout_with_upstream_subject_insertion():
    """Regression: sstart + (pos - qstart) assumed an ungapped HSP.

    An extra base in the subject upstream of the site shifted the read-out by
    one -- the characteristic ONT error, and enough to flip serotype 2 to 1/2.
    """
    #            query pos: 1234 567
    qseq = "ACGT-ACGTAC"
    sseq = "ACGTAACGTGC"
    #  query position 9 is the 'G' in sseq (index 9), not index 8
    assert triplet_at_query_pos(qseq, sseq, 1, 1, 11, 9)["base"] == "G"


def test_readout_with_upstream_subject_deletion():
    assert triplet_at_query_pos("ACGTAACGTGC", "ACGT-ACGTGC", 1, 1, 10, 10)["base"] == "G"


def test_readout_minus_strand_needs_no_complementing():
    """BLAST already reports sseq on the query strand for a minus-strand hit."""
    t = triplet_at_query_pos("ACGTACGT", "ACGTGCGT", 1, 200, 193, 5)
    assert t["strand"] == "-"
    assert t["base"] == "G"
    assert t["contig_pos"] == 196


def test_readout_reports_deletion_at_the_site():
    t = triplet_at_query_pos("ACGTACGT", "ACGT-CGT", 1, 1, 7, 5)
    assert t["base"] == "-"
    assert t["contig_pos"] is None
    assert t["triplet_status"] == "DELETED"


def test_readout_off_the_end():
    t = triplet_at_query_pos("ACGT", "ACGT", 1, 1, 4, 99)
    assert t["contig_pos"] is None
    assert t["triplet_status"] == "INCOMPLETE"


# --- whole-triplet read-out ---------------------------------------------

def test_triplet_ungapped_plus_strand():
    qseq = "AATGGCC"
    sseq = "AATGGCC"
    t = triplet_at_query_pos(qseq, sseq, 1, 101, 107, 5)
    assert t["triplet"] == "TGG"
    assert t["triplet_status"] == "OK"
    assert t["positions"] == [103, 104, 105]
    assert t["contig_pos"] == 105
    assert t["strand"] == "+"


def test_triplet_minus_strand_is_not_complemented_twice():
    """sseq already arrives on the query strand; only coordinates run back."""
    qseq = "AATGGCC"
    sseq = "AATGGCC"
    t = triplet_at_query_pos(qseq, sseq, 1, 207, 201, 5)
    assert t["strand"] == "-"
    assert t["triplet"] == "TGG"
    assert t["triplet_status"] == "OK"
    assert t["positions"] == [205, 204, 203]


def test_triplet_agg_is_not_reported_as_tgg():
    """The whole point: reading only the wobble base made AGG look like TGG."""
    t = triplet_at_query_pos("AATGGCC", "AAAGGCC", 1, 1, 7, 5)
    assert t["base"] == "G", "the wobble base alone is still G"
    assert t["triplet"] == "AGG"
    assert t["triplet_status"] == "UNEXPECTED_CODON"


def test_triplet_with_deleted_diagnostic_base():
    t = triplet_at_query_pos("AATGGCC", "AATG-CC", 1, 1, 6, 5)
    assert t["triplet"] == "TG-"
    assert t["triplet_status"] == "DELETED"


def test_triplet_with_deletion_elsewhere_in_the_codon():
    t = triplet_at_query_pos("AATGGCC", "AA-GGCC", 1, 1, 6, 5)
    assert t["triplet"] == "-GG"
    assert t["triplet_status"] == "DELETED"


def test_triplet_with_ambiguous_base():
    t = triplet_at_query_pos("AATGGCC", "AATNGCC", 1, 1, 7, 5)
    assert t["triplet"] == "TNG"
    assert t["triplet_status"] == "AMBIGUOUS"


def test_triplet_partially_off_the_alignment():
    """The HSP reaches the diagnostic base but not the start of its codon."""
    #  alignment starts at query position 4, diagnostic position 5
    t = triplet_at_query_pos("GGCC", "GGCC", 4, 1, 4, 5)
    assert t["triplet_status"] == "INCOMPLETE"
    assert t["triplet"] == "?GG", "the unreached position is marked, not guessed"


def test_triplet_broken_by_a_subject_insertion_inside_the_codon():
    """An inserted subject base between two codon positions means those three
    bases are not a codon of the subject's own gene."""
    #  query:   A A T - G G   (query gap = subject insertion)
    qseq = "AAT-GCC"
    sseq = "AATAGCC"
    t = triplet_at_query_pos(qseq, sseq, 1, 1, 7, 5)
    assert t["triplet"] == "TGC"
    assert t["triplet_status"] == "NON_CONTIGUOUS"
    assert t["positions"] == [3, 5, 6]


def test_triplet_upstream_indel_does_not_shift_the_codon():
    """An indel before the codon shifts coordinates but not the read-out."""
    #  subject has one extra base at query position 2/3 boundary
    qseq = "AA-TGGCC"
    sseq = "AAATGGCC"
    t = triplet_at_query_pos(qseq, sseq, 1, 1, 8, 5)
    assert t["triplet"] == "TGG"
    assert t["triplet_status"] == "OK"
    assert t["positions"] == [4, 5, 6]
