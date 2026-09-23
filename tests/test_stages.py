import pytest
from unittest.mock import patch
from pathlib import Path

from swineotype.stages import (
    interpret_resolver,
    stage1_score,
    stage2_resolver_call,
    triplet_at_query_pos,
)
from helpers import REF_1_14, REF_2_12, SUIS

STAGE2_CFG = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}


def ungapped_row(qseqid, pos, triplet, qlen=1000, contig="s1", sstart=1, bitscore=2000):
    """A 13-column resolver row whose subject carries `triplet` at pos-2..pos.

    The query carries the canonical TGG so the row is a plausible alignment;
    only the subject varies.
    """
    qseq = "A" * (pos - 3) + "TGG" + "A" * (qlen - pos)
    sseq = "A" * (pos - 3) + triplet + "A" * (qlen - pos)
    send = sstart + qlen - 1
    return "\t".join([qseqid, contig, "100", str(qlen), str(qlen), "0", str(bitscore),
                      "1", str(qlen), str(sstart), str(send), qseq, sseq])


BASE_CFG = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
            "keep_debug": False, "tmp_dir": "tmp",
            "pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"}}


def cfg(**over):
    c = dict(BASE_CFG)
    c.update(over)
    return c


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
@patch("swineotype.stages.parse_whitelist_headers")
def test_stage1_score(mock_parse_whitelist_headers, mock_make_db_if_needed, mock_run_blast, mock_ensure_tool):
    mock_parse_whitelist_headers.return_value = (
        {"q1": "1", "q2": "14", "q3": "2"},
        {"q1": "wzx", "q2": "wzy", "q3": "wzx"},
        {"1": SUIS, "14": SUIS, "2": SUIS},
    )
    mock_make_db_if_needed.return_value = "db_prefix"
    mock_run_blast.return_value = (
        "q1\ts1\t90\t100\t100\t0\t1000\t1\t100\t1\t100\n"
        "q2\ts2\t95\t100\t100\t0\t2000\t1\t100\t1\t100\n"
        "q3\ts3\t80\t100\t100\t0\t200\t1\t100\t1\t100\n"
    )

    result = stage1_score("assembly.fasta", "whitelist.fasta", 4, Path("run_dir"), cfg())

    assert result["top"] == "14"
    assert result["second"] == "1"
    # 1 and 14 are ONE family, so this is a within-family tie, not competition.
    assert result["family_top"] == "1_vs_14"
    assert result["family_decisive"] is True


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage1_duplicate_copies_do_not_inflate_score(mock_make_db, mock_run_blast, mock_ensure_tool):
    """Regression: a duplicated gene copy used to count twice.

    Two HSPs covering the SAME query interval are two copies in the assembly,
    not two halves of a split gene, so only one may be counted. Left summing,
    the inflated score pushed the wrong serotype to the top of the ranking and
    Stage 2 was then pointed at the wrong resolver pair.
    """
    with patch("swineotype.stages.parse_whitelist_headers",
               return_value=({"dup": "1", "single": "2"}, {"dup": "wzx", "single": "wzx"},
                             {"1": SUIS, "2": SUIS})):
        mock_make_db.return_value = "db_prefix"
        mock_run_blast.return_value = (
            "dup\tc1\t100\t100\t100\t0\t1000\t1\t100\t1\t100\n"     # copy 1
            "dup\tc2\t100\t100\t100\t0\t1000\t1\t100\t1\t100\n"     # copy 2, same query span
            "single\tc3\t100\t100\t100\t0\t1500\t1\t100\t1\t100\n"
        )
        result = stage1_score("a.fasta", "w.fasta", 4, Path("run_dir"), cfg())

    assert result["scores"]["1"] == 1000.0, "duplicate copy must not be counted twice"
    assert result["top"] == "2"


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage1_split_gene_still_counts_both_parts(mock_make_db, mock_run_blast, mock_ensure_tool):
    """Two HSPs covering DIFFERENT query intervals are a gene split across
    contigs and must both count, or fragmented assemblies get under-scored."""
    with patch("swineotype.stages.parse_whitelist_headers",
               return_value=({"split": "9"}, {"split": "wzy"}, {"9": SUIS})):
        mock_make_db.return_value = "db_prefix"
        mock_run_blast.return_value = (
            "split\tc1\t100\t50\t100\t0\t400\t1\t50\t1\t50\n"
            "split\tc2\t100\t50\t100\t0\t400\t51\t100\t1\t50\n"
        )
        result = stage1_score("a.fasta", "w.fasta", 4, Path("run_dir"), cfg())

    assert result["scores"]["9"] == 800.0
    assert result["allele_evidence"]["split"]["split"] is True, \
        "a cross-contig chain must be flagged as such"
    assert result["allele_evidence"]["split"]["contigs"] == ["c1", "c2"]


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_1_vs_14(mock_make_db, mock_run_blast, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_1_14, 492, "TGG")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "1_vs_14")

    assert result["ref_id"] == REF_1_14
    assert result["base"] == "G"
    assert result["triplet"] == "TGG"
    assert result["triplet_status"] == "OK"
    assert interpret_resolver(result, STAGE2_CFG) == "14"


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_2_vs_1_2(mock_make_db, mock_run_blast, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_2_12, 483, "TGT")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "2_vs_1_2")

    assert result["ref_id"] == REF_2_12
    assert result["base"] == "T"
    assert result["triplet"] == "TGT"
    assert interpret_resolver(result, STAGE2_CFG) == "1/2"


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_ignores_references_outside_the_allowed_pair(mock_make_db, mock_run_blast, mock_ensure_tool):
    """The four resolver refs are 98-99.5% identical, so the pair filter is the
    only thing stopping a cross-pair call. Make sure it actually filters."""
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_1_14, 492, "TGG")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "2_vs_1_2")
    assert result is None


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
    status = "OK" if triplet in ("TGG", "TGC", "TGT") else "UNEXPECTED_CODON"
    assert interpret_resolver(ev(REF_2_12, triplet, status), {}) is None


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
