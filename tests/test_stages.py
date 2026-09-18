import pytest
from unittest.mock import patch
from pathlib import Path

from swineotype.stages import (
    base_at_query_pos,
    interpret_resolver,
    stage1_score,
    stage2_resolver_call,
)

# Real reference conventions, so nobody reads these fixtures and infers the
# wrong polarity: G at the diagnostic site -> Trp161 -> serotype 2 or 14.
REF_1_14 = "cps14K|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1"
REF_2_12 = "cps2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2"

SUIS = "Streptococcus suis"

STAGE2_CFG = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}


def ungapped_row(qseqid, pos, base, qlen=1000, contig="s1"):
    """A 13-column resolver BLAST row whose subject carries `base` at `pos`."""
    qseq = "A" * qlen
    sseq = "A" * (pos - 1) + base + "A" * (qlen - pos)
    return "\t".join([qseqid, contig, "100", str(qlen), str(qlen), "0", "2000",
                      "1", str(qlen), "1", str(qlen), qseq, sseq])


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

    result = stage1_score("assembly.fasta", "whitelist.fasta", 4, Path("run_dir"), {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100, "ambig_set": {"1", "14", "2", "1/2"}, "keep_debug": False, "tmp_dir": "tmp"})

    assert result["top"] == "14"
    assert result["second"] == "1"
    assert result["decisive"] is True
    assert result["must_stage2_for_pair"] is True


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
        cfg = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
               "ambig_set": {"1", "14", "2", "1/2"}, "keep_debug": False, "tmp_dir": "tmp"}
        result = stage1_score("a.fasta", "w.fasta", 4, Path("run_dir"), cfg)

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
        cfg = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
               "ambig_set": set(), "keep_debug": False, "tmp_dir": "tmp"}
        result = stage1_score("a.fasta", "w.fasta", 4, Path("run_dir"), cfg)

    assert result["scores"]["9"] == 800.0


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_1_vs_14(mock_make_db, mock_run_blast, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_1_14, 492, "G")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "1_vs_14")

    assert result["ref_id"] == REF_1_14
    assert result["base"] == "G"
    assert interpret_resolver(result, STAGE2_CFG) == "14"


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_2_vs_1_2(mock_make_db, mock_run_blast, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_2_12, 483, "T")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "2_vs_1_2")

    assert result["ref_id"] == REF_2_12
    assert result["base"] == "T"
    assert interpret_resolver(result, STAGE2_CFG) == "1/2"


@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_ignores_references_outside_the_allowed_pair(mock_make_db, mock_run_blast, mock_ensure_tool):
    """The four resolver refs are 98-99.5% identical, so the pair filter is the
    only thing stopping a cross-pair call. Make sure it actually filters."""
    mock_make_db.return_value = "db_prefix"
    mock_run_blast.return_value = ungapped_row(REF_1_14, 492, "G")

    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), STAGE2_CFG, "2_vs_1_2")
    assert result is None


@pytest.mark.parametrize(
    "ref_id, base, expected_serotype",
    [
        (REF_1_14, "G", "14"),
        (REF_1_14, "C", "1"),
        (REF_1_14, "T", "1"),
        (REF_1_14, "A", None),
        (REF_1_14, "-", None),
        (REF_2_12, "G", "2"),
        (REF_2_12, "C", "1/2"),
        (REF_2_12, "T", "1/2"),
        (REF_2_12, "A", None),
        (REF_2_12, "-", None),
    ],
)
def test_interpret_resolver_logic(ref_id, base, expected_serotype):
    assert interpret_resolver({"ref_id": ref_id, "base": base}, {}) == expected_serotype


# --- gap-aware read-out -------------------------------------------------

def test_base_at_query_pos_ungapped():
    qseq = "ACGTACGTAC"
    sseq = "ACGTGCGTAC"
    base, pos, strand = base_at_query_pos(qseq, sseq, 1, 101, 110, 5)
    assert (base, pos, strand) == ("G", 105, "+")


def test_base_at_query_pos_with_upstream_subject_insertion():
    """Regression: sstart + (pos - qstart) assumed an ungapped HSP.

    An extra base in the subject upstream of the site shifted the read-out by
    one -- the characteristic ONT error, and enough to flip serotype 2 to 1/2.
    """
    #            query pos: 1234 567
    qseq = "ACGT-ACGTAC"
    sseq = "ACGTAACGTGC"
    #  query position 9 is the 'G' in sseq (index 9), not index 8
    base, _, _ = base_at_query_pos(qseq, sseq, 1, 1, 11, 9)
    assert base == "G"


def test_base_at_query_pos_with_upstream_subject_deletion():
    qseq = "ACGTAACGTGC"
    sseq = "ACGT-ACGTGC"
    base, _, _ = base_at_query_pos(qseq, sseq, 1, 1, 10, 10)
    assert base == "G"


def test_base_at_query_pos_minus_strand_needs_no_complementing():
    """BLAST already reports sseq on the query strand for a minus-strand hit."""
    qseq = "ACGTACGT"
    sseq = "ACGTGCGT"
    base, pos, strand = base_at_query_pos(qseq, sseq, 1, 200, 193, 5)
    assert strand == "-"
    assert base == "G"
    assert pos == 196


def test_base_at_query_pos_reports_deletion_at_the_site():
    qseq = "ACGTACGT"
    sseq = "ACGT-CGT"
    base, pos, _ = base_at_query_pos(qseq, sseq, 1, 1, 7, 5)
    assert base == "-"
    assert pos is None


def test_base_at_query_pos_off_the_end():
    base, pos, _ = base_at_query_pos("ACGT", "ACGT", 1, 1, 4, 99)
    assert base is None
