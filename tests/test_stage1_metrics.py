"""Coverage, identity and score must describe the SAME evidence.

The previous aggregation computed each of the three from a different set of
HSPs: coverage was unioned over every HSP of the allele, identity was
length-weighted over every HSP including overlapping copies, and only the
score used a greedy non-overlapping subset. A full-length 100%-identity match
therefore passed on its own, but failed once five overlapping 80%-identity
copies were added -- they dragged the reported identity below the threshold
while contributing nothing at all to the score.

Aggregation policy now ("best consistent copy, then an explicit split chain"):

  1. HSPs are partitioned into physical copies: same contig, same strand, no
     query overlap, collinear in query and subject.
  2. The allele's evidence is the single highest-scoring copy. Coverage,
     identity and score all come from that copy.
  3. Only if that copy alone falls short of min_cov is a chain across copies
     attempted, adding HSPs on disjoint query intervals best-score-first. The
     result is flagged `split` and carries every contig it used, because
     fragments from unrelated locations are not by themselves one intact gene.
"""
import pytest

from swineotype.stages import allele_evidence, group_into_copies

FULL = 1000


def hsp(qstart, qend, sseqid="c1", sstart=None, pident=100.0, bitscore=None, strand="+"):
    length = qend - qstart + 1
    sstart = qstart if sstart is None else sstart
    send = sstart + length - 1 if strand == "+" else sstart - length + 1
    return {"qstart": qstart, "qend": qend, "sseqid": sseqid, "sstart": sstart,
            "send": send, "strand": strand, "pident": pident, "length": length,
            "qlen": FULL, "bitscore": bitscore if bitscore is not None else length * 1.8}


# --- copy partitioning -------------------------------------------------

def test_overlapping_hsps_on_one_contig_are_separate_copies():
    copies = group_into_copies([hsp(1, 1000), hsp(1, 1000, sstart=5000)])
    assert len(copies) == 2


def test_disjoint_collinear_hsps_on_one_contig_are_one_copy():
    copies = group_into_copies([hsp(1, 500), hsp(501, 1000, sstart=501)])
    assert len(copies) == 1


def test_hsps_on_different_contigs_are_different_copies():
    copies = group_into_copies([hsp(1, 500), hsp(501, 1000, sseqid="c2", sstart=1)])
    assert len(copies) == 2


def test_non_collinear_hsps_are_different_copies():
    """Query advances but subject goes backwards: a rearrangement, not a gene."""
    copies = group_into_copies([hsp(1, 500, sstart=5000), hsp(501, 1000, sstart=1)])
    assert len(copies) == 2


# --- metric consistency ------------------------------------------------

def test_a_full_length_perfect_hit_passes_on_its_own():
    ev = allele_evidence([hsp(1, 1000)], FULL, 0.8)
    assert ev["coverage"] == pytest.approx(1.0)
    assert ev["identity"] == pytest.approx(100.0)
    assert ev["split"] is False
    assert ev["n_copies"] == 1


def test_weaker_overlapping_copies_do_not_dilute_a_full_length_hit():
    """Regression, and the exact shape described in the review: five 80%
    duplicates must leave the 100% match's metrics untouched."""
    good = hsp(1, 1000, pident=100.0, bitscore=1800)
    weak = [hsp(1, 1000, sseqid=f"c{i}", pident=80.0, bitscore=700) for i in range(2, 7)]

    alone = allele_evidence([good], FULL, 0.8)
    with_duplicates = allele_evidence([good] + weak, FULL, 0.8)

    assert with_duplicates["identity"] == pytest.approx(alone["identity"])
    assert with_duplicates["coverage"] == pytest.approx(alone["coverage"])
    assert with_duplicates["score"] == pytest.approx(alone["score"])
    assert with_duplicates["identity"] == pytest.approx(100.0)
    assert with_duplicates["n_copies"] == 6, "the copies are still counted, just not merged"


def test_equivalent_duplicates_are_counted_once():
    a = hsp(1, 1000, sseqid="c1", bitscore=1800)
    b = hsp(1, 1000, sseqid="c2", bitscore=1800)
    ev = allele_evidence([a, b], FULL, 0.8)
    assert ev["score"] == pytest.approx(1800)
    assert ev["coverage"] == pytest.approx(1.0)
    assert ev["n_copies"] == 2


def test_overlapping_hsps_of_one_copy_do_not_double_count_coverage():
    ev = allele_evidence([hsp(1, 600, bitscore=1000), hsp(400, 1000, sstart=4000, bitscore=900)],
                         FULL, 0.8)
    assert ev["coverage"] == pytest.approx(0.6), "only the winning copy's span"
    assert ev["score"] == pytest.approx(1000)


def test_non_overlapping_fragments_are_chained_only_when_needed_and_flagged():
    """Split-contig support is retained, with its limitation made explicit."""
    ev = allele_evidence([hsp(1, 500, sseqid="c1", bitscore=900),
                          hsp(501, 1000, sseqid="c2", sstart=1, bitscore=900)],
                         FULL, 0.8)
    assert ev["coverage"] == pytest.approx(1.0)
    assert ev["score"] == pytest.approx(1800)
    assert ev["split"] is True, "a cross-contig chain is not proof of one intact gene"
    assert ev["contigs"] == ["c1", "c2"]


def test_a_sufficient_single_copy_is_not_chained_with_unrelated_fragments():
    """One copy already clears min_cov, so a stray fragment elsewhere neither
    adds to the score nor implies a longer gene."""
    ev = allele_evidence([hsp(1, 1000, sseqid="c1", bitscore=1800),
                          hsp(1, 200, sseqid="c9", bitscore=300)],
                         FULL, 0.8)
    assert ev["split"] is False
    assert ev["score"] == pytest.approx(1800)
    assert ev["contigs"] == ["c1"]


def test_identity_comes_from_the_same_hsps_as_the_score():
    """The invariant the old code broke: one evidence set, three metrics."""
    strong = hsp(1, 800, pident=99.0, bitscore=1400)
    weak = hsp(1, 800, sseqid="c2", pident=70.0, bitscore=400)
    ev = allele_evidence([strong, weak], FULL, 0.99)
    assert ev["score"] == pytest.approx(1400)
    assert ev["identity"] == pytest.approx(99.0)
    assert ev["aligned_len"] == 800


def test_subject_coordinates_and_orientation_are_retained():
    ev = allele_evidence([hsp(1, 1000, sseqid="c1", sstart=5000, strand="-")], FULL, 0.8)
    assert ev["strands"] == ["-"]
    assert ev["subject_spans"] == [("c1", 4001, 5000)]
