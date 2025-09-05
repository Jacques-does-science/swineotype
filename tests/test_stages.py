from unittest.mock import patch
from pathlib import Path
from swineotype.stages import stage1_score

from swineotype.stages import stage2_resolver_call, interpret_resolver

@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
@patch("swineotype.stages.parse_whitelist_headers")
def test_stage1_score(mock_parse_whitelist_headers, mock_make_db_if_needed, mock_run_blast, mock_ensure_tool):
    mock_parse_whitelist_headers.return_value = (
        {"q1": "1", "q2": "14", "q3": "2"},
        {"q1": "wzx", "q2": "wzy", "q3": "wzx"},
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
@patch("swineotype.stages.run")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_1_vs_14(mock_make_db, mock_run_blast, mock_run, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    qseqid = "cps14K|pair=1_vs_14|pos=481|baseA=T|A=14|B=1"
    blast_out = f"{qseqid}\ts1\t100\t1000\t1000\t0\t2000\t1\t1000\t1\t1000"
    mock_run_blast.return_value = blast_out
    mock_run.return_value = ">s1:1-1000\nT"

    config = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}
    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), config, "1_vs_14")

    assert result["ref_id"] == qseqid
    assert result["base"] == "T"

    final_sero = interpret_resolver(result, config)
    assert final_sero == "1"

@patch("swineotype.stages.ensure_tool")
@patch("swineotype.stages.run")
@patch("swineotype.stages.run_blast")
@patch("swineotype.stages.make_db_if_needed")
def test_stage2_resolver_call_2_vs_1_2(mock_make_db, mock_run_blast, mock_run, mock_ensure_tool):
    mock_make_db.return_value = "db_prefix"
    qseqid = "cps2K|pair=2_vs_1_2|pos=481|baseA=G|A=1/2|B=2"
    blast_out = f"{qseqid}\ts1\t100\t1000\t1000\t0\t2000\t1\t1000\t1\t1000"
    mock_run_blast.return_value = blast_out
    mock_run.return_value = ">s1:1-1000\nG"

    config = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}
    result = stage2_resolver_call("assembly.fa", "resolver.fa", 4, Path("run_dir"), config, "2_vs_1_2")

    assert result["ref_id"] == qseqid
    assert result["base"] == "G"

    final_sero = interpret_resolver(result, config)
    assert final_sero == "2"
