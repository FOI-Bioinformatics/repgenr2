"""FastTree builder: the thread budget reaches FastTreeMP through OMP_NUM_THREADS."""

from __future__ import annotations

import logging
from pathlib import Path

import repgenr.treebuilders.fasttree as mod
from repgenr.treebuilders.base import TreeParams


def test_threads_reach_fasttree_as_omp_env(monkeypatch, tmp_path: Path) -> None:
    seen: dict = {}

    def run_tool(caps, command, *, logger, tool_env=None, stdout_path=None, **kwargs):
        seen["env"] = dict(tool_env or {})
        Path(stdout_path).write_text("(a,b);", encoding="utf-8")
        return 0

    monkeypatch.setattr(mod, "run_tool", run_tool)
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n>b\nAG\n", encoding="utf-8")
    mod.FasttreeBuilder().build(msa, tmp_path / "out", TreeParams(threads=5), logging.getLogger())
    assert seen["env"].get("OMP_NUM_THREADS") == "5"
