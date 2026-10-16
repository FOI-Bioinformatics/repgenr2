"""Tests for the in-stage parallel_map helper."""

from __future__ import annotations

import threading
import time

import pytest

from repgenr.core.executors import parallel_map


def test_order_preserved() -> None:
    assert parallel_map(lambda x: x * 2, [1, 2, 3, 4], workers=3) == [2, 4, 6, 8]


def test_sequential_when_single_worker() -> None:
    assert parallel_map(str, [1, 2, 3], workers=1) == ["1", "2", "3"]


def test_empty() -> None:
    assert parallel_map(lambda x: x, [], workers=4) == []


def test_exception_propagates() -> None:
    def boom(x):
        if x == 2:
            raise ValueError("boom")
        return x

    with pytest.raises(ValueError, match="boom"):
        parallel_map(boom, [1, 2, 3], workers=2)


def test_actually_concurrent() -> None:
    # 4 tasks that each sleep; with >=4 workers they overlap, so max concurrency > 1.
    active = 0
    peak = 0
    lock = threading.Lock()

    def task(_):
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        time.sleep(0.05)
        with lock:
            active -= 1
        return 1

    parallel_map(task, list(range(4)), workers=4)
    assert peak >= 2


def test_parallel_each_delivers_on_the_calling_thread() -> None:
    from repgenr.core.executors import parallel_each

    caller = threading.get_ident()
    seen: list[tuple[int, int, int]] = []
    parallel_each(
        lambda x: x * 2, [1, 2, 3, 4], 3, lambda x, r: seen.append((x, r, threading.get_ident()))
    )
    assert sorted((x, r) for x, r, _ in seen) == [(1, 2), (2, 4), (3, 6), (4, 8)]
    assert {t for _, _, t in seen} == {caller}


def test_parallel_each_sequential_and_empty() -> None:
    from repgenr.core.executors import parallel_each

    seen: list[int] = []
    parallel_each(str, [1, 2], 1, lambda x, _r: seen.append(x))
    parallel_each(str, [], 4, lambda x, _r: seen.append(x))
    assert seen == [1, 2]


def test_parallel_each_failure_keeps_running_results_and_skips_queued() -> None:
    """The item running beside a failure finishes and is delivered; queued
    items never start; the error is raised afterwards."""
    from repgenr.core.executors import parallel_each

    started: list[int] = []
    delivered: list[int] = []
    gate = threading.Event()

    def work(x: int) -> int:
        started.append(x)
        if x == 1:
            gate.wait(5)  # still running when item 2 fails
            return x
        if x == 2:
            time.sleep(0.05)
            gate.set()
            raise ValueError("boom")
        return x

    with pytest.raises(ValueError, match="boom"):
        parallel_each(work, [1, 2, 3, 4, 5], 2, lambda x, _r: delivered.append(x))
    assert delivered == [1]
    assert sorted(started) == [1, 2]
