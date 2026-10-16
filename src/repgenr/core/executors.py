"""Lightweight in-stage parallelism.

The parallel work in RepGenR2 stages is subprocess-bound (progressiveMauve,
dRep, skani, ...), so a thread pool is the right tool: the GIL is released while
each external tool runs, and threads avoid the pickling constraints that
``multiprocessing`` imposes on closures/adapters. Cross-stage and HPC
parallelism is handled separately by the Nextflow layer.
"""

from __future__ import annotations

import logging
import threading
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed


def parallel_map[T, R](
    fn: Callable[[T], R],
    items: Iterable[T],
    workers: int,
    *,
    logger: logging.Logger | None = None,
) -> list[R]:
    """Apply ``fn`` to each item, up to ``workers`` at a time, preserving order.

    Runs sequentially when ``workers <= 1`` or there is at most one item. The
    first worker exception (or a termination signal) cancels the queued items
    and propagates once the running ones settle.
    """
    items_list: Sequence[T] = list(items)
    if workers <= 1 or len(items_list) <= 1:
        return [fn(item) for item in items_list]

    results: list[R | None] = [None] * len(items_list)
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        future_to_index = {pool.submit(fn, item): i for i, item in enumerate(items_list)}
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
    except BaseException:
        # A failure or a termination signal: queued items must not start;
        # the running ones settle (their tools are stopped by the handler).
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    return results  # type: ignore[return-value]


def parallel_each[T, R](
    fn: Callable[[T], R],
    items: Iterable[T],
    workers: int,
    on_done: Callable[[T, R], None],
) -> None:
    """Apply ``fn`` to each item, up to ``workers`` at a time, and hand each
    result to ``on_done`` on the calling thread as it completes.

    For work whose results must be recorded as they arrive (an index written
    after every genome), so an interrupted run keeps what finished. On a
    worker failure the queued items are cancelled, the running ones settle and
    their results still reach ``on_done``, and the first error is raised. A
    termination signal cancels the queued items and propagates once the
    running ones settle. Runs sequentially when ``workers <= 1``.
    """
    items_list: Sequence[T] = list(items)
    if workers <= 1 or len(items_list) <= 1:
        for item in items_list:
            on_done(item, fn(item))
        return

    # Set by the first failing worker, so a queued item that a freed worker
    # picks up before this thread sees the failure does not start.
    stop = threading.Event()

    def guarded(item: T) -> R:
        if stop.is_set():
            raise _Skipped
        try:
            return fn(item)
        except BaseException:
            stop.set()
            raise

    pool = ThreadPoolExecutor(max_workers=workers)
    error: BaseException | None = None
    try:
        future_to_item = {pool.submit(guarded, item): item for item in items_list}
        for future in as_completed(future_to_item):
            if future.cancelled():
                continue
            exc = future.exception()
            if isinstance(exc, _Skipped):
                continue
            if exc is not None:
                if error is None:
                    error = exc
                    for pending in future_to_item:
                        pending.cancel()
                continue
            on_done(future_to_item[future], future.result())
    except BaseException:
        stop.set()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    if error is not None:
        raise error


class _Skipped(Exception):
    """A queued item not started because another item failed."""
