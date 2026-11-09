"""Support values of tree nodes."""

from __future__ import annotations


def parse_support(label: str | None) -> float | None:
    """The support value of a node label, None when it holds none.

    IQ-TREE run with both ``-alrt`` and ``-B`` writes ``SH-aLRT/UFBoot``
    (for example ``87.3/98``); the last field, the ultrafast bootstrap, is
    read, so the value is the one a ``-B``-only run writes. A label with a
    single field is read as it is.
    """
    if label is None:
        return None
    try:
        return float(label.rsplit("/", 1)[-1])
    except ValueError:
        return None
