"""Shared HTTP client with retry/backoff for external services.

RepGenR2 talks to flaky, externally-operated services (the GTDB API and metadata
download, NCBI Entrez). Routing them through one process-wide ``requests.Session``
with a urllib3 retry policy gives every caller the same transient-failure
resilience -- retrying 429/5xx with exponential backoff and honoring
``Retry-After`` -- and the same clean error surface: a failed request raises
:class:`WorkdirError` naming the URL, never a bare ``requests`` traceback. The
streaming :func:`download` additionally verifies the byte count against
``Content-Length`` and writes through a ``.part`` file, so an interrupted
transfer never leaves a truncated file that later parses as if complete.
"""

from __future__ import annotations

import hashlib
import logging
import os
import re
from collections.abc import Iterator
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlsplit

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .errors import WorkdirError

_log = logging.getLogger(__name__)

# (connect, read) seconds. A server that does not accept a connection within
# the connect timeout is unreachable (blocked network, dead proxy); with the
# retries below that is reported in about two minutes instead of the twelve
# that a single 120 s timeout per attempt took. Reads keep the long timeout.
_CONNECT_TIMEOUT = 15
_READ_TIMEOUT = 120
_DEFAULT_TIMEOUT: tuple[float, float] = (_CONNECT_TIMEOUT, _READ_TIMEOUT)
Timeout = float | tuple[float, float]
_CHUNK = 1 << 20  # 1 MiB streaming chunks

_RETRY = Retry(
    total=5,
    backoff_factor=1.0,  # 0s, 1s, 2s, 4s, 8s between attempts
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=frozenset({"GET", "HEAD"}),
    respect_retry_after_header=True,
    raise_on_status=False,
)


class HTTPStatusError(WorkdirError):
    """The server answered with an error status (4xx/5xx after retries).

    ``status`` lets a caller tell "this resource does not exist" (404) from a
    network or server failure, e.g. to report an unknown taxon as user input.
    """

    def __init__(self, message: str, status: int):
        super().__init__(message)
        self.status = status


def _request_error(prefix: str, url: str, exc: requests.RequestException) -> WorkdirError:
    response = getattr(exc, "response", None)
    if isinstance(exc, requests.HTTPError) and response is not None:
        return HTTPStatusError(f"{prefix}: {url} ({exc})", response.status_code)
    return WorkdirError(f"{prefix}: {url} ({exc})")


@lru_cache(maxsize=1)
def session() -> requests.Session:
    """Process-wide session with the retry policy mounted on http(s)."""
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=_RETRY)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "repgenr"})
    return s


# The NCBI Datasets API host that the ``datasets`` CLI talks to. The CLI has
# no connect timeout of its own: on a blocked network each attempt waited about
# 8.5 minutes, three attempts in all, before the stage failed.
NCBI_DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/version"
_PROBE_TIMEOUT: tuple[float, float] = (_CONNECT_TIMEOUT, 30)


def _probe_session() -> requests.Session:
    """A session without the retry policy, for the single reachability request.

    ``trust_env`` is left on, so HTTPS_PROXY, HTTP_PROXY and NO_PROXY apply as
    they do for :func:`session` and for the tools that use the same network.
    """
    return requests.Session()


# Set to 1 to skip :func:`require_reachable`, e.g. where the probe's request
# is refused by a network policy that still lets the tool itself through.
SKIP_PROBE_ENV = "REPGENR_SKIP_NET_PROBE"


def require_reachable(url: str, *, what: str, timeout: Timeout = _PROBE_TIMEOUT) -> None:
    """Raise :class:`WorkdirError` (exit 3) when the host of ``url`` cannot be reached.

    One small GET request without retries, before handing the network to an
    external tool that retries slowly on its own. Only a failure to connect (a
    refused or reset connection, a failed name lookup, an unreachable proxy, a
    connect timeout) counts as unreachable. Any HTTP answer, an error status
    included, shows the host is reachable; so do a TLS error, which may only
    mean that requests does not trust a certificate the system store holds
    (the tool may use the system store), and a read timeout on a slow answer.
    ``REPGENR_SKIP_NET_PROBE=1`` skips the request.
    """
    if os.environ.get(SKIP_PROBE_ENV, "").strip().lower() in ("1", "true", "yes"):
        return
    host = urlsplit(url).hostname or url
    s = _probe_session()
    try:
        s.get(url, timeout=timeout, allow_redirects=False)
    except requests.exceptions.SSLError as exc:
        _log.debug("Reachability probe of %s: TLS error counted as reachable (%s)", host, exc)
    except requests.ConnectionError as exc:
        # ConnectTimeout and ProxyError are ConnectionErrors; SSLError is
        # handled above. ReadTimeout is not a ConnectionError.
        raise WorkdirError(
            f"Cannot reach {host}, which {what} needs ({exc}). Check the network "
            "connection and, behind a proxy, the HTTPS_PROXY and NO_PROXY settings "
            "(and REQUESTS_CA_BUNDLE for a proxy with its own certificate authority). "
            f"{SKIP_PROBE_ENV}=1 skips this check."
        ) from exc
    except requests.RequestException as exc:
        _log.debug("Reachability probe of %s: %s counted as reachable", host, exc)
    finally:
        s.close()


def _get(url: str, *, params: dict | None, timeout: Timeout) -> requests.Response:
    try:
        resp = session().get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        return resp
    except requests.RequestException as exc:
        raise _request_error("HTTP request failed", url, exc) from exc


def get_json(url: str, *, params: dict | None = None, timeout: Timeout = _DEFAULT_TIMEOUT) -> dict:
    """GET ``url`` and parse JSON; raises :class:`WorkdirError` on failure."""
    resp = _get(url, params=params, timeout=timeout)
    try:
        return resp.json()
    except ValueError as exc:
        raise WorkdirError(f"Expected JSON from {url} but could not parse it ({exc})") from exc


def get_text(url: str, *, params: dict | None = None, timeout: Timeout = _DEFAULT_TIMEOUT) -> str:
    """GET ``url`` and return the response body as text (status-checked)."""
    return _get(url, params=params, timeout=timeout).text


def get_bytes(
    url: str, *, params: dict | None = None, timeout: Timeout = _DEFAULT_TIMEOUT
) -> bytes:
    """GET ``url`` and return the response body as bytes (status-checked).

    The body is returned as sent: a file served as ``application/gzip``
    without a ``Content-Encoding`` header stays compressed.
    """
    return _get(url, params=params, timeout=timeout).content


def iter_lines(
    url: str, *, params: dict | None = None, timeout: Timeout = _DEFAULT_TIMEOUT
) -> Iterator[str]:
    """GET ``url`` and yield the body line by line as it arrives (status-checked).

    Lines end at newline only (a trailing carriage return is removed). The
    body is never held whole, so a large tabular answer costs the memory
    of what the caller keeps. A connection lost mid-body raises
    :class:`WorkdirError`, as a failed request does.
    """
    try:
        with session().get(url, params=params, timeout=timeout, stream=True) as resp:
            resp.raise_for_status()
            if resp.encoding is None:
                resp.encoding = "utf-8"
            # Split on newline only: str.splitlines would also split on
            # U+0085, U+2028, \v, \f and \x1c-\x1e, which may occur in a value.
            for line in resp.iter_lines(chunk_size=_CHUNK, decode_unicode=True, delimiter="\n"):
                yield line.rstrip("\r")
    except requests.RequestException as exc:
        raise _request_error("HTTP request failed", url, exc) from exc


def download(
    url: str,
    dest: Path,
    *,
    timeout: Timeout = _DEFAULT_TIMEOUT,
    logger: logging.Logger | None = None,
    md5: str | None = None,
) -> Path:
    """Stream ``url`` to ``dest``, verifying the size and writing atomically.

    Writes to ``<dest>.part`` and renames on success; a transfer that drops
    short of the server's ``Content-Length`` is deleted and raised as a
    :class:`WorkdirError` rather than left as a silently-truncated file.

    With ``md5`` (hex, any case) the digest is updated per chunk as the bytes
    are written and compared at the end, so the file is not read a second
    time. A mismatch deletes the ``.part`` file and raises the error that
    :func:`verify_md5` raises.
    """
    dest = Path(dest)
    tmp = dest.with_name(dest.name + ".part")
    written = 0
    expected = 0
    digest = hashlib.md5() if md5 else None
    try:
        with session().get(url, stream=True, timeout=timeout) as resp:
            resp.raise_for_status()
            expected = int(resp.headers.get("Content-Length", 0) or 0)
            with open(tmp, "wb") as fo:
                for chunk in resp.iter_content(chunk_size=_CHUNK):
                    fo.write(chunk)
                    written += len(chunk)
                    if digest is not None:
                        digest.update(chunk)
    except requests.RequestException as exc:
        tmp.unlink(missing_ok=True)
        raise _request_error("Download failed", url, exc) from exc

    if expected and written != expected:
        tmp.unlink(missing_ok=True)
        raise WorkdirError(f"Incomplete download: {url} got {written} of {expected} bytes.")
    if digest is not None and md5 is not None:
        actual = digest.hexdigest()
        if actual != md5.lower():
            tmp.unlink(missing_ok=True)
            raise WorkdirError(_md5_mismatch(dest.name, md5, actual))
    tmp.replace(dest)
    if logger is not None:
        logger.info("Downloaded %s (%d bytes)", dest.name, written)
    return dest


_MD5_RE = re.compile(r"^([0-9a-fA-F]{32})\s+\.?/?(.+)$")


def verify_md5_manifest(path: Path, manifest_url: str, *, logger: logging.Logger) -> bool:
    """Verify ``path`` against an md5sum-format manifest published beside it.

    Returns True when the checksum matches. An unavailable manifest or a file
    not listed in it is logged and skipped (False) -- upstream layouts vary --
    but an actual mismatch raises :class:`WorkdirError`: the download is
    corrupt and must not be parsed.
    """
    try:
        text = get_text(manifest_url)
    except WorkdirError as exc:
        logger.warning("Checksum manifest unavailable (%s); skipping verification", exc)
        return False
    expected = None
    for line in text.splitlines():
        m = _MD5_RE.match(line.strip())
        if m and Path(m.group(2)).name == path.name:
            expected = m.group(1).lower()
            break
    if expected is None:
        logger.warning("%s not listed in %s; skipping verification", path.name, manifest_url)
        return False
    verify_md5(path, expected)
    logger.info("Checksum verified for %s", path.name)
    return True


def verify_md5(path: Path, expected: str) -> None:
    """Raise :class:`WorkdirError` unless ``path`` hashes to ``expected`` (hex, any case)."""
    digest = hashlib.md5()
    with open(path, "rb") as fo:
        for chunk in iter(lambda: fo.read(_CHUNK), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected.lower():
        raise WorkdirError(_md5_mismatch(path.name, expected, actual))


def _md5_mismatch(name: str, expected: str, actual: str) -> str:
    return (
        f"Checksum mismatch for {name}: expected {expected.lower()}, got {actual}. "
        "The download is corrupt; delete it and re-run."
    )
