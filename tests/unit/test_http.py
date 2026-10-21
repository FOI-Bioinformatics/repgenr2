"""Tests for the shared HTTP client (no network).

The session is replaced with a fake so the JSON/text/download helpers and their
error surfaces are exercised deterministically: status errors and unparseable
bodies become WorkdirError, and a short download is rejected with no leftover
file.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
import requests

from repgenr.core import http
from repgenr.core.errors import WorkdirError

_LOG = logging.getLogger("test")
# Taken at import, before the autouse fixture in tests/conftest.py stubs it.
_REAL_PROBE_SESSION = http._probe_session


class _FakeResp:
    def __init__(self, *, json_data=None, text="", content=b"", headers=None, status_exc=None):
        self._json = json_data
        self.text = text
        self._content = content
        self.headers = headers or {}
        self._status_exc = status_exc

    @property
    def content(self):
        return self._content

    def raise_for_status(self):
        if self._status_exc is not None:
            raise self._status_exc

    def json(self):
        if isinstance(self._json, ValueError):
            raise self._json
        return self._json

    def iter_content(self, chunk_size=1):
        for i in range(0, len(self._content), chunk_size):
            yield self._content[i : i + chunk_size]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _FakeSession:
    def __init__(self, resp):
        self._resp = resp

    def get(self, url, **kw):
        return self._resp


def _patch(monkeypatch, resp):
    monkeypatch.setattr(http, "session", lambda: _FakeSession(resp))


def test_get_json_ok(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(json_data={"ok": 1}))
    assert http.get_json("https://x/y") == {"ok": 1}


def test_get_json_bad_body(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(json_data=ValueError("nope")))
    with pytest.raises(WorkdirError, match="could not parse"):
        http.get_json("https://x/y")


def test_get_text_ok(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(text="hello"))
    assert http.get_text("https://x/y") == "hello"


def test_get_bytes_returns_the_body_as_sent(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(content=b"\x1f\x8bgz"))
    assert http.get_bytes("https://x/y") == b"\x1f\x8bgz"


def test_status_error_becomes_workdir_error(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.HTTPError("500")))
    with pytest.raises(WorkdirError, match="HTTP request failed"):
        http.get_text("https://x/y")


def test_download_ok(monkeypatch, tmp_path: Path) -> None:
    body = b"A" * 2048
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    dest = tmp_path / "f.gz"
    http.download("https://x/f.gz", dest)
    assert dest.read_bytes() == body
    assert not (tmp_path / "f.gz.part").exists()  # temp renamed away


def test_download_truncated_is_rejected(monkeypatch, tmp_path: Path) -> None:
    # server promises 4096 bytes but the body is short -> reject, leave nothing
    _patch(monkeypatch, _FakeResp(content=b"A" * 10, headers={"Content-Length": "4096"}))
    dest = tmp_path / "f.gz"
    with pytest.raises(WorkdirError, match="Incomplete download"):
        http.download("https://x/f.gz", dest)
    assert not dest.exists()
    assert not (tmp_path / "f.gz.part").exists()


def test_download_request_error(monkeypatch, tmp_path: Path) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.ConnectionError("reset")))
    with pytest.raises(WorkdirError, match="Download failed"):
        http.download("https://x/f.gz", tmp_path / "f.gz")
    assert not (tmp_path / "f.gz.part").exists()


def _count_reads(monkeypatch) -> list[str]:
    """Record every opening of a file for reading (the md5 must not need one)."""
    import builtins

    reads: list[str] = []
    real_open = builtins.open

    def counting_open(file, mode="r", *args, **kwargs):  # noqa: ANN001, ANN202
        if "r" in mode:
            reads.append(str(file))
        return real_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", counting_open)
    return reads


def test_download_with_md5_hashes_while_writing(monkeypatch, tmp_path: Path) -> None:
    import hashlib

    body = bytes(range(256)) * 9000  # more than one 1 MiB chunk
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    reads = _count_reads(monkeypatch)
    verified: list[Path] = []
    monkeypatch.setattr(http, "verify_md5", lambda path, md5: verified.append(path))
    dest = tmp_path / "f.gz"
    digest = hashlib.md5(body).hexdigest().upper()  # any case
    assert http.download("https://x/f.gz", dest, md5=digest) == dest
    # The file was opened for writing only: never read back, never verified again.
    assert reads == []
    assert verified == []
    assert dest.read_bytes() == body
    assert not (tmp_path / "f.gz.part").exists()


def test_download_with_a_wrong_md5_leaves_nothing(monkeypatch, tmp_path: Path) -> None:
    body = b"A" * 2048
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    dest = tmp_path / "f.gz"
    with pytest.raises(WorkdirError) as info:
        http.download("https://x/f.gz", dest, md5="0" * 32)
    # The same text verify_md5 gives.
    assert str(info.value).startswith(f"Checksum mismatch for f.gz: expected {'0' * 32}, got ")
    assert "The download is corrupt; delete it and re-run." in str(info.value)
    assert not dest.exists()
    assert not (tmp_path / "f.gz.part").exists()


def test_download_without_md5_does_not_hash(monkeypatch, tmp_path: Path) -> None:
    import hashlib

    body = b"A" * 2048
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    real_md5 = hashlib.md5
    made: list[int] = []

    def md5(*args, **kwargs):  # noqa: ANN202
        made.append(1)
        return real_md5(*args, **kwargs)

    monkeypatch.setattr(http.hashlib, "md5", md5)
    dest = tmp_path / "f.gz"
    http.download("https://x/f.gz", dest)
    assert dest.read_bytes() == body and made == []


def test_verify_md5_and_download_give_the_same_mismatch_text(monkeypatch, tmp_path: Path) -> None:
    body = b"payload"
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    with pytest.raises(WorkdirError) as during:
        http.download("https://x/f.gz", tmp_path / "f.gz", md5="ab" * 16)
    local = tmp_path / "f.gz"
    local.write_bytes(body)
    with pytest.raises(WorkdirError) as after:
        http.verify_md5(local, "ab" * 16)
    assert str(during.value) == str(after.value)


# --- MD5 manifest verification ------------------------------------------------


def test_verify_md5_manifest_ok(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    import hashlib

    digest = hashlib.md5(b"payload").hexdigest()
    manifest = f"{digest}  ./data.tsv.gz\nffff  ./other.txt\n"
    monkeypatch.setattr(http, "get_text", lambda url, **k: manifest)
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is True


def test_verify_md5_manifest_mismatch_raises(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    manifest = "0" * 32 + "  ./data.tsv.gz\n"
    monkeypatch.setattr(http, "get_text", lambda url, **k: manifest)
    with pytest.raises(WorkdirError, match="hecksum"):
        http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG)


def test_verify_md5_manifest_unavailable_skips(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")

    def boom(url, **k):
        raise WorkdirError("404")

    monkeypatch.setattr(http, "get_text", boom)
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is False


def test_verify_md5_manifest_unlisted_file_skips(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    monkeypatch.setattr(http, "get_text", lambda url, **k: "aa" * 16 + "  ./other.txt\n")
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is False


def test_verify_md5_accepts_the_right_digest_and_rejects_a_wrong_one(tmp_path) -> None:
    import hashlib

    from repgenr.core.errors import WorkdirError
    from repgenr.core.http import verify_md5

    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(b"@r1\nACGT\n+\nIIII\n")
    good = hashlib.md5(path.read_bytes()).hexdigest()
    verify_md5(path, good.upper())  # case-insensitive, no return value needed
    with pytest.raises(WorkdirError, match="Checksum mismatch"):
        verify_md5(path, "0" * 32)


def _http_error(status: int) -> requests.HTTPError:
    response = requests.Response()
    response.status_code = status
    return requests.HTTPError(f"{status} error", response=response)


@pytest.mark.parametrize("status", [404, 503])
def test_status_error_carries_the_status(monkeypatch, tmp_path: Path, status: int) -> None:
    """A caller can tell an unknown resource (404) from a failing server."""
    _patch(monkeypatch, _FakeResp(status_exc=_http_error(status)))
    with pytest.raises(http.HTTPStatusError) as info:
        http.get_json("https://x/y")
    assert info.value.status == status
    with pytest.raises(http.HTTPStatusError) as info:
        http.download("https://x/y", tmp_path / "f")
    assert info.value.status == status


def test_connection_error_is_not_a_status_error(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.ConnectionError("refused")))
    with pytest.raises(WorkdirError) as info:
        http.get_json("https://x/y")
    assert not isinstance(info.value, http.HTTPStatusError)


def test_requests_use_a_short_connect_timeout(monkeypatch, tmp_path: Path) -> None:
    """A blocked network is reported in minutes: through an unreachable proxy
    one GTDB API request waited 481 s (six 120 s connect attempts)."""
    seen = []

    class _Recording(_FakeSession):
        def get(self, url, **kw):
            seen.append(kw["timeout"])
            return super().get(url, **kw)

    body = b"x"
    resp = _FakeResp(json_data={}, content=body, headers={"Content-Length": "1"})
    monkeypatch.setattr(http, "session", lambda: _Recording(resp))
    http.get_json("https://x/y")
    http.download("https://x/y", tmp_path / "f")
    assert seen and all(isinstance(t, tuple) and t[0] <= 30 and t[1] >= 120 for t in seen)


class _ProbeSession:
    def __init__(self, exc: Exception | None = None):
        self.exc = exc
        self.calls: list[tuple[str, dict]] = []

    def get(self, url, **kw):
        self.calls.append((url, kw))
        if self.exc is not None:
            raise self.exc
        return None

    def close(self) -> None:
        pass


def _status_response(status: int) -> requests.Response:
    resp = requests.Response()
    resp.status_code = status
    resp.url = http.NCBI_DATASETS_URL
    return resp


@pytest.mark.parametrize("status", [404, 503])
def test_an_error_status_counts_as_reachable(monkeypatch, status: int) -> None:
    """A server that answers, even with 4xx/5xx, is reachable; datasets handles the rest."""

    class _Answers(_ProbeSession):
        def get(self, url, **kw):
            super().get(url, **kw)
            return _status_response(status)

    probe = _Answers()
    monkeypatch.setattr(http, "_probe_session", lambda: probe)
    http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    assert len(probe.calls) == 1


@pytest.mark.parametrize(
    "exc",
    [
        # requests trusts its own CA bundle, datasets (Go) the system store: a
        # proxy CA in the system store alone fails here but not in datasets.
        requests.exceptions.SSLError("certificate verify failed"),
        requests.ReadTimeout("read timed out"),
    ],
)
def test_tls_errors_and_read_timeouts_count_as_reachable(monkeypatch, exc) -> None:
    monkeypatch.setattr(http, "_probe_session", lambda: _ProbeSession(exc))
    http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")


def test_an_unreachable_proxy_is_unreachable(monkeypatch) -> None:
    exc = requests.exceptions.ProxyError("Unable to connect to proxy")
    monkeypatch.setattr(http, "_probe_session", lambda: _ProbeSession(exc))
    with pytest.raises(WorkdirError, match="REQUESTS_CA_BUNDLE") as info:
        http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    assert http.SKIP_PROBE_ENV in str(info.value)


def test_the_probe_can_be_switched_off(monkeypatch) -> None:
    probe = _ProbeSession(requests.ConnectionError("refused"))
    monkeypatch.setattr(http, "_probe_session", lambda: probe)
    monkeypatch.setenv(http.SKIP_PROBE_ENV, "1")
    http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    assert probe.calls == []


def test_require_reachable_passes_when_the_host_answers(monkeypatch) -> None:
    probe = _ProbeSession()
    monkeypatch.setattr(http, "_probe_session", lambda: probe)
    http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    [(url, kw)] = probe.calls
    assert url == http.NCBI_DATASETS_URL
    # One attempt with the shared 15 s connect timeout.
    assert kw["timeout"][0] == 15


@pytest.mark.parametrize(
    "exc",
    [requests.ConnectTimeout("timed out"), requests.ConnectionError("refused")],
)
def test_require_reachable_names_the_host_when_it_does_not_answer(monkeypatch, exc) -> None:
    """datasets on a blocked network made three attempts of about 8.5 minutes each."""
    monkeypatch.setattr(http, "_probe_session", lambda: _ProbeSession(exc))
    with pytest.raises(WorkdirError, match="api.ncbi.nlm.nih.gov") as info:
        http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    assert info.value.exit_code == 3
    assert "HTTPS_PROXY" in str(info.value)


def test_probe_session_reads_the_proxy_environment(monkeypatch) -> None:
    """The probe goes through the same proxies as datasets and core.http."""
    monkeypatch.setenv("HTTPS_PROXY", "http://proxy.example:3128")
    monkeypatch.delenv("NO_PROXY", raising=False)
    monkeypatch.delenv("no_proxy", raising=False)
    s = _REAL_PROBE_SESSION()
    try:
        settings = s.merge_environment_settings(http.NCBI_DATASETS_URL, {}, None, None, None)
    finally:
        s.close()
    assert settings["proxies"]["https"] == "http://proxy.example:3128"


def _streamed(body: bytes, *, status: int = 200, chunks: int = 7) -> requests.Response:
    """A real requests.Response whose body arrives in small chunks."""
    import io

    resp = requests.Response()
    resp.status_code = status
    resp.url = "https://x/y"
    resp.raw = io.BufferedReader(io.BytesIO(body), buffer_size=chunks)  # type: ignore[assignment]
    return resp


def test_iter_lines_yields_the_body_line_by_line(monkeypatch) -> None:
    seen: dict = {}

    class _Session:
        def get(self, url, **kw):
            seen.update(kw)
            return _streamed("a\tb\nÅ1\t2\n3\t\n".encode())

    monkeypatch.setattr(http, "session", lambda: _Session())
    # A trailing newline ends in an empty line, which the TSV reader skips.
    assert list(http.iter_lines("https://x/y", params={"q": 1})) == ["a\tb", "Å1\t2", "3\t", ""]
    assert seen["stream"] is True and seen["params"] == {"q": 1}


def test_iter_lines_status_and_connection_errors_are_workdir_errors(monkeypatch) -> None:
    _patch(monkeypatch, _streamed(b"", status=400))
    with pytest.raises(http.HTTPStatusError, match="HTTP request failed"):
        list(http.iter_lines("https://x/y"))

    class _Cut:
        status_code = 200
        encoding = "utf-8"

        def raise_for_status(self):
            pass

        def iter_lines(self, **kw):
            yield "a\tb"
            raise requests.exceptions.ChunkedEncodingError("connection broken")

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    _patch(monkeypatch, _Cut())
    with pytest.raises(WorkdirError, match="connection broken"):
        list(http.iter_lines("https://x/y"))


def test_iter_lines_splits_on_newline_only(monkeypatch) -> None:
    """Characters str.splitlines treats as line ends stay inside a value."""
    name = "Bacillus\x85a b\vc\fd\x1ce\x1df\x1eg"
    body = f"run_accession\tscientific_name\r\nR1\t{name}\r\nR2\tplain\n".encode()
    _patch(monkeypatch, _streamed(body))
    assert list(http.iter_lines("https://x/y")) == [
        "run_accession\tscientific_name",
        f"R1\t{name}",
        "R2\tplain",
        "",
    ]
