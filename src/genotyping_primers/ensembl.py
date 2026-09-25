"""Fetching the target stretch from a reference genome.

A minimal Ensembl REST client, built on ``urllib`` so the package needs nothing
beyond numpy. It covers the two things a deletion design needs: the sequence of
a region, and the name of the assembly that region came from.

Nothing here is required. Point ``--sequence`` at a local FASTA and the whole
module stays unused.
"""

from __future__ import annotations

import json
import logging
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

from .config import ENSEMBL_REST_URL

logger = logging.getLogger(__name__)

USER_AGENT = (
    "genotyping-primer-designer/0.1 "
    "(https://github.com/takeisue1270/Genotyping_Primer_Designer_for_Deletion)"
)

#: ``11:69000000-69001000``, ``chr11:69,000,000..69,001,000``, ``11:69000000+1000``.
_LOCUS_RE = re.compile(
    r"^(?:(?P<region>[\w.\-]+)\s*:\s*)?"
    r"(?P<start>[\d,_]+)"
    r"\s*(?P<op>-|\.\.|:|\+)\s*"
    r"(?P<end>[\d,_]+)$"
)


class EnsemblError(RuntimeError):
    """Raised when Ensembl refuses a request or returns something unusable."""


def parse_locus(spec: str) -> tuple[str, int, int]:
    """``"11:69,000,000-69,001,000"`` -> ``("11", 69000000, 69001000)``.

    The region name is optional, so ``"101-600"`` parses as ``("", 101, 600)``
    and addresses the supplied sequence directly. Both coordinates are 1-based
    inclusive, matching how a genome browser and a SnapGene map both report a
    selection. A ``+`` separator means "this many bases from the start", which
    is how deletions are often specified in a protocol.

    An end before the start is returned as given, not refused: on a plasmid it
    means the span crosses the origin. Only the caller knows whether the target
    is circular, so only the caller can reject it.
    """
    match = _LOCUS_RE.match(spec.strip())
    if not match:
        raise ValueError(
            f"cannot read {spec!r} as a region; expected forms are "
            "'101-600', '101..600', '101+500' or '11:69000000-69001000'"
        )
    region = (match.group("region") or "").strip()
    if region.lower().startswith("chr") and len(region) > 3:
        region = region[3:]
    start = int(match.group("start").replace(",", "").replace("_", ""))
    end = int(match.group("end").replace(",", "").replace("_", ""))
    if match.group("op") == "+":
        end = start + end - 1
    if start < 1 or end < 1:
        raise ValueError(f"{spec!r}: coordinates are 1-based, so both must be >= 1")
    return region, start, end


@dataclass
class EnsemblClient:
    """A throttled JSON/plain-text client for the Ensembl REST API."""

    server: str = ENSEMBL_REST_URL
    max_requests_per_second: float = 12.0
    timeout: float = 60.0
    max_retries: int = 4

    def __post_init__(self) -> None:
        self._min_interval = 1.0 / self.max_requests_per_second
        self._last_request = 0.0
        self._cache: dict[str, str] = {}

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request = time.monotonic()

    def _request(self, endpoint: str, content_type: str, max_retries: int | None = None) -> str:
        url = f"{self.server}{endpoint}"
        request = urllib.request.Request(
            url,
            headers={
                "Content-Type": content_type,
                "Accept": content_type,
                "User-Agent": USER_AGENT,
            },
        )

        attempts = self.max_retries if max_retries is None else max_retries
        last_error: Exception | None = None
        for attempt in range(attempts):
            self._throttle()
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read().decode("utf-8")
            except urllib.error.HTTPError as exc:
                last_error = exc
                if exc.code == 429:
                    wait = float(exc.headers.get("Retry-After", 1.0))
                    logger.warning("Ensembl rate limit hit; waiting %.1f s", wait)
                    time.sleep(wait)
                    continue
                if exc.code in (500, 502, 503, 504):
                    wait = 2.0**attempt
                    logger.warning("Ensembl returned %s; retrying in %.0f s", exc.code, wait)
                    time.sleep(wait)
                    continue
                body = exc.read().decode("utf-8", errors="replace")[:400]
                raise EnsemblError(f"{exc.code} {exc.reason} for {url}\n{body}") from exc
            except urllib.error.URLError as exc:
                last_error = exc
                wait = 2.0**attempt
                logger.warning(
                    "Network error talking to Ensembl (%s); retrying in %.0f s", exc.reason, wait
                )
                time.sleep(wait)

        raise EnsemblError(f"giving up on {url} after {attempts} attempts: {last_error}")

    def assembly_name(self, species: str) -> str:
        """Default assembly name for a species, e.g. ``GRCm39``; ``""`` if unavailable.

        Best-effort: the ``/info/assembly`` endpoints are the flakiest part of
        the public service, and they only supply a label for the report.
        """
        try:
            info = json.loads(self._request(f"/info/assembly/{species}", "application/json", 2))
        except (EnsemblError, json.JSONDecodeError) as exc:
            logger.info("assembly name unavailable for %s: %s", species, exc)
            return ""
        if not isinstance(info, dict):
            return ""
        return str(info.get("default_coord_system_version") or info.get("assembly_name") or "")

    def sequence_region(self, species: str, region: str, start: int, end: int) -> str:
        """Fetch ``region:start..end`` (1-based inclusive) from the + strand."""
        if start < 1:
            raise ValueError("start must be >= 1")
        if end < start:
            raise ValueError("end must be >= start")
        key = f"{species}/{region}:{start}-{end}"
        if key not in self._cache:
            locator = f"{region}:{start}..{end}:1"
            text = self._request(f"/sequence/region/{species}/{locator}", "text/plain")
            self._cache[key] = "".join(text.split()).upper()
        return self._cache[key]


__all__ = ["EnsemblClient", "EnsemblError", "parse_locus"]
