"""Packaged, source-linked evidence catalogue for mathematical topics.

The manifest records implementation scope; it is not a certificate that a
complete mathematical field is implemented or formally proved.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from importlib.resources import files
import json
from typing import Any


@lru_cache(maxsize=1)
def _cached_evidence() -> dict[str, Any]:
    """Validate the packaged JSON once; never expose this cached object."""
    content = files("calc").joinpath("topic_manifest.json").read_text(encoding="utf-8")
    parsed = json.loads(content)
    if not isinstance(parsed, dict) or parsed.get("schema_version") != 1:
        raise ValueError("unsupported topic evidence manifest schema")
    topics = parsed.get("topics")
    if not isinstance(topics, list) or not all(isinstance(item, dict) for item in topics):
        raise ValueError("invalid topic evidence records")
    return parsed


def load_evidence() -> dict[str, Any]:
    """Return a fresh, mutable manifest snapshot isolated from other callers."""
    return deepcopy(_cached_evidence())


def get_evidence(slug: str) -> dict[str, Any]:
    """Return a deep copy of one record; raise KeyError for unknown slugs."""
    for record in _cached_evidence()["topics"]:
        if record.get("slug") == slug:
            return deepcopy(record)
    raise KeyError(slug)

