"""Packaged, source-linked evidence catalogue for mathematical topics.

The manifest records implementation scope; it is not a certificate that a
complete mathematical field is implemented or formally proved.
"""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files
import json
from typing import Any


@lru_cache(maxsize=1)
def load_evidence() -> dict[str, Any]:
    """Load the single packaged JSON manifest, checking its top-level shape."""
    content = files("calc").joinpath("topic_manifest.json").read_text(encoding="utf-8")
    parsed = json.loads(content)
    if not isinstance(parsed, dict) or parsed.get("schema_version") != 1:
        raise ValueError("unsupported topic evidence manifest schema")
    topics = parsed.get("topics")
    if not isinstance(topics, list) or not all(isinstance(item, dict) for item in topics):
        raise ValueError("invalid topic evidence records")
    return parsed


def get_evidence(slug: str) -> dict[str, Any]:
    """Return one evidence record; raise KeyError for unknown slugs."""
    for record in load_evidence()["topics"]:
        if record.get("slug") == slug:
            return dict(record)
    raise KeyError(slug)

