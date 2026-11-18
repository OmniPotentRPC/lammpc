#!/usr/bin/env python3
"""Check local and external links in the maintained documentation."""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*"))]
MARKDOWN_LINK = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
BARE_URL = re.compile(r"https?://[^\s>)\]]+")
DOI_URL = re.compile(r"https://doi\.org/(10\.[0-9]{4,9}/[^\s>)\]]+)")


def request(url: str) -> bytes:
    headers = {"User-Agent": "lammpc-doc-link-check/1"}
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=30
            ) as response:
                if response.status >= 400:
                    raise RuntimeError(f"HTTP {response.status}")
                return response.read()
        except (OSError, RuntimeError, urllib.error.URLError) as error:
            last_error = error
            if attempt != 2:
                time.sleep(2**attempt)
    raise RuntimeError(str(last_error))


def check_local(source: Path, target: str) -> None:
    path_text = urllib.parse.unquote(target.split("#", 1)[0])
    if not path_text:
        return
    candidate = (source.parent / path_text).resolve()
    if not candidate.exists():
        raise RuntimeError(f"{source.relative_to(ROOT)}: missing {target}")


def main() -> int:
    external: set[str] = set()
    dois: set[str] = set()
    errors: list[str] = []

    for source in SOURCES:
        if not source.is_file():
            continue
        text = source.read_text(encoding="utf-8")
        targets = set(MARKDOWN_LINK.findall(text))
        targets.update(BARE_URL.findall(text))
        for target in targets:
            if target.startswith(("http://", "https://")):
                if not target.startswith("https://doi.org/"):
                    external.add(target)
            elif not target.startswith(("mailto:", "#")):
                try:
                    check_local(source, target)
                except RuntimeError as error:
                    errors.append(str(error))
        dois.update(DOI_URL.findall(text))

    for url in sorted(external):
        try:
            request(url)
            print(f"ok {url}")
        except RuntimeError as error:
            errors.append(f"{url}: {error}")

    for doi in sorted(dois):
        url = f"https://doi.org/api/handles/{urllib.parse.quote(doi, safe='/')}"
        try:
            payload = json.loads(request(url))
            if payload.get("responseCode") != 1:
                raise RuntimeError(f"DOI response code {payload.get('responseCode')}")
            print(f"ok doi:{doi}")
        except (RuntimeError, ValueError) as error:
            errors.append(f"doi:{doi}: {error}")

    if errors:
        for error in errors:
            print(f"error {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
