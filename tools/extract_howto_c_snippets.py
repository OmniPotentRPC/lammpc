#!/usr/bin/env python3
"""Extract C code fences from docs/how-to.md into a build-time smoke test.

Only fenced ```c blocks in the input Markdown file are extracted. Each snippet
is wrapped in a function with dummy 1-byte buffers and a NULL session so
renamed symbols, changed parameter lists, or stale struct fields fail the
build. Other doc pages (such as docs/hard-boundary.org) and live LAMMPS
session execution are outside this check.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

SESSION_DECL_RE = re.compile(r"\bLammpcSession\s*\*\s*session\b")
OUTPUT_DECL_RE = re.compile(r"\bvoid\s*\*\s*output\b")


def extract_c_blocks(markdown: str) -> list[tuple[int, str]]:
    blocks: list[tuple[int, str]] = []
    in_c_fence = False
    start_line = 0
    lines: list[str] = []

    for lineno, line in enumerate(markdown.splitlines(), start=1):
        stripped = line.strip()
        if not in_c_fence:
            if stripped == "```c":
                in_c_fence = True
                start_line = lineno + 1
                lines = []
        else:
            if stripped == "```":
                blocks.append((start_line, "\n".join(lines)))
                in_c_fence = False
            else:
                lines.append(line)

    if in_c_fence:
        raise RuntimeError(f"unclosed ```c fence starting at line {start_line - 1}")
    if not blocks:
        raise RuntimeError("no ```c fences found")
    return blocks


def render_c(input_path: Path, blocks: list[tuple[int, str]]) -> str:
    source_label = input_path.as_posix().replace("\\", "\\\\").replace('"', '\\"')
    parts = [
        "/* Generated from docs/how-to.md at build time; do not edit. */",
        '#include "lammpc.h"',
        "",
        "#include <stdint.h>",
        "#include <stdio.h>",
        "#include <stdlib.h>",
        "",
    ]

    for index, (start_line, body) in enumerate(blocks):
        declares_session = bool(SESSION_DECL_RE.search(body))
        declares_output = bool(OUTPUT_DECL_RE.search(body))
        parts.extend(
            [
                f"static int snippet_{index}(void) {{",
                "  const uint8_t config_bytes[1] = {0};",
                "  const size_t config_size = sizeof(config_bytes);",
                "  const uint8_t input[1] = {0};",
                "  const size_t input_size = sizeof(input);",
                "  const double point[3] = {0.0, 0.0, 0.0};",
                "  const double normal[3] = {0.0, 0.0, 1.0};",
                "  const size_t coordinate_count = 3;",
                "  (void)config_bytes;",
                "  (void)config_size;",
                "  (void)input;",
                "  (void)input_size;",
                "  (void)point;",
                "  (void)normal;",
                "  (void)coordinate_count;",
            ]
        )
        if not declares_session:
            parts.extend(["  LammpcSession *session = NULL;", "  (void)session;"])
        parts.append(f'#line {start_line} "{source_label}"')
        for line in body.splitlines():
            parts.append(f"  {line}" if line else "")
        if declares_output:
            parts.extend(["  (void)output_size;", "  free(output);"])
        if declares_session:
            parts.append("  lammpc_session_destroy(session);")
        parts.extend(["  return 0;", "}", ""])

    parts.extend(
        [
            "int main(void) {",
            "  if (lammpc_version() == NULL)",
            "    return 1;",
        ]
    )
    for index in range(len(blocks)):
        parts.append(f"  (void)snippet_{index}();")
    parts.extend(["  return 0;", "}"])
    return "\n".join(parts) + "\n"


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "usage: extract_howto_c_snippets.py INPUT_MD OUTPUT_C",
            file=sys.stderr,
        )
        return 2
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    blocks = extract_c_blocks(input_path.read_text(encoding="utf-8"))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_c(input_path, blocks), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
