#!/usr/bin/env python3
"""Extract C code fences from docs/how-to.md into a build-time smoke test.

Only fenced C blocks in the input Markdown file are extracted. Each snippet is
wrapped in a function with valid PotentialConfig and ForceInput buffers and its
return status is checked in main() (expecting 0 when LAMMPS is linked and 1 on
the stub build) so renamed symbols, changed parameter lists, stale struct
fields, or failing calls fail the test. Other doc pages (such as
docs/hard-boundary.org) are outside this check.
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


CONFIG_BYTES_C = """static const uint8_t kConfigBytes[] = {
    0x00, 0x00, 0x00, 0x00, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x02, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x09, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x21, 0x00, 0x00, 0x00, 0x32, 0x00, 0x00, 0x00, 0x21, 0x00, 0x00, 0x00,
    0x3a, 0x00, 0x00, 0x00, 0x21, 0x00, 0x00, 0x00, 0x5a, 0x00, 0x00, 0x00,
    0x25, 0x00, 0x00, 0x00, 0x0e, 0x00, 0x00, 0x00, 0x2d, 0x00, 0x00, 0x00,
    0x0c, 0x00, 0x00, 0x00, 0x2d, 0x00, 0x00, 0x00, 0x0d, 0x00, 0x00, 0x00,
    0x2d, 0x00, 0x00, 0x00, 0x32, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x6d, 0x65, 0x74, 0x61, 0x6c, 0x00, 0x00, 0x00, 0x61, 0x74, 0x6f, 0x6d,
    0x69, 0x63, 0x00, 0x00, 0x6c, 0x6a, 0x2f, 0x63, 0x75, 0x74, 0x20, 0x36,
    0x2e, 0x30, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0x62, 0x00, 0x00, 0x00, 0x31, 0x20, 0x31, 0x20, 0x30, 0x2e, 0x32, 0x20,
    0x32, 0x2e, 0x35, 0x00, 0x00, 0x00, 0x00, 0x00, 0x12, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xd3, 0x4d, 0x62, 0x10, 0x58, 0xf9, 0x43, 0x40,
    0x70, 0x20, 0x70, 0x20, 0x70, 0x00, 0x00, 0x00,
};"""

FORCE_INPUT_BYTES_C = """static const uint8_t kForceInputBytes[] = {
    0x00, 0x00, 0x00, 0x00, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x02, 0x00, 0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x11, 0x00, 0x00, 0x00,
    0x1d, 0x00, 0x00, 0x00, 0x19, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00,
    0x19, 0x00, 0x00, 0x00, 0x4d, 0x00, 0x00, 0x00, 0x39, 0x00, 0x00, 0x00,
    0x4a, 0x00, 0x00, 0x00, 0x3d, 0x00, 0x00, 0x00, 0x1a, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,
    0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x24, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x40,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x40, 0x61, 0x6e, 0x67, 0x73,
    0x74, 0x72, 0x6f, 0x6d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x65, 0x56, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
};"""


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
        CONFIG_BYTES_C,
        "",
        FORCE_INPUT_BYTES_C,
        "",
    ]

    snippet_takes_session: list[bool] = []
    for index, (start_line, body) in enumerate(blocks):
        declares_session = bool(SESSION_DECL_RE.search(body))
        declares_output = bool(OUTPUT_DECL_RE.search(body))
        takes_session = not declares_session
        snippet_takes_session.append(takes_session)
        param_sig = "LammpcSession *session" if takes_session else "void"
        parts.extend(
            [
                f"static int snippet_{index}({param_sig}) {{",
                "  const uint8_t *config_bytes = kConfigBytes;",
                "  const size_t config_size = sizeof(kConfigBytes);",
                "  const uint8_t *input = kForceInputBytes;",
                "  const size_t input_size = sizeof(kForceInputBytes);",
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
        if takes_session:
            parts.append("  (void)session;")
        parts.append(f'#line {start_line} "{source_label}"')
        for line in body.splitlines():
            parts.append(f"  {line}" if line else "")
        if declares_output:
            parts.extend(
                [
                    "  free(output);",
                    "  if (output_size == 0)",
                    "    return 1;",
                ]
            )
        if declares_session:
            parts.append("  lammpc_session_destroy(session);")
        parts.extend(["  return 0;", "}", ""])

    parts.extend(
        [
            "int main(void) {",
            "  if (lammpc_version() == NULL)",
            "    return 1;",
            "  const int available = lammpc_available();",
            "  const int expected = available ? 0 : 1;",
            "  LammpcSession *session = NULL;",
            "  if (available) {",
            "    session = lammpc_session_create_from_config(",
            "        kConfigBytes, sizeof(kConfigBytes));",
            "    if (session == NULL)",
            "      return 1;",
            "  }",
            "  int status = 0;",
        ]
    )
    for index, takes_session in enumerate(snippet_takes_session):
        call_arg = "session" if takes_session else ""
        parts.extend(
            [
                f"  if (snippet_{index}({call_arg}) != expected)",
                "    status = 1;",
            ]
        )
    parts.extend(
        [
            "  lammpc_session_destroy(session);",
            "  return status;",
            "}",
        ]
    )
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
