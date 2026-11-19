#!/usr/bin/env python3
"""Extract C code fences from the documentation into build-time tests.

Given docs/how-to.md and an output path, every fenced C block is wrapped in a
function with valid PotentialConfig and ForceInput buffers and its return
status is checked in main() (expecting 0 when LAMMPS is linked and 1 on the
stub build) so renamed symbols, changed parameter lists, stale struct fields,
or failing calls fail the test.

With `--check-docs ROOT` every C block in `README.md` and in the Markdown and
Org pages under `docs/` is compiled with `-Wall -Wextra -Werror` against the
public headers. Blocks in how-to.md are fragments and get the same wrapper as
above. Blocks in every other page must be complete translation units that
include their own headers. Compiling only, nothing is linked or run.
"""

from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

FRAGMENT_PAGES = frozenset({"how-to.md"})

SESSION_DECL_RE = re.compile(r"\bLammpcSession\s*\*\s*session\b")
OUTPUT_DECL_RE = re.compile(r"\bvoid\s*\*\s*output\b")
OPEN_FENCE_RE = re.compile(r"^(?P<fence>`{3,}|~{3,})[ \t]*(?P<info>.*)$")
C_INFO_RE = re.compile(
    r"^(?:[cC](?:$|[ \t].*|\{.*\}$)|\{\s*\.[cC](?:[ \t#][^}]*)?\})$"
)
ORG_BEGIN_RE = re.compile(
    r"^#\+begin_src[ \t]+(?P<lang>[^\s]+)(?:[ \t].*)?$", re.IGNORECASE
)
ORG_END_RE = re.compile(r"^#\+end_src$", re.IGNORECASE)


def _is_closing_fence(stripped: str, fence_marker: str) -> bool:
    if fence_marker == "#+end_src":
        return bool(ORG_END_RE.match(stripped))
    char = fence_marker[0]
    return len(stripped) >= len(fence_marker) and stripped == char * len(stripped)


def extract_c_blocks(
    markdown: str, require_blocks: bool = True
) -> list[tuple[int, str]]:
    blocks: list[tuple[int, str]] = []
    active_fence: str | None = None
    in_c_fence = False
    start_line = 0
    lines: list[str] = []

    for lineno, line in enumerate(markdown.splitlines(), start=1):
        stripped = line.strip()
        if active_fence is None:
            match = OPEN_FENCE_RE.match(stripped)
            if match is not None:
                active_fence = match.group("fence")
                info = match.group("info").strip()
                in_c_fence = bool(C_INFO_RE.match(info))
                if in_c_fence:
                    start_line = lineno + 1
                    lines = []
                continue
            org_match = ORG_BEGIN_RE.match(stripped)
            if org_match is not None:
                active_fence = "#+end_src"
                in_c_fence = org_match.group("lang").lower() == "c"
                if in_c_fence:
                    start_line = lineno + 1
                    lines = []
                continue
        else:
            if _is_closing_fence(stripped, active_fence):
                if in_c_fence:
                    blocks.append((start_line, "\n".join(lines)))
                active_fence = None
                in_c_fence = False
            elif in_c_fence:
                lines.append(line)

    if active_fence is not None and in_c_fence:
        raise RuntimeError(f"unclosed ```c fence starting at line {start_line - 1}")
    if not blocks and require_blocks:
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


def doc_pages(root: Path) -> list[Path]:
    pages = [root / "README.md"]
    for pattern in ("*.md", "*.org"):
        pages.extend(sorted((root / "docs").glob(pattern)))
    return pages


def compile_source(
    cc: list[str], include_dir: Path, source: str, label: str
) -> str | None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "snippet.c"
        path.write_text(source, encoding="utf-8")
        result = subprocess.run(
            [
                *cc,
                "-std=c11",
                "-Wall",
                "-Wextra",
                "-Werror",
                f"-I{include_dir}",
                "-c",
                str(path),
                "-o",
                str(Path(tmp) / "snippet.o"),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    if result.returncode == 0:
        return None
    return f"{label}\n{result.stderr}"


def check_docs(root: Path, cc: list[str], include_dir: Path) -> int:
    failures: list[str] = []
    checked = 0
    for page in doc_pages(root):
        label_path = page.relative_to(root).as_posix()
        blocks = extract_c_blocks(
            page.read_text(encoding="utf-8"), require_blocks=False
        )
        for start_line, body in blocks:
            checked += 1
            label = f"{label_path}:{start_line}"
            if page.name in FRAGMENT_PAGES:
                source = render_c(Path(label_path), [(start_line, body)])
            else:
                source = f'#line {start_line} "{label_path}"\n{body}\n'
            error = compile_source(cc, include_dir, source, label)
            print(f"{'FAIL' if error else 'ok  '} {label}")
            if error:
                failures.append(error)
    if checked == 0:
        print("no C snippets found in the documentation", file=sys.stderr)
        return 1
    for error in failures:
        print(error, file=sys.stderr)
    print(f"{checked - len(failures)} of {checked} documentation C snippets compile")
    return 1 if failures else 0


def run_self_test() -> int:
    sample = "\n".join(
        [
            "# Sample",
            "```c",
            "int plain = 1;",
            "```",
            "```c {.numberLines}",
            "int with_attr = 2;",
            "```",
            "```c{linenos=true}",
            "int attached_attr = 3;",
            "```",
            "```{.c}",
            "int pandoc_plain = 4;",
            "```",
            "``` {.c .numberLines #id}",
            "int pandoc_attr = 5;",
            "```",
            "~~~c",
            "int tilde_c = 6;",
            "~~~",
            "#+begin_src c",
            "int org_c = 7;",
            "#+end_src",
            "```sh",
            "echo not_c",
            "```",
            "```text",
            "plain text",
            "```",
            "```cpp",
            "int not_c_cpp = 0;",
            "```",
            "```c++",
            "int not_c_cxx = 0;",
            "```",
            "```c#",
            "int not_c_sharp = 0;",
            "```",
            "```cmake",
            "project(not_c)",
            "```",
            "```capnp",
            "@0x1234;",
            "```",
            "```{.cpp}",
            "int not_c_pandoc_cpp = 0;",
            "```",
            "````markdown",
            "```c",
            "int nested_inside_non_c = 0;",
            "```",
            "````",
            "```",
            "untagged block",
            "```",
            "#+begin_src sh",
            "echo org_sh",
            "#+end_src",
        ]
    )
    extracted = extract_c_blocks(sample)
    expected = [
        (3, "int plain = 1;"),
        (6, "int with_attr = 2;"),
        (9, "int attached_attr = 3;"),
        (12, "int pandoc_plain = 4;"),
        (15, "int pandoc_attr = 5;"),
        (18, "int tilde_c = 6;"),
        (21, "int org_c = 7;"),
    ]
    if extracted != expected:
        raise AssertionError(f"expected {expected!r}, got {extracted!r}")

    for non_c in (
        "```sh\necho hi\n```\n",
        "```cpp\nint x;\n```\n",
        "```c++\nint x;\n```\n",
        "```{.cpp}\nint x;\n```\n",
        "```\nint x;\n```\n",
    ):
        try:
            extract_c_blocks(non_c)
        except RuntimeError:
            pass
        else:
            raise AssertionError(f"non-C block unexpectedly extracted: {non_c!r}")

    try:
        extract_c_blocks("```{.c}\nint x = 0;\n")
    except RuntimeError:
        pass
    else:
        raise AssertionError("unclosed C fence did not raise RuntimeError")

    if extract_c_blocks("no code here\n", require_blocks=False):
        raise AssertionError("page without C fences produced blocks")

    return 0


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        return run_self_test()
    if len(sys.argv) >= 2 and sys.argv[1] == "--check-docs":
        args = sys.argv[2:]
        cc: list[str] = []
        include_dir: Path | None = None
        root: Path | None = None
        while args:
            flag = args.pop(0)
            if flag == "--cc" and args:
                cc.append(args.pop(0))
            elif flag == "--include-dir" and args:
                include_dir = Path(args.pop(0)).resolve()
            elif flag == "--root" and args:
                root = Path(args.pop(0)).resolve()
            else:
                break
        else:
            if cc and include_dir is not None and root is not None:
                return check_docs(root, cc, include_dir)
        print(
            "usage: extract_howto_c_snippets.py --check-docs --root DIR "
            "--include-dir DIR --cc COMPILER [--cc ARG ...]",
            file=sys.stderr,
        )
        return 2
    if len(sys.argv) != 3:
        print(
            "usage: extract_howto_c_snippets.py INPUT_MD OUTPUT_C | "
            "--self-test | --check-docs ...",
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
