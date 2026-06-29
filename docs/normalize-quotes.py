#!/usr/bin/env python3
"""
Replace straight ASCII quotes with typographic ("smart") quotes in all
Markdown files under docs/src/, leaving fenced code blocks and inline
math untouched.

Replacement rules
-----------------
"  →  " (U+201C, left double quotation mark)
        when preceded by whitespace or opening markup: space, tab, (, [, {, *, _
"  →  " (U+201D, right double quotation mark)
        in all other positions (after a word char or closing punctuation)
'  →  ' (U+2019, right single quotation mark / apostrophe)
        only when flanked by two word characters, i.e. contractions and
        possessives (don't, it's, ring's, …); all other apostrophes are left
        alone so as not to corrupt single-quoted strings in code snippets or
        prime notation (b', d') inside math blocks

Blocks that are skipped entirely
---------------------------------
- Fenced blocks: anything between opening ``` and closing ``` lines
- Inline math:   anything between unescaped $ … $ on the same line

In addition, any U+2019 that slipped into a math context (e.g. introduced
by a smart-quote editor) is reverted to ASCII ' so KaTeX can render it.

Run from the repo root:
    python3 docs/normalize-quotes.py
or via:
    make normalize-quotes
"""

import glob
import os
import re

# Characters that, when immediately preceding a ", indicate an *opening* quote.
# Notably excludes ,  ;  :  .  ?  !  )  ]  }  which precede *closing* quotes
# in patterns like  "word,"  "word."  "word?"  (trailing punct inside quotes).
OPENING_BEFORE = set(' \t\n([{*_')
RSQUOTE = '\u2019'   # U+2019 RIGHT SINGLE QUOTATION MARK


def smartify(path: str) -> bool:
    """Rewrite *path* in-place with smart quotes. Returns True if changed."""
    with open(path) as f:
        lines = f.readlines()

    result = []
    in_fence = False
    in_math_fence = False

    for line in lines:
        stripped = line.rstrip('\n')

        # Toggle fenced-block state on ``` lines (```math, ```julia, ``` …)
        if re.match(r'^\s*```', stripped):
            if not in_fence:
                in_math_fence = bool(re.match(r'^\s*```math\b', stripped))
            in_fence = not in_fence
            result.append(line)
            continue

        if in_fence:
            # Revert U+2019 in math fences back to ASCII '---smart-quote
            # editors can corrupt prime notation (b', c', d') and KaTeX
            # rejects U+2019 as an unknown symbol.
            if in_math_fence:
                result.append(line.replace(RSQUOTE, "'"))
            else:
                result.append(line)
            continue

        # Process prose character by character, toggling over inline $…$ math
        out = []
        in_math = False
        s = stripped

        for i, ch in enumerate(s):
            # Toggle inline-math state on unescaped $
            if ch == '$' and (i == 0 or s[i - 1] != '\\'):
                in_math = not in_math
                out.append(ch)

            elif in_math:
                # Revert any U+2019 that slipped into inline math back to ASCII.
                out.append("'" if ch == RSQUOTE else ch)

            elif ch == '"':
                prev = s[i - 1] if i > 0 else ' '
                out.append('“' if prev in OPENING_BEFORE else '”')

            elif ch == "'":
                prev = s[i - 1] if i > 0 else ''
                nxt  = s[i + 1] if i < len(s) - 1 else ''
                if prev.isalpha() and nxt.isalpha():
                    out.append(RSQUOTE)   # apostrophe / right single quote
                else:
                    out.append(ch)

            else:
                out.append(ch)

        result.append(''.join(out) + '\n')

    new_content = ''.join(result)
    with open(path) as f:
        old_content = f.read()

    if new_content != old_content:
        with open(path, 'w') as f:
            f.write(new_content)
        return True
    return False


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(__file__), 'src')
    changed = []
    for path in sorted(glob.glob(os.path.join(root, '**', '*.md'), recursive=True)):
        if smartify(path):
            changed.append(os.path.relpath(path))

    if changed:
        for p in changed:
            print(f'  normalized: {p}')
        print(f'\n{len(changed)} file(s) updated')
    else:
        print('all files already use smart quotes')
