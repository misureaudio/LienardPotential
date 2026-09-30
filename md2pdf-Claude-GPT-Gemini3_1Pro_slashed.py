import sys
import subprocess
import re
import os
import tempfile
import argparse


# Enhanced Pandoc defaults with configurable options
def get_pandoc_defaults_yaml(orientation='portrait', paper_format='a4', margin='1in'):
    """Generate YAML defaults with configurable paper format and orientation"""

    # Paper format configurations - use geometry package format
    paper_configs = {
        'a4': 'a4paper',
        'a3': 'a3paper', 
        'letter': 'letterpaper',
        'legal': 'legalpaper',
        'a5': 'a5paper'
    }

    # Geometry settings based on orientation
    paper_setting = paper_configs.get(paper_format.lower(), 'a4paper')
    if orientation.lower() == 'landscape':
        geometry = f"{paper_setting}, landscape, margin={margin}"
    else:
        geometry = f"{paper_setting}, margin={margin}"

    return f"""
# Pandoc defaults file for md2pdf conversion
# ------------------------------------------
# Paper: {paper_format.upper()}, Orientation: {orientation.title()}

standalone: true
listings: true
pdf-engine: xelatex
highlight-style: tango

variables:
  mainfont: "DejaVu Serif"
  monofont: "DejaVu Sans Mono"
  geometry: "{geometry}"

  header-includes: |
    \\usepackage{{xcolor}}
    \\usepackage{{graphicx}}
    \\usepackage{{adjustbox}}
    \\usepackage{{listings}}
    \\usepackage{{longtable}}
    \\usepackage{{ltablex}}
    \\usepackage{{pifont}}
    \\usepackage{{slashed}}
    \\keepXColumns
    \\usepackage{{booktabs}}
    \\usepackage{{array}}
    % Better text wrapping and table handling
    \\usepackage{{ragged2e}}
    \\renewcommand{{\\arraystretch}}{{1.3}}
    % Configure listings for better text wrapping
    \\lstset{{
        language=Python,
        basicstyle=\\ttfamily\\small,
        breaklines=true,
        breakatwhitespace=true,
        showstringspaces=false,
        commentstyle=\\color{{green!60!black}},
        keywordstyle=\\color{{blue}},
        stringstyle=\\color{{purple}},
        postbreak=\\mbox{{\\textcolor{{red}}{{\\ensuremath{{\\hookrightarrow}}}}}},
        columns=flexible,
        keepspaces=true,
        tabsize=2,
        frame=single,
        framesep=5pt,
        xleftmargin=10pt,
        xrightmargin=10pt
    }}
    % Better table column types for text wrapping
    \\newcolumntype{{L}}[1]{{>{{\\RaggedRight\\arraybackslash}}p{{#1}}}}
    \\newcolumntype{{C}}[1]{{>{{\\centering\\arraybackslash}}p{{#1}}}}
    \\newcolumntype{{R}}[1]{{>{{\\RaggedLeft\\arraybackslash}}p{{#1}}}}
    % Allow tables to break across pages
    \\setlength{{\\LTpre}}{{0pt}}
    \\setlength{{\\LTpost}}{{0pt}}
"""

# BIG GEMINI MOD 4 TABLES
def preprocess_markdown_for_tables(input_md: str, orientation='portrait') -> str:
    """
    Finds Markdown pipe tables and converts them to LaTeX tabularx.
    Fixed to allow LaTeX math ($, ^, _, \) to pass through unescaped.
    """

    def escape_latex_text(text):
        """
        Escapes special LaTeX characters in table cells, preserving inline math spans.
        Splits the cell on $...$ regions and only escapes the plain-text portions.
        """
        # parts = re.split(r'(\$[^$]+\$)', text)
        parts = re.split(r'(\$[^$]+\$|\\\([^)]*\\\)|\\\[[^\]]*\\\])', text) # GPT5.x
        result = []
        for part in parts:
            '''
            if part.startswith('$') and part.endswith('$') and len(part) > 1:
                # Math span — pass through, only escape table-breaking &
                inner = part[1:-1].replace('&', r'\&')
                result.append('$' + inner + '$')
            '''
            if ((part.startswith('$') and part.endswith('$')) or
                (part.startswith(r'\(') and part.endswith(r'\)')) or
                    (part.startswith(r'\[') and part.endswith(r'\]'))):
                inner = part
                inner = inner.replace('&', r'\&')
                # result.append(inner)
                result.append('$' + inner + '$')
            else:
                # Plain text — full LaTeX escaping
                replacements = [
                    ('\\', r'\textbackslash{}'),
                    ('&', r'\&'),
                    ('%', r'\%'),
                    ('$', r'\$'),
                    ('#', r'\#'),
                    ('_', r'\_'),
                    ('{', r'\{'),
                    ('}', r'\}'),
                    ('~', r'\textasciitilde{}'),
                    ('^', r'\textasciicircum{}'),
                ]
                for old, new in replacements:
                    part = part.replace(old, new)
                result.append(part)

        text = ''.join(result)
        # Handle Markdown **bold** -> \textbf{bold}
        text = re.sub(r'\*\*(.*?)\*\*', r'\\textbf{\1}', text)
        return text.strip()

    # Regex to find the table
    table_pattern = r'(?:^\|(?:.*\|)+.*\n)(?:^\|(?: *[:–-]+ *\|)+.*\n)(?:(?:^\|(?:.*\|)+.*\n?)*)'
    # table_pattern = r'(^\|.*\|\n^\|[ :\-|]+\|\n(?:^\|.*\|\n?)*)'  #G PT5.x

    def convert_table_to_latex(match):
        table_md = match.group(0).strip()
        lines = table_md.split('\n')

        if len(lines) < 2:
            return table_md

        header_line = lines[0]
        row_lines = lines[2:]

        def split_row(row_md):
            """Split a markdown table row on | but not inside $...$ math spans."""
            cells = []
            current = ''
            in_math = False
            for ch in row_md.strip().strip('|'):
                if ch == '$':
                    in_math = not in_math
                    current += ch
                elif ch == '|' and not in_math:
                    cells.append(current.strip())
                    current = ''
                else:
                    current += ch
            if current.strip():
                cells.append(current.strip())
            return cells

        header_cells_md = split_row(header_line)
        num_cols = len(header_cells_md)
        if num_cols == 0:
            return table_md

        col_spec = "|".join(["X"] * num_cols)

        latex_table = [
            r'\begin{tabularx}{\linewidth}{|' + col_spec + r'|}',
            r'\hline'
        ]

        def render_cell(text):
            """
            Render a table cell for tabularx X columns.
            Math-heavy cells (>50% math content) use adjustbox{max width=\hsize}
            which correctly respects the column width, not page width.
            Mixed/plain cells with inline math use escape_latex_text which
            preserves $...$ spans as-is.
            """
            stripped = text.strip()
            math_spans = re.findall(r'\$[^$]+\$', stripped)
            math_chars = sum(len(s) for s in math_spans)
            if math_spans and math_chars > 0.5 * max(len(stripped), 1):
                # Math-heavy cell: strip outer $ delimiters and wrap entire
                # expression in adjustbox so it scales to column width
                inner = stripped.replace('$', ' ').strip()
                return r'\adjustbox{max width=\hsize}{$\displaystyle ' + inner + r'$}'
            else:
                # Mixed or plain: escape text, preserve inline math spans
                return escape_latex_text(stripped)

        # Header
        header_cells_tex = [r'\textbf{' + escape_latex_text(cell) + '}' for cell in header_cells_md]
        latex_table.append(' & '.join(header_cells_tex) + r' \\ \hline')

        # Rows
        for row_md in row_lines:
            if not row_md.strip(): continue
            row_cells_md = split_row(row_md)
            while len(row_cells_md) < num_cols:
                row_cells_md.append('')
            row_cells_md = row_cells_md[:num_cols]
            row_cells_tex = [render_cell(cell) for cell in row_cells_md]
            latex_table.append(' & '.join(row_cells_tex) + r' \\ \hline')

        latex_table.append(r'\end{tabularx}')

        # CHANGED: Wrap the result in a raw LaTeX block for Pandoc
        return '\n\n```{=latex}\n' + '\n'.join(latex_table) + '\n```\n\n'

    return re.sub(table_pattern, convert_table_to_latex, input_md, flags=re.MULTILINE)


def md2pdf(input_md: str, orientation='portrait', paper_format='a4', margin='1in') -> bytes:
    """
    Convert Markdown to PDF with configurable orientation and paper format.

    Args:
        input_md: Markdown content as string
        orientation: 'portrait' or 'landscape' 
        paper_format: 'a4', 'a3', 'letter', 'legal', 'a5'
        margin: Margin specification (e.g., '1in', '2cm', '0.5in')
    """

    # --- Replace LaTeX-style math delimiters with TeX-style ones
    def math_delim_to_dollars(latex_style):
        latex_style = re.sub(
            r'\\\((.*?)\\\)',
            r' $\1$ ',
            latex_style,
            flags=re.DOTALL
        )
        latex_style = re.sub(
            r'\\\[(.*?)\\\]',
            r' $$\1$$ ',
            latex_style,
            flags=re.DOTALL
        )
        return re.sub(
            r'\n\s*\\[\[\]]\s*\n',
            '\n$$\n',
            latex_style,
            flags=re.DOTALL
        )

    input_md = math_delim_to_dollars(input_md)

    # --- Enhanced Preprocessing ---
    # Fix Unicode escape sequences
    input_md = re.sub(
        r'\\[uU]([0-9a-fA-F]{4})',
        lambda m: chr(int(m.group(1), 16)),
        input_md
    )

    # Fix double subscripts in math: e.g. v_h_x -> v_{h,x}
    # Operates on both inline $...$ and display $$...$$ math
    def fix_double_subscripts(math_content):
        """
        Detect and fix double subscript patterns like a_b_c -> a_{b,c}
        Handles simple tokens: letters, digits, braced groups.
        """
        # Pattern: something_token_token where tokens are simple (not already braced)
        # e.g. v_h_x -> v_{h,x}  u_h^{n}_k -> leave alone (superscript between)
        # Match: base _simple _simple  where simple = single char or {group}
        def merge_subs(m):
            base = m.group(1)
            sub1 = m.group(2)
            sub2 = m.group(3)
            # strip braces if present
            s1 = sub1.strip('{}')
            s2 = sub2.strip('{}')
            return f'{base}_{{{s1},{s2}}}'

        # Only fix bare double subscripts (no superscript between them)
        fixed = re.sub(
            r'([A-Za-z0-9}])_([A-Za-z0-9]|\{[^}]+\})_([A-Za-z0-9]|\{[^}]+\})',
            merge_subs,
            math_content
        )
        return fixed

    def fix_unbraced_command_subscripts(math_content):
        """
        Fix _\command{arg} -> _{\command{arg}} and similarly for ^.
        Handles multi-brace commands like \frac{a}{b}.
        """
        result = []
        i = 0
        s = math_content
        n = len(s)
        trigger = re.compile(r'([_^])\\([A-Za-z]+)')
        while i < n:
            m = trigger.search(s, i)
            if not m:
                result.append(s[i:])
                break
            result.append(s[i:m.start()])
            op = m.group(1)
            cmd = m.group(2)
            j = m.end()
            args = ''
            while j < n and s[j] == '{':
                depth = 0
                k = j
                while k < n:
                    if s[k] == '{':
                        depth += 1
                    elif s[k] == '}':
                        depth -= 1
                        if depth == 0:
                            k += 1
                            break
                    k += 1
                args += s[j:k]
                j = k
            if args:
                result.append(f'{op}{{\\{cmd}{args}}}')
            else:
                result.append(f'{op}{{\\{cmd}}}')
            i = j
        return ''.join(result)

    def apply_fix_to_math_spans(text):
        def fix_all(math_content):
            math_content = fix_double_subscripts(math_content)
            math_content = fix_unbraced_command_subscripts(math_content)
            return math_content
        text = re.sub(
            r'\$\$(.+?)\$\$',
            lambda m: '$$' + fix_all(m.group(1)) + '$$',
            text, flags=re.DOTALL
        )
        text = re.sub(
            r'(?<!\$)\$([^$]+?)\$(?!\$)',
            lambda m: '$' + fix_all(m.group(1)) + '$',
            text
        )
        return text

    input_md = apply_fix_to_math_spans(input_md)

    # Remove newlines between $$ and LaTeX formulas
    input_md = re.sub(
        r'\$\$([^$]+)\$\$',
        lambda m: '$$' + m.group(1).replace("\n", " ") + '$$',
        input_md
    )

    '''
    # Fix inline math and handle pipe characters in tables
    def fix_math_in_tables(match: re.Match) -> str:
        content = match.group(1).strip()
        fixed_content = re.sub(r'(?<!\\)\|', r' \\vert ', content)
        return f'${fixed_content}$'

    input_md = re.sub(
        r'(?<!\$)\$ *([^$]*?) *?\$(?!\$)',
        fix_math_in_tables,
        input_md
    )
    '''

    # Escape dollar signs in text
    '''
    pattern = r"""
        (?<!\\)(
            \$\$(?!\s).*?(?<!\s)\$\$
        )|(?<!\\)(
            \$(?!\s)(?!\d).*?(?<!\s)\$
        )|(?<!\\)(\$)
    """
    '''
    # Gemini correction
    pattern = r"""
        (?<!\\)(
            \$\$(?!\s).*?(?<!\s)\$\$
        )|(?<!\\)(
            \$(?!\s).*?(?<!\s)\$
        )|(?<!\\)(\$)
    """

    def replacement_logic(match):
        if match.group(1): return match.group(1)
        if match.group(2): return match.group(2)
        if match.group(3): return r'\$'
        return match.group(0)

    # Table preprocessing MUST run before dollar escape to protect
    # inline math inside table cells from being mangled by the escaper
    input_md = preprocess_markdown_for_tables(input_md, orientation)

    input_md = re.sub(pattern, replacement_logic, input_md, flags=re.VERBOSE | re.DOTALL)
    input_md = input_md.replace('```<end_code>', '```')

    # --- Pandoc Execution ---
    temp_defaults_path = None
    try:
        # Create temporary defaults file with specified configuration
        defaults_content = get_pandoc_defaults_yaml(orientation, paper_format, margin)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.yaml', encoding='utf-8') as defaults_file:
            defaults_file.write(defaults_content)
            temp_defaults_path = defaults_file.name
        # BIG GEMINI MOD 4 TABLES
        # markdowntype = 'commonmark+task_lists+strikeout+pipe_tables+autolink_bare_uris+tex_math_dollars'
        markdowntype = 'commonmark+task_lists+strikeout+pipe_tables+autolink_bare_uris+tex_math_dollars+raw_attribute'

        command = [
            "pandoc",
            f"--defaults={temp_defaults_path}",
            "-f", markdowntype,
            "-t", "pdf",
            "-o", "-"
        ]

        process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout_bytes, stderr_bytes = process.communicate(input=input_md.encode('utf-8'))

    finally:
        if temp_defaults_path and os.path.exists(temp_defaults_path):
            os.remove(temp_defaults_path)

    # Handle Pandoc output
    if process.returncode != 0:
        raise RuntimeError(
            f"Pandoc failed with return code {process.returncode}. "
            f"Error: {stderr_bytes.decode('utf-8', errors='replace')}"
        )

    return stdout_bytes


def main():
    """Enhanced command-line interface with orientation and paper format options"""
    parser = argparse.ArgumentParser(
        description='Convert Markdown to PDF with configurable layout options',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python md2pdf.py input.md                           # Default: A4 portrait
  python md2pdf.py input.md -o landscape              # A4 landscape  
  python md2pdf.py input.md -p a3 -o landscape        # A3 landscape
  python md2pdf.py input.md -p letter -m 0.5in        # US Letter, small margins
  python md2pdf.py - -o landscape < input.md          # Read from stdin, landscape
        """
    )

    parser.add_argument('input', nargs='?', default='-',
                        help='Input markdown file (use "-" or omit for stdin)')
    parser.add_argument('-o', '--orientation', choices=['portrait', 'landscape'], 
                        default='portrait', help='Page orientation (default: portrait)')
    parser.add_argument('-p', '--paper', choices=['a4', 'a3', 'letter', 'legal', 'a5'],
                        default='a4', help='Paper format (default: a4)')
    parser.add_argument('-m', '--margin', default='1in',
                        help='Page margins (default: 1in, e.g., "0.5in", "2cm")')
    parser.add_argument('--output', '-w', 
                        help='Output PDF file (default: stdout)')

    args = parser.parse_args()

    # Read input
    if args.input == '-':
        input_data = sys.stdin.read()
    else:
        try:
            with open(args.input, 'rt', encoding='utf-8') as f:
                input_data = f.read()
        except FileNotFoundError:
            print(f'*** Error: File "{args.input}" not found', file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f'*** Error reading file: {e}', file=sys.stderr)
            sys.exit(1)

    # Convert to PDF
    try:
        output_pdf_bytes = md2pdf(input_data, args.orientation, args.paper, args.margin)
    except Exception as e:
        print(f'*** Error during conversion: {e}', file=sys.stderr)
        sys.exit(1)

    # Write output
    if args.output:
        try:
            with open(args.output, 'wb') as f:
                f.write(output_pdf_bytes)
            print(f'PDF written to {args.output}', file=sys.stderr)
        except Exception as e:
            print(f'*** Error writing output file: {e}', file=sys.stderr)
            sys.exit(1)
    else:
        sys.stdout.buffer.write(output_pdf_bytes)
        sys.stdout.buffer.flush()

    sys.exit(0)


if __name__ == '__main__':
    main()
