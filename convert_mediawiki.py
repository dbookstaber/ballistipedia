#!/usr/bin/env python3
"""
Convert MediaWiki XML export to MkDocs-compatible Markdown.
Handles: math, refs/footnotes, images, media links, tables, wikilinks, formatting.
"""

import xml.etree.ElementTree as ET
import re
import os
import sys
import hashlib
import urllib.request
import urllib.parse
import time

NS = '{http://www.mediawiki.org/xml/export-0.10/}'

# Pages to skip (redirects and meta pages handled separately)
SKIP_TITLES = {'Main Page'}


def slugify(title):
    """Convert a wiki title to a filename slug."""
    slug = title.replace(' ', '-').replace('?', '').replace("'", '').replace('"', '')
    slug = slug.replace('ρ', 'rho')
    slug = re.sub(r'[^\w\-.]', '-', slug)
    slug = re.sub(r'-+', '-', slug).strip('-')
    return slug.lower()


def parse_xml(xml_path):
    """Parse the MediaWiki XML export and return pages and redirects."""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    pages = {}
    redirects = {}

    for page in root.findall(f'{NS}page'):
        title = page.find(f'{NS}title').text
        ns = page.find(f'{NS}ns').text
        if ns != '0':
            continue

        redirect = page.find(f'{NS}redirect')
        if redirect is not None:
            redirects[title] = redirect.get('title')
            continue

        if title in SKIP_TITLES:
            continue

        revision = page.find(f'{NS}revision')
        text_elem = revision.find(f'{NS}text')
        text = text_elem.text if text_elem is not None and text_elem.text else ''
        if text.startswith('#REDIRECT'):
            target = re.search(r'\[\[(.+?)\]\]', text)
            if target:
                redirects[title] = target.group(1)
            continue

        pages[title] = text

    return pages, redirects


def make_section_slug(section):
    """Convert a MediaWiki section name to a Markdown-compatible anchor slug."""
    # First decode MediaWiki-style URL encoding (.XX hex codes)
    section = re.sub(r'\.([0-9A-Fa-f]{2})', lambda m: chr(int(m.group(1), 16)), section)
    # Also decode percent-encoding
    section = urllib.parse.unquote(section)
    sec_slug = section.lower().replace(' ', '-').replace('_', '-')
    # Strip everything except alphanumeric, hyphens
    sec_slug = re.sub(r'[^a-z0-9\-]', '', sec_slug)
    # Collapse multiple hyphens and strip leading/trailing hyphens
    sec_slug = re.sub(r'-{2,}', '-', sec_slug)
    sec_slug = sec_slug.strip('-')
    return sec_slug


def resolve_link(title, redirects, pages):
    """Resolve a wiki title to its final slug, following redirects."""
    seen = set()
    current = title
    while current in redirects and current not in seen:
        seen.add(current)
        current = redirects[current]
    if current in pages:
        if current == 'Home':
            return 'index'
        return slugify(current)
    return None


def convert_wikitext(title, text, redirects, pages, all_files, all_media):
    """Convert MediaWiki wikitext to MkDocs Markdown."""
    # Track footnotes
    footnotes = []
    footnote_counter = [0]
    named_refs = {}  # name -> footnote number

    def replace_ref(match):
        full_tag = match.group(0)
        content = match.group(1) if match.group(1) else ''

        # Check for named ref
        name_match = re.search(r'name\s*=\s*"?([^"\s>]+)"?', full_tag)
        ref_name = name_match.group(1) if name_match else None

        # Check if self-closing (ref with name but no content, referencing earlier ref)
        is_self_closing = content.strip() == '' or full_tag.rstrip().endswith('/>')

        if ref_name and ref_name in named_refs and (is_self_closing or not content.strip()):
            # Back-reference to existing named ref
            return f'[^{named_refs[ref_name]}]'

        footnote_counter[0] += 1
        n = footnote_counter[0]

        if ref_name:
            named_refs[ref_name] = n

        # Convert wiki links inside ref
        ref_md = convert_inline_links(content, redirects, pages, all_files, all_media)
        footnotes.append((n, ref_md))
        return f'[^{n}]'

    def convert_inline_links(text, redirects, pages, all_files, all_media):
        """Convert [[...]] and [...] links inline."""
        # Category links (display) - these are category pages, link as best we can
        def category_display_link(m):
            cat_name = m.group(1).strip()
            display = m.group(2)
            # Check if there's a page with this name
            slug = resolve_link(cat_name, redirects, pages)
            if slug:
                return f'[{display}]({slug}.md)'
            # Otherwise just display as text
            return display
        text = re.sub(
            r'\[\[:Category:([^\]|]+)\|([^\]]+)\]\]',
            category_display_link,
            text
        )
        # Category tags (strip)
        text = re.sub(r'\[\[Category:[^\]]+\]\]', '', text)

        # Media links [[Media:...]]
        def media_link(m):
            fname = m.group(1).strip()
            display = m.group(2) if m.group(2) else fname
            safe_fname = fname.replace(' ', '%20')
            all_media.add(fname)
            return f'[{display}](media/{safe_fname})'

        text = re.sub(r'\[\[Media:([^\]|]+?)(?:\|([^\]]+))?\]\]', media_link, text)

        # File links [[File:...]] are handled separately
        def file_link(m):
            return convert_file_tag(m.group(0), all_files)
        text = re.sub(r'\[\[File:[^\]]+\]\]', file_link, text)

        # Internal links with section [[Page#Section|Display]]
        def internal_link(m):
            full = m.group(1)
            display = m.group(2) if m.group(2) else None

            # Parse page and section
            if '#' in full:
                page_part, section = full.split('#', 1)
            else:
                page_part, section = full, None

            page_part = page_part.replace('_', ' ').strip()

            if not page_part and section:
                # Same-page section link
                sec_slug = make_section_slug(section)
                display_text = display or section.replace('_', ' ')
                return f'[{display_text}](#{sec_slug})'

            slug = resolve_link(page_part, redirects, pages)
            if slug is None:
                # External or broken link - just display text
                return display or page_part

            display_text = display or page_part
            if section:
                sec_slug = make_section_slug(section)
                return f'[{display_text}]({slug}.md#{sec_slug})'
            return f'[{display_text}]({slug}.md)'

        text = re.sub(r'\[\[([^\]|]+?)(?:\|([^\]]+?))?\]\]', internal_link, text)

        # External links [url display]
        text = re.sub(
            r'\[([a-z]+://[^\s\]]+)\s+([^\]]+)\]',
            r'[\2](\1)',
            text
        )
        # Bare external links [url]
        text = re.sub(r'\[([a-z]+://[^\s\]]+)\]', r'<\1>', text)

        return text

    def convert_file_tag(tag, all_files):
        """Convert [[File:...]] to markdown image."""
        inner = tag[7:-2]  # Strip [[File: and ]]
        parts = [p.strip() for p in inner.split('|')]
        filename = parts[0].strip()
        # Guard against newlines or extraneous whitespace in filename
        filename = filename.split('\n')[0].strip()
        safe_fname = filename.replace(' ', '%20')

        # Check if this is a non-image file (source code, spreadsheets, etc.)
        ext = os.path.splitext(filename)[1].lower()
        image_exts = {'.png', '.jpg', '.jpeg', '.gif', '.svg', '.bmp', '.webp'}

        if ext not in image_exts:
            # Treat as a downloadable media file, not an image
            all_media.add(filename)
            return f'[{filename}](media/{safe_fname})'
        all_files.add(filename)

        safe_fname = filename.replace(' ', '%20')

        # Parse options
        caption = ''
        width = ''
        is_thumb = False
        align = ''

        for p in parts[1:]:
            if p.endswith('px'):
                width = p
            elif p in ('thumb', 'thumbnail', 'frame'):
                is_thumb = True
            elif p in ('right', 'left', 'center', 'none'):
                align = p
            else:
                # Last non-option part is typically the caption
                caption = p

        # Clean caption (may have wiki markup)
        if caption:
            caption = re.sub(r"<math>(.+?)</math>", r'$\1$', caption)
            caption = caption.replace("'''", '**').replace("''", '*')

        alt = caption or filename
        alt = re.sub(r'[<>]', '', alt)

        if is_thumb or caption:
            # Use figure with caption
            result = f'\n![{alt}](images/{safe_fname})'
            if caption:
                result += f'\n{{ .figure-caption }}\n*{caption}*'
            result += '\n'
            return result
        else:
            return f'![{alt}](images/{safe_fname})'

    # --- Main conversion ---

    # Protect File/Media links from entity decoding that might break parsing
    # Temporarily replace them with placeholders
    file_placeholders = []
    def save_file_link(m):
        idx = len(file_placeholders)
        file_placeholders.append(m.group(0))
        return f'__FILE_PLACEHOLDER_{idx}__'
    text = re.sub(r'\[\[(?:File|Media):[^\]]+\]\]', save_file_link, text)

    # Decode HTML entities
    text = text.replace('&#039;', "'")
    text = text.replace('&quot;', '"')
    text = text.replace('&amp;', '&')
    text = text.replace('&lt;', '<')
    text = text.replace('&gt;', '>')
    text = text.replace('&mdash;', '—')
    text = text.replace('&ndash;', '–')
    text = text.replace('&nbsp;', ' ')

    # Restore File/Media links with entities decoded only within them
    for idx, original in enumerate(file_placeholders):
        decoded = original.replace('&#039;', "'").replace('&quot;', '"')
        decoded = decoded.replace('&amp;', '&').replace('&lt;', '<').replace('&gt;', '>')
        text = text.replace(f'__FILE_PLACEHOLDER_{idx}__', decoded)

    # Remove navigation boilerplate at top and bottom
    text = re.sub(r'<p style="text-align:right"><B>Previous:</B>.*?</p>', '', text)
    text = re.sub(r'<p style="text-align:right"><B>Next:</B>.*?</p>', '', text)

    # Remove trailing <BR/><HR/> etc
    text = re.sub(r'(?:<BR\s*/?>|<br\s*/?>|<HR\s*/?>|<hr\s*/?>)\s*$', '', text, flags=re.IGNORECASE)

    # Handle <ref>...</ref> -> footnotes (including named refs)
    text = re.sub(r'<ref[^>]*>(.*?)</ref>', replace_ref, text, flags=re.DOTALL)
    # Handle self-closing named refs <ref name="..." />
    text = re.sub(r'<ref\s+name\s*=\s*"?([^"\s/>]+)"?\s*/>', lambda m: f'[^{named_refs.get(m.group(1), "?")}]', text)
    # Handle <references /> -> will be replaced with footnote list
    text = re.sub(r'<references\s*/?>', '<!-- FOOTNOTES -->', text)

    # Convert math tags
    # Block math: only when <math>...</math> is the SOLE content on the line (possibly with : indent)
    def block_math(m):
        latex = m.group(1).strip()
        return f'\n$$\n{latex}\n$$\n'

    text = re.sub(r'^:?[ ]*<math>((?:(?!</math>).)*)</math>[ ]*$', block_math, text, flags=re.MULTILINE)

    # Inline math (also handle multiline math tags)
    def inline_math_replace(m):
        content = m.group(1).replace('\n', ' ').strip()
        return f'${content}$'
    text = re.sub(r'<math>(.*?)</math>', inline_math_replace, text, flags=re.DOTALL)

    # Convert wiki-style lists BEFORE header conversion (since both use # character)
    # Wiki lists: # for ordered, * for unordered, #* for mixed nesting
    lines = text.split('\n')
    converted_lines = []
    for line in lines:
        stripped = line.strip()
        # Sub-items: #* or ## or ** 
        if re.match(r'^#\*\s*', stripped):
            converted_lines.append('    - ' + re.sub(r'^#\*\s*', '', stripped))
            continue
        if re.match(r'^##\s*', stripped):
            converted_lines.append('    1. ' + re.sub(r'^##\s*', '', stripped))
            continue
        if re.match(r'^\*\*\s*', stripped):
            converted_lines.append('    - ' + re.sub(r'^\*\*\s*', '', stripped))
            continue
        # Top-level ordered list (# at start of line, not part of == header syntax)
        if re.match(r'^#\s*', stripped) and not re.match(r'^#\s*=', stripped):
            converted_lines.append('1. ' + re.sub(r'^#\s*', '', stripped))
            continue
        # Top-level unordered list
        if re.match(r'^\*\s*', stripped):
            converted_lines.append('- ' + re.sub(r'^\*\s*', '', stripped))
            continue
        converted_lines.append(line)
    text = '\n'.join(converted_lines)

    # Convert section headers (MediaWiki uses = for headers)
    # First handle headers that contain wiki links
    def convert_header(m):
        level = len(m.group(1))
        content = m.group(2).strip()
        prefix = '#' * level + ' '
        return prefix + content

    text = re.sub(r'^(={1,4})\s*(.+?)\s*={1,4}', convert_header, text, flags=re.MULTILINE)

    # Convert bold/italic (wiki-style)
    text = re.sub(r"'''(.+?)'''", r'**\1**', text)
    text = re.sub(r"''(.+?)''", r'*\1*', text)

    # Strip <small> and </small> tags
    text = re.sub(r'</?small>', '', text, flags=re.IGNORECASE)

    # Convert HTML formatting
    text = re.sub(r'<code>(.*?)</code>', r'`\1`', text, flags=re.DOTALL)
    text = re.sub(r'<sub>(.*?)</sub>', r'~\1~', text)  # MkDocs Material supports this
    text = re.sub(r'<sup>(.*?)</sup>', r'^\1^', text)
    text = re.sub(r'<[Bb]>(.*?)</[Bb]>', r'**\1**', text)
    text = re.sub(r'<[Ii]>(.*?)</[Ii]>', r'*\1*', text)
    text = re.sub(r'<blockquote>(.*?)</blockquote>', r'> \1', text, flags=re.DOTALL)

    # Clean up HTML tags - but preserve <div id="..."> as markdown anchor
    text = re.sub(r'<br\s*/?>', '  \n', text, flags=re.IGNORECASE)
    text = re.sub(r'<[Hh][Rr]\s*/?>', '\n---\n', text)
    text = re.sub(r'<center>(.*?)</center>', r'\1', text, flags=re.DOTALL)
    # Preserve div id anchors as markdown-compatible anchors (lowercase for consistency)
    text = re.sub(r'<div\s+id="([^"]+)"\s*>\s*</div>', lambda m: f'<a id="{m.group(1).lower()}"></a>', text)
    text = re.sub(r'<div[^>]*>(.*?)</div>', r'\1', text, flags=re.DOTALL)

    # Remove TOC markers
    text = re.sub(r'\{?\|align=right\s*\|__TOC__\s*\|\}?', '', text)
    text = re.sub(r'__TOC__', '', text)
    text = re.sub(r'__NOTOC__', '', text)

    # Handle draft notice tables
    text = re.sub(
        r'\{?\|\s*class="wikitable"\s*\|\s*\[\[File:Bullseye\.jpg\|50px\]\]\s*(.+?)\s*\|\}?',
        r'!!! warning\n    \1\n',
        text,
        flags=re.DOTALL
    )

    # Convert wiki tables to markdown tables
    text = convert_tables(text)

    # Convert links
    text = convert_inline_links(text, redirects, pages, all_files, all_media)

    # Handle indented lines (: prefix)
    lines = text.split('\n')
    converted_lines = []

    for line in lines:
        stripped = line.strip()

        # Indented lines (: prefix) that aren't already converted
        depth = 0
        temp = stripped
        while temp.startswith(':'):
            depth += 1
            temp = temp[1:]
        if depth > 0:
            temp = temp.strip()
            if temp.startswith('$$') or temp.startswith('$'):
                converted_lines.append(temp)
            else:
                # Use blockquote for indentation
                converted_lines.append('> ' * depth + temp)
            continue

        converted_lines.append(line)

    text = '\n'.join(converted_lines)

    # Add footnotes at the end
    if footnotes:
        # Replace the footnotes placeholder or append
        footnote_text = '\n\n'
        for n, content in footnotes:
            # Clean up the content
            content = content.strip()
            footnote_text += f'[^{n}]: {content}\n\n'

        if '<!-- FOOTNOTES -->' in text:
            text = text.replace('<!-- FOOTNOTES -->', footnote_text)
        else:
            text += footnote_text

    # Remove leftover <!-- FOOTNOTES --> if no footnotes
    text = text.replace('<!-- FOOTNOTES -->', '')

    # Ensure blank line before list starts (Python-Markdown requires this)
    # Only insert when the previous line is NOT a list item and NOT blank
    text = re.sub(r'(\n[^\n\-\*\d >][^\n]*)\n((?:[-*]|\d+\.) )', r'\1\n\n\2', text)

    # Clean up excessive blank lines
    text = re.sub(r'\n{4,}', '\n\n\n', text)

    # Clean up leading/trailing whitespace
    text = text.strip() + '\n'

    return text


def convert_tables(text):
    """Convert MediaWiki tables to Markdown tables."""
    result = []
    lines = text.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith('{|'):
            # Start of wiki table
            table_lines = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('|}'):
                table_lines.append(lines[i])
                i += 1
            # i now points to |} line
            md_table = wiki_table_to_md(table_lines)
            result.append(md_table)
        else:
            result.append(line)
        i += 1
    return '\n'.join(result)


def wiki_table_to_md(lines):
    """Convert wiki table lines to markdown table."""
    rows = []
    current_row = []
    has_header = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('|-'):
            if current_row:
                rows.append(current_row)
                current_row = []
        elif stripped.startswith('!'):
            has_header = True
            # Header cells
            cells = re.split(r'\|\||!!', stripped[1:])
            current_row = [c.strip() for c in cells]
        elif stripped.startswith('|'):
            cells = re.split(r'\|\|', stripped[1:])
            current_row.extend([c.strip() for c in cells])

    if current_row:
        rows.append(current_row)

    if not rows:
        return ''

    # Build markdown table
    if has_header and len(rows) > 0:
        header = rows[0]
        data_rows = rows[1:]
    else:
        header = ['' for _ in rows[0]] if rows else []
        data_rows = rows

    ncols = max(len(r) for r in [header] + data_rows) if rows else 0
    if ncols == 0:
        return ''

    # Pad rows
    header = header + [''] * (ncols - len(header))
    data_rows = [r + [''] * (ncols - len(r)) for r in data_rows]

    md = '| ' + ' | '.join(header) + ' |\n'
    md += '| ' + ' | '.join(['---'] * ncols) + ' |\n'
    for row in data_rows:
        md += '| ' + ' | '.join(row) + ' |\n'

    return md


def get_mediawiki_image_url(filename, base_url='http://ballistipedia.com'):
    """Compute the URL for a MediaWiki uploaded file using the MD5 hash scheme."""
    # MediaWiki stores files at /images/X/XY/Filename
    # where X = first char of MD5, XY = first two chars of MD5 of the filename
    name_for_hash = filename.replace(' ', '_')
    md5 = hashlib.md5(name_for_hash.encode('utf-8')).hexdigest()
    encoded_name = urllib.parse.quote(name_for_hash)
    return f'{base_url}/images/{md5[0]}/{md5[:2]}/{encoded_name}'


def download_file(url, dest_path, retries=2):
    """Download a file from a URL to a local path."""
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=30) as response:
                data = response.read()
                os.makedirs(os.path.dirname(dest_path), exist_ok=True)
                with open(dest_path, 'wb') as f:
                    f.write(data)
                return True
        except Exception as e:
            if attempt < retries:
                time.sleep(1)
            else:
                print(f'  FAILED to download {url}: {e}')
                return False
    return False


def download_all_files(files, media_files, output_dir, base_url='http://ballistipedia.com'):
    """Download all referenced images and media files from the wiki."""
    images_dir = os.path.join(output_dir, 'images')
    media_dir = os.path.join(output_dir, 'media')
    os.makedirs(images_dir, exist_ok=True)
    os.makedirs(media_dir, exist_ok=True)

    all_downloads = []

    for fname in sorted(files):
        url = get_mediawiki_image_url(fname, base_url)
        dest = os.path.join(images_dir, fname)
        all_downloads.append(('image', fname, url, dest))

    for fname in sorted(media_files):
        url = get_mediawiki_image_url(fname, base_url)
        dest = os.path.join(media_dir, fname)
        all_downloads.append(('media', fname, url, dest))

    print(f'\nDownloading {len(all_downloads)} files from {base_url}...')
    success = 0
    failed = []

    for kind, fname, url, dest in all_downloads:
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f'  [SKIP] {kind}: {fname} (already exists)')
            success += 1
            continue
        print(f'  [{kind}] {fname}...', end=' ')
        if download_file(url, dest):
            size = os.path.getsize(dest)
            print(f'OK ({size} bytes)')
            success += 1
        else:
            failed.append((kind, fname, url))

    print(f'\nDownloaded: {success}/{len(all_downloads)}')
    if failed:
        print(f'Failed ({len(failed)}):')
        for kind, fname, url in failed:
            print(f'  {kind}: {fname} -> {url}')

    return failed


def main():
    if len(sys.argv) < 2:
        print('Usage: python convert_mediawiki.py <export.xml> [output_dir]')
        sys.exit(1)

    xml_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'docs'

    print(f'Parsing {xml_path}...')
    pages, redirects = parse_xml(xml_path)
    print(f'Found {len(pages)} pages and {len(redirects)} redirects')

    # Print redirect map
    print('\nRedirects:')
    for src, dst in sorted(redirects.items()):
        print(f'  {src} -> {dst}')

    all_files = set()
    all_media = set()

    # Convert all pages
    print(f'\nConverting pages to {output_dir}...')
    os.makedirs(output_dir, exist_ok=True)

    # Create supporting files for MkDocs
    js_dir = os.path.join(output_dir, 'javascripts')
    css_dir = os.path.join(output_dir, 'stylesheets')
    os.makedirs(js_dir, exist_ok=True)
    os.makedirs(css_dir, exist_ok=True)

    with open(os.path.join(js_dir, 'mathjax.js'), 'w', encoding='utf-8') as f:
        f.write("""window.MathJax = {
  tex: {
    inlineMath: [["\\\\(", "\\\\)"]],
    displayMath: [["\\\\[", "\\\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise()
})
""")

    with open(os.path.join(css_dir, 'extra.css'), 'w', encoding='utf-8') as f:
        f.write(""".figure-caption {
  text-align: center;
  font-style: italic;
  font-size: 0.9em;
  color: #666;
}

img {
  max-width: 100%;
  height: auto;
}

.md-typeset table {
  font-size: 0.8em;
}
""")

    for title, text in sorted(pages.items()):
        slug = slugify(title)
        if title == 'Home':
            filename = 'index.md'
        else:
            filename = f'{slug}.md'

        filepath = os.path.join(output_dir, filename)
        md = convert_wikitext(title, text, redirects, pages, all_files, all_media)

        # Add page title as H1 if not already starting with one
        if not md.startswith('# '):
            md = f'# {title}\n\n{md}'

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(md)
        print(f'  {title} -> {filename}')

    # Post-processing: fix known broken cross-references from source wiki
    print('\nApplying cross-reference fixups...')
    fixups = {
        'index.md': [
            ('range-statistics.md#extreme-spread', 'describing-precision.md#extreme-spread'),
            ('precision-models.md#how-large-a-sample-do-we-need', 'closed-form-precision.md#how-large-a-sample-do-we-need'),
        ],
        'closed-form-precision.md': [
            ('closed-form-precision.md#spread-measures', 'closed-form-precision.md#using'),
        ],
        'faq.md': [
            ('closed-form-precision.md#symmetric-bivariate-normal-rayleigh-distribution',
             'closed-form-precision.md#symmetric-bivariate-normal-shots-imply-rayleigh-distributed-distances'),
            ('precision-models.md#bessel-correction-factor', 'closed-form-precision.md#bessel-correction-factor'),
            ('precision-models.md#correction-factors', 'closed-form-precision.md#correction-factors'),
        ],
        'measuring-tools.md': [
            ('describing-precision.md#mean-radius', 'describing-precision.md#mean-radius-mr'),
        ],
    }
    for fname, replacements in fixups.items():
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            with open(fpath, 'r', encoding='utf-8') as f:
                content = f.read()
            for old, new in replacements:
                content = content.replace(old, new)
            with open(fpath, 'w', encoding='utf-8') as f:
                f.write(content)

    # Add missing anchor in herb-references.md (#nuttall1975 -> nuttall1975a)
    herb_path = os.path.join(output_dir, 'herb-references.md')
    if os.path.exists(herb_path):
        with open(herb_path, 'r', encoding='utf-8') as f:
            content = f.read()
        content = content.replace(
            '<a id="nuttall1975a">',
            '<a id="nuttall1975"></a><a id="nuttall1975a">'
        )
        with open(herb_path, 'w', encoding='utf-8') as f:
            f.write(content)

    print(f'\nReferenced images: {len(all_files)}')
    for f in sorted(all_files):
        print(f'  {f}')
    print(f'\nReferenced media files: {len(all_media)}')
    for f in sorted(all_media):
        print(f'  {f}')

    # Download files
    failed = download_all_files(all_files, all_media, output_dir)

    # Summary
    print(f'\n=== CONVERSION COMPLETE ===')
    print(f'Pages converted: {len(pages)}')
    print(f'Output directory: {output_dir}')
    if failed:
        print(f'WARNING: {len(failed)} files failed to download')

    return 0


if __name__ == '__main__':
    sys.exit(main())
