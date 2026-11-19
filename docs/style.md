# Documentation style

HTML pages share `assets/docs.css`: Arial, 17 px body text, a white background, and regular-weight headings. GitHub supplies its own font for Markdown pages.

Use the same font in figure and animation labels. Keep titles short, give each figure a clear purpose, and define its colors and measurements nearby. Put acquisition details and measurement definitions in collapsed sections.

Scientific images, display scales, and measured values are preserved when changing presentation. Coregistration animations keep their original frames and timing.

Build the HTML pages from the repository root:

```bash
python tools/build_documentation.py
```
