"""Shared presentation settings for figures and self-contained HTML reports."""
from pathlib import Path

def figure_font_path():
    from matplotlib import font_manager
    for family in ['Arial','Liberation Sans','DejaVu Sans']:
        try:return font_manager.findfont(family,fallback_to_default=False)
        except ValueError:pass
    return font_manager.findfont('sans-serif')

def configure_matplotlib():
    import matplotlib
    from matplotlib import font_manager
    family=font_manager.FontProperties(fname=figure_font_path()).get_name()
    matplotlib.rcParams.update({'font.family':family,'font.size':11,'axes.titleweight':'normal','figure.titleweight':'normal'})

def report_css():
    path=Path(__file__).resolve().parents[1]/'docs/assets/docs.css'
    css=path.read_text() if path.exists() else 'body{font:17px/1.65 Arial,sans-serif;background:white;color:#26333c}main{max-width:1140px;margin:auto}'
    return css+'\nheader{max-width:1140px;margin:auto;background:white;border:0}section{border:0;border-top:1px solid #e1e5e8;border-radius:0;padding:20px 0;margin:24px 0}img{border:0}a.button{border-radius:2px}'
