"""Generate the ChalMolDB mark (assets/ and static/ icon + logo SVGs).

The mark is a complete fused bicyclic ring (two hexagons sharing one edge, ten atoms, eleven bonds)
with S, Se and Te atoms, drawn in the same shaded-sphere style as the hero molecule of the home page.
"""

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = 34
H = R * math.sqrt(3) / 2
CX_L, CY = 48, 70
CX_R = CX_L + 2 * H


def hexagon(cx):
    return [(cx + R * math.cos(math.radians(a)), CY - R * math.sin(math.radians(a)))
            for a in (90, 30, -30, -90, -150, 150)]


L, RR = hexagon(CX_L), hexagon(CX_R)
ATOMS = {"a0": L[0], "a1": L[1], "a2": L[2], "a3": L[3], "a4": L[4], "a5": L[5],
         "b0": RR[0], "b1": RR[1], "b2": RR[2], "b3": RR[3]}
BONDS = [("a0", "a1"), ("a1", "a2"), ("a2", "a3"), ("a3", "a4"), ("a4", "a5"), ("a5", "a0"),
         ("a1", "b0"), ("b0", "b1"), ("b1", "b2"), ("b2", "b3"), ("b3", "a2")]
HETERO = {"a0": ("S", 13, 13), "b0": ("Se", 14, 11.5), "b3": ("Te", 14, 11.5)}
WIDTH = round(CX_R + H + 18)
HEIGHT = 140


def defs(p):
    return f"""<defs>
<radialGradient id="{p}C" cx="32%" cy="25%" r="72%"><stop stop-color="#ffffff"/><stop offset=".48" stop-color="#dbe1e7"/><stop offset="1" stop-color="#8b96a1"/></radialGradient>
<radialGradient id="{p}S" cx="30%" cy="24%" r="72%"><stop stop-color="#fff5b6"/><stop offset=".42" stop-color="#e7bd39"/><stop offset="1" stop-color="#ad7e00"/></radialGradient>
<radialGradient id="{p}Se" cx="30%" cy="24%" r="72%"><stop stop-color="#c6fff6"/><stop offset=".42" stop-color="#39aa9e"/><stop offset="1" stop-color="#0d6b63"/></radialGradient>
<radialGradient id="{p}Te" cx="30%" cy="24%" r="72%"><stop stop-color="#eee6f7"/><stop offset=".42" stop-color="#846b9f"/><stop offset="1" stop-color="#55406d"/></radialGradient>
<linearGradient id="{p}Bond" x1="0" x2="1"><stop stop-color="#99a4af"/><stop offset=".5" stop-color="#cbd2d9"/><stop offset="1" stop-color="#7f8a95"/></linearGradient>
<filter id="{p}Shadow" x="-30%" y="-30%" width="160%" height="160%"><feDropShadow dx="0" dy="2.5" stdDeviation="2.5" flood-color="#173b5b" flood-opacity=".18"/></filter>
</defs>"""


def mark(p, dx=0, dy=0):
    n = lambda v: f"{v:.1f}"  # noqa: E731
    out = [f'<g transform="translate({dx} {dy})" filter="url(#{p}Shadow)" stroke-linecap="round">',
           f'<g stroke="url(#{p}Bond)" stroke-width="6">']
    out += [f'<path d="M{n(ATOMS[a][0])} {n(ATOMS[a][1])} {n(ATOMS[b][0])} {n(ATOMS[b][1])}"/>' for a, b in BONDS]
    out.append('</g><g stroke="#ffffff" stroke-width="1.6">')
    for key, (x, y) in ATOMS.items():
        if key in HETERO:
            el, radius, _ = HETERO[key]
            out.append(f'<circle cx="{n(x)}" cy="{n(y)}" r="{radius}" fill="url(#{p}{el})"/>')
        else:
            out.append(f'<circle cx="{n(x)}" cy="{n(y)}" r="7.5" fill="url(#{p}C)"/>')
    out.append('</g><g font-family="Arial,Helvetica,sans-serif" font-weight="800" text-anchor="middle" '
               'dominant-baseline="central" fill="#ffffff">')
    for key, (el, _, size) in HETERO.items():
        x, y = ATOMS[key]
        out.append(f'<text x="{n(x)}" y="{n(y + 0.5)}" font-size="{size}">{el}</text>')
    out.append("</g></g>")
    return "\n".join(out)


def main():
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-label="ChalMolDB molecular mark">'
    icon = f"{head}\n{defs('cm')}\n{mark('cm')}\n</svg>\n"
    (ROOT / "assets/chalmoldb_icon.svg").write_text(icon)
    tile = f'<rect width="{WIDTH}" height="{HEIGHT}" rx="20" fill="#F8FAFC"/>'
    (ROOT / "static/chalmoldb_icon.svg").write_text(f"{head}\n{tile}\n{defs('cm')}\n{mark('cm')}\n</svg>\n")
    logo = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 150" role="img" aria-labelledby="t d">
<title id="t">ChalMolDB</title><desc id="d">Chalcogen Molecular Database</desc>
{defs('lg')}
{mark('lg', 8, 5)}
<text x="205" y="69" font-family="Inter,Arial,Helvetica,sans-serif" font-size="52" font-weight="800" letter-spacing="-1.6" fill="#153A5B">ChalMolDB</text>
<text x="208" y="101" font-family="Inter,Arial,Helvetica,sans-serif" font-size="18" font-weight="500" letter-spacing=".3" fill="#6B7582">Chalcogen Molecular Database</text>
</svg>
"""
    for folder in ("assets", "static"):
        (ROOT / folder / "chalmoldb_logo.svg").write_text(logo)
    print("viewBox", WIDTH, HEIGHT)


if __name__ == "__main__":
    main()
