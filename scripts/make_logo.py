"""Generate the ChalMolDB mark and the home-page hero molecule.

Both show the same chemically valid structure: a thiophene–selenophene–tellurophene trimer
(2,2':5',2''-ter-chalcogenophene) in the anti conformation, the kind of chalcogenophene chain found
in the database's donor units. Each chalcogen has exactly two ring bonds; every ring carbon has three
bonds (ring + inter-ring bond or implicit H). Outputs:

  assets/ and static/chalmoldb_icon.svg, chalmoldb_logo.svg
  assets/hero_molecule.svg (the <g> block embedded in app_v08.py is regenerated with --hero)
"""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHALCOGENS = ("S", "Se", "Te")


def terchalcogenophene(rr, link):
    """Atoms {name: (element, x, y)} and bonds of the anti S/Se/Te trimer; screen y points down."""
    atoms, bonds = {}, []
    cx, cy = 0.0, 0.0
    for i, el in enumerate(CHALCOGENS):
        up = i % 2 == 0
        # X, Ca (right alpha), Cb, Cc, Cd (left alpha), in ring order
        angles = [90, 18, -54, -126, 162] if up else [-90, -18, 54, 126, 198]
        names = [f"X{i}", f"Ca{i}", f"Cb{i}", f"Cc{i}", f"Cd{i}"]
        for name, ang in zip(names, angles):
            x = cx + rr * math.cos(math.radians(ang))
            y = cy - rr * math.sin(math.radians(ang))
            atoms[name] = (el if name.startswith("X") else "C", x, y)
        bonds += [(names[k], names[(k + 1) % 5]) for k in range(5)]
        if i > 0:
            bonds.append((f"Ca{i - 1}", f"Cd{i}"))
        # next ring centre: along the line through this ring's right alpha carbon
        right = atoms[f"Ca{i}"]
        ux, uy = right[1] - cx, right[2] - cy
        norm = math.hypot(ux, uy)
        step = 2 * rr + link  # centre-to-centre along the alpha-alpha line
        cx, cy = cx + ux / norm * step, cy + uy / norm * step
    # centre the drawing on (0, 0)
    xs = [a[1] for a in atoms.values()]
    ys = [a[2] for a in atoms.values()]
    ox, oy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    atoms = {k: (e, x - ox, y - oy) for k, (e, x, y) in atoms.items()}
    return atoms, bonds


def defs(p, shadow=(2.5, 2.5, 0.18)):
    dy, sd, op = shadow
    return f"""<defs>
<radialGradient id="{p}C" cx="32%" cy="25%" r="72%"><stop stop-color="#ffffff"/><stop offset=".48" stop-color="#dbe1e7"/><stop offset="1" stop-color="#8b96a1"/></radialGradient>
<radialGradient id="{p}S" cx="30%" cy="24%" r="72%"><stop stop-color="#fff5b6"/><stop offset=".42" stop-color="#e7bd39"/><stop offset="1" stop-color="#ad7e00"/></radialGradient>
<radialGradient id="{p}Se" cx="30%" cy="24%" r="72%"><stop stop-color="#c6fff6"/><stop offset=".42" stop-color="#39aa9e"/><stop offset="1" stop-color="#0d6b63"/></radialGradient>
<radialGradient id="{p}Te" cx="30%" cy="24%" r="72%"><stop stop-color="#eee6f7"/><stop offset=".42" stop-color="#846b9f"/><stop offset="1" stop-color="#55406d"/></radialGradient>
<linearGradient id="{p}Bond" x1="0" x2="1"><stop stop-color="#99a4af"/><stop offset=".5" stop-color="#cbd2d9"/><stop offset="1" stop-color="#7f8a95"/></linearGradient>
<filter id="{p}Shadow" x="-30%" y="-30%" width="160%" height="160%"><feDropShadow dx="0" dy="{dy}" stdDeviation="{sd}" flood-color="#173b5b" flood-opacity="{op}"/></filter>
</defs>"""


def is_double(a, b):
    """Kekule structure of each ring: X-Ca, Ca=Cb, Cb-Cc, Cc=Cd, Cd-X; inter-ring bonds single."""
    if a[-1] != b[-1]:
        return False
    return {a[:2], b[:2]} in ({"Ca", "Cb"}, {"Cc", "Cd"})


def molecule(p, rr, link, bond_w, c_r, x_r, font, ox, oy, edge=1.6, double_bonds=False):
    atoms, bonds = terchalcogenophene(rr, link)
    n = lambda v: f"{v:.1f}"  # noqa: E731
    out = [f'<g transform="translate({n(ox)} {n(oy)})" filter="url(#{p}Shadow)" stroke-linecap="round">',
           f'<g stroke="url(#{p}Bond)" stroke-width="{bond_w}">']
    for a, b in bonds:
        _, x1, y1 = atoms[a]
        _, x2, y2 = atoms[b]
        if double_bonds and is_double(a, b):
            # two thinner parallel sticks
            length = math.hypot(x2 - x1, y2 - y1)
            nx, ny = -(y2 - y1) / length, (x2 - x1) / length
            off, w = 0.6 * bond_w, 0.45 * bond_w
            for sgn in (1, -1):
                out.append(f'<path stroke-width="{w:.1f}" d="M{n(x1 + sgn * off * nx)} {n(y1 + sgn * off * ny)} '
                           f'{n(x2 + sgn * off * nx)} {n(y2 + sgn * off * ny)}"/>')
            continue
        out.append(f'<path d="M{n(x1)} {n(y1)} {n(x2)} {n(y2)}"/>')
    out.append(f'</g><g stroke="#ffffff" stroke-width="{edge}">')
    for el, x, y in atoms.values():
        out.append(f'<circle cx="{n(x)}" cy="{n(y)}" r="{x_r if el != "C" else c_r}" fill="url(#{p}{el})"/>')
    out.append('</g><g font-family="Arial,Helvetica,sans-serif" font-weight="800" text-anchor="middle" '
               'dominant-baseline="central" fill="#ffffff">')
    for el, x, y in atoms.values():
        if el != "C":
            out.append(f'<text x="{n(x)}" y="{n(y + 0.04 * font)}" font-size="{font if el == "S" else font * 0.88:.1f}">{el}</text>')
    out.append("</g></g>")
    return "\n".join(out), atoms, bonds


def check_valence(atoms, bonds):
    degree = {k: 0 for k in atoms}
    for a, b in bonds:
        degree[a] += 1
        degree[b] += 1
    for k, (el, _, _) in atoms.items():
        if el != "C":
            assert degree[k] == 2, (k, degree[k])
        else:
            assert degree[k] in (2, 3), (k, degree[k])  # 2 = CH
    assert sum(1 for e, _, _ in atoms.values() if e != "C") == 3 and len(bonds) == 17


def main():
    # icon: 160 x 120
    icon_body, atoms, bonds = molecule("cm", rr=17, link=17, bond_w=4.5, c_r=5, x_r=11.5, font=11.5, ox=80, oy=60)
    check_valence(atoms, bonds)
    head = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 160 120" role="img" aria-label="ChalMolDB molecular mark">'
    (ROOT / "assets/chalmoldb_icon.svg").write_text(f"{head}\n{defs('cm')}\n{icon_body}\n</svg>\n")
    tile = '<rect width="160" height="120" rx="20" fill="#F8FAFC"/>'
    (ROOT / "static/chalmoldb_icon.svg").write_text(f"{head}\n{tile}\n{defs('cm')}\n{icon_body}\n</svg>\n")

    logo_body, _, _ = molecule("lg", rr=17, link=17, bond_w=4.5, c_r=5, x_r=11.5, font=11.5, ox=95, oy=75)
    logo = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 150" role="img" aria-labelledby="t d">
<title id="t">ChalMolDB</title><desc id="d">Chalcogen Molecular Database</desc>
{defs('lg')}
{logo_body}
<text x="205" y="69" font-family="Inter,Arial,Helvetica,sans-serif" font-size="52" font-weight="800" letter-spacing="-1.6" fill="#153A5B">ChalMolDB</text>
<text x="208" y="101" font-family="Inter,Arial,Helvetica,sans-serif" font-size="18" font-weight="500" letter-spacing=".3" fill="#6B7582">Chalcogen Molecular Database</text>
</svg>
"""
    for folder in ("assets", "static"):
        (ROOT / folder / "chalmoldb_logo.svg").write_text(logo)

    hero_body, _, _ = molecule("hero", rr=56, link=58, bond_w=13, c_r=19, x_r=36, font=23, ox=380, oy=280,
                               edge=3, double_bonds=True)
    hero = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 760 540">\n'
            f'{defs("hero", (9, 10, 0.20))}\n{hero_body}\n</svg>\n')
    (ROOT / "assets/hero_molecule.svg").write_text(hero)
    if "--hero" in sys.argv:
        print(hero_body)


if __name__ == "__main__":
    main()
