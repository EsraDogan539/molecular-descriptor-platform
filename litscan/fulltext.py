import re, urllib.request
from pathlib import Path
out = Path("docs/litscan/fulltext"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (ChalMolDB literature check; mailto:esra.6615@gmail.com)"}
for key, pmc in {"helicenes_jcc2025": "PMC12223471", "se_te_bicyclic_omega2026": "PMC13000587"}.items():
    with urllib.request.urlopen(urllib.request.Request(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML", headers=UA), timeout=120) as r:
        xml = r.read().decode("utf-8", "replace")
    (out / f"{key}.xml").write_text(xml)
    txt = re.sub(r"<table-wrap", "\n\n[TABLE] <table-wrap", xml)
    txt = re.sub(r"</(p|title|tr|caption)>", "\n", txt)
    txt = re.sub(r"</t[dh]>", " | ", txt)
    txt = re.sub(r"<[^>]+>", "", txt)
    (out / f"{key}.txt").write_text(re.sub(r"\n\s*\n+", "\n\n", txt))
    # figures: list graphic hrefs
    figs = re.findall(r'<graphic[^>]*xlink:href="([^"]+)"', xml)
    (out / f"{key}_figures.txt").write_text("\n".join(figs))
    for f in figs[:12]:
        for ext in ("", ".jpg"):
            try:
                url = f"https://europepmc.org/articles/{pmc}/bin/{f}{ext}"
                with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                    (out / f"{key}_{f}{ext or ''}").write_bytes(r.read()); break
            except Exception:
                pass
print("ok")
