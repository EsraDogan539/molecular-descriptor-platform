"""Second reference check: candidate DOIs (Crossref), the ChalMolDB Zenodo DOI (DataCite), and the open-access full
texts of Haciefendioglu & Yildirim (JCIM 2025, PMC12152970; ACS Omega 2022, PMC9631416) from Europe PMC."""
import json, sys, time, urllib.parse, urllib.request
from pathlib import Path
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
H = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:esra.6615@gmail.com)"}
def fetch(url, js=True):
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=40) as r:
                b = r.read(); return json.loads(b) if js else b.decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404: return None
            time.sleep(3)
        except Exception: time.sleep(3)
    return "ERROR"
res = {}
for doi in (Path(__file__).parent / "refdata/candidates2.txt").read_text().split():
    m = fetch("https://api.crossref.org/works/" + urllib.parse.quote(doi))
    if m in (None, "ERROR"): res[doi] = {"status": "not_found" if m is None else "error"}; continue
    m = m["message"]
    res[doi] = {"status": "ok", "title": (m.get("title") or [""])[0], "journal": (m.get("short-container-title") or m.get("container-title") or [""])[0],
                "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0], "volume": m.get("volume"), "page": m.get("page") or m.get("article-number"),
                "authors": [(a.get("family"), a.get("given")) for a in m.get("author", [])]}
    time.sleep(0.1)
z = fetch("https://api.datacite.org/dois/10.5281/zenodo.23558298")
res["10.5281/zenodo.23558298"] = {"status": "ok", "title": z["data"]["attributes"]["titles"][0]["title"], "creators": [c["name"] for c in z["data"]["attributes"]["creators"]],
                                  "year": z["data"]["attributes"]["publicationYear"]} if isinstance(z, dict) else {"status": str(z)}
json.dump(res, open(OUT / "refs_verified2.json", "w"), indent=1, ensure_ascii=False)
for pmc in ["PMC12152970", "PMC9631416"]:
    t = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML", js=False)
    (OUT / f"{pmc}.xml").write_text(t if isinstance(t, str) else "NONE")
print("done", sum(v.get("status") == "ok" for v in res.values()), "/", len(res))
