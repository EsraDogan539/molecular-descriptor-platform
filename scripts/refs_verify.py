"""Verify references against Crossref: (A) source DOIs of the 227 recent polymers in ChalMolDB,
(B) DOIs cited in the thesis manuscript, (C) candidate DOIs, (D) bibliographic searches (top 3 hits)."""
import glob, json, sys, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd
ROOT, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
H = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:esra.6615@gmail.com)"}
def get(url):
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=30) as r: return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404: return None
            time.sleep(2 + 3 * k)
        except Exception: time.sleep(2 + 3 * k)
    return "ERROR"
def meta(m):
    if m in (None, "ERROR"): return {"status": "not_found" if m is None else "error"}
    m = m["message"]
    return {"status": "ok", "title": (m.get("title") or [""])[0], "journal": (m.get("container-title") or [""])[0],
            "short_journal": (m.get("short-container-title") or [""])[0], "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "volume": m.get("volume"), "issue": m.get("issue"), "page": m.get("page"), "article_number": m.get("article-number"),
            "authors": [(a.get("family"), a.get("given")) for a in m.get("author", [])], "type": m.get("type"), "publisher": m.get("publisher")}
def doi_meta(d): return meta(get("https://api.crossref.org/works/" + urllib.parse.quote(d)))
res = {"A": {}, "B": {}, "C": {}, "D": {}}
src = ROOT / "data/sources/liu_2025_nanoscale/new_227"
d = pd.concat([pd.read_csv(f).rename(columns=str.lower).assign(file=Path(f).name) for f in glob.glob(str(src / "new-polymer-*.csv"))])
for ref, g in d.groupby("reference"):
    key = ref.split("doi.org/")[-1] if "doi.org/" in ref else ref
    res["A"][key] = {"n_polymers": int(len(g)), "names": list(g["name"].astype(str)), "groups": sorted(set(g.file)), **(doi_meta(key) if "doi.org/" in ref else {"status": "no_doi"})}
    time.sleep(0.1)
print("A", sum(v["status"] == "ok" for v in res["A"].values()), "/", len(res["A"]), flush=True)
for part, fn in [("B", "ms_dois.txt"), ("C", "candidates.txt")]:
    for doi in (Path(__file__).parent / "refdata" / fn).read_text().split():
        res[part][doi] = doi_meta(doi); time.sleep(0.1)
    print(part, sum(v["status"] == "ok" for v in res[part].values()), "/", len(res[part]), flush=True)
for q in (Path(__file__).parent / "refdata/searches.txt").read_text().strip().splitlines():
    r = get("https://api.crossref.org/works?rows=5&select=DOI,title,container-title,issued,author,is-referenced-by-count&query.bibliographic=" + urllib.parse.quote(q))
    items = [] if r in (None, "ERROR") else r["message"]["items"]
    res["D"][q] = [{"doi": i["DOI"], "title": (i.get("title") or [""])[0], "journal": (i.get("container-title") or [""])[0],
                    "year": (i.get("issued", {}).get("date-parts") or [[None]])[0][0], "cites": i.get("is-referenced-by-count"),
                    "first_author": (i.get("author") or [{}])[0].get("family")} for i in items]
json.dump(res, open(OUT / "refs_verified.json", "w"), indent=1, ensure_ascii=False)
print("done")
