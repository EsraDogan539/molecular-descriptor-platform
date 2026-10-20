import json, time, urllib.parse, urllib.request
UA = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:hkayi@ankara.edu.tr)"}
def get(url):
    try: return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30))
    except Exception as e: return {"error": str(e)}
def summ(it):
    return {"doi": it.get("DOI"), "title": (it.get("title") or [""])[0][:160], "journal": (it.get("container-title") or [""])[0],
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0], "volume": it.get("volume"), "issue": it.get("issue"),
            "page": it.get("page"), "article": it.get("article-number"),
            "authors": [f"{x.get('given','')} {x.get('family','')}".strip() for x in it.get("author", [])][:12]}
out = {}
for doi in ["10.1038/sdata.2016.18", "10.1007/s10822-011-9440-2", "10.1186/s13321-015-0069-3", "10.1039/c7sc02664a"]:
    d = get("https://api.crossref.org/works/" + doi); out[doi] = summ(d.get("message", {})) if "message" in d else d
    time.sleep(1)
for name, q, extra in [
    ("ref12", "Kayi Ozkilinc Kaya Journal of the Turkish Chemical Society Section B 2017", ""),
    ("ref14", "Kayi Sen Ozkilinc chalcogen quinoxaline Journal of Molecular Modeling 2024", ""),
    ("ref16", "Sushko Online chemical modeling environment OCHEM web platform data storage model development", ""),
    ("ref39", "Pedregosa Scikit-learn Machine Learning in Python", "")]:
    d = get("https://api.crossref.org/works?rows=5&query.bibliographic=" + urllib.parse.quote(q) + extra)
    out[name] = [summ(it) for it in d.get("message", {}).get("items", [])]
    time.sleep(1)
d = get("https://api.datacite.org/dois/10.5281/zenodo.591637")
a = d.get("data", {}).get("attributes", {})
out["rdkit"] = {"titles": a.get("titles"), "creators": [c.get("name") for c in (a.get("creators") or [])][:5], "publisher": a.get("publisher"), "year": a.get("publicationYear"), "err": d.get("error")}
json.dump(out, open("docs/refcheck2.json", "w"), ensure_ascii=False, indent=1)
