import json, time, urllib.request
UA = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:hkayi@ankara.edu.tr)"}
def get(url):
    err = None
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30))
        except Exception as e: err = str(e); time.sleep(3)
    return {"error": err}
def summ(it):
    return {"doi": it.get("DOI"), "title": (it.get("title") or [""])[0], "journal": (it.get("container-title") or [""])[0],
            "short": (it.get("short-container-title") or [""])[0],
            "issued": (it.get("issued", {}).get("date-parts") or [[None]])[0], "print": (it.get("published-print", {}).get("date-parts") or [[None]])[0],
            "volume": it.get("volume"), "issue": it.get("issue"), "page": it.get("page"), "article": it.get("article-number"),
            "authors": [[x.get("given", ""), x.get("family", ""), x.get("suffix", "")] for x in it.get("author", [])]}
D = {"LIU": "10.1039/d4nr03702b", "SAEKI": "10.1021/acs.jpclett.8b00635", "KUMAR": "10.1002/jcc.70175",
     "BUYUK": "10.1021/acsomega.5c04257", "CEMP": "10.1186/s13321-026-01295-2", "OFET": "10.1021/jacs.6c11486",
     "STJOHN": "10.1063/1.5099132"}
out = {}
for k, doi in D.items():
    d = get("https://api.crossref.org/works/" + doi)
    out[k] = [summ(d["message"])] if "message" in d else d
    time.sleep(1)
json.dump(out, open("docs/newrefs2.json", "w"), ensure_ascii=False, indent=1)
print({k: (v[0]["title"][:60] if isinstance(v, list) else v) for k, v in out.items()})
