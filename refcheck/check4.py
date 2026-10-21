import json, time, urllib.parse, urllib.request
UA = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:hkayi@ankara.edu.tr)"}
def get(url):
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
Q = {
 "CEP": "Hachmann Harvard Clean Energy Project large-scale computational screening design organic photovoltaics World Community Grid J. Phys. Chem. Lett. 2011",
 "HOPV15": "Lopez Pyzer-Knapp Harvard organic photovoltaic dataset Scientific Data 2016",
 "OE62": "Stuke Atomic structures and orbital energies of 61,489 crystal-forming organic molecules Scientific Data 2020",
 "DTU": "Jorgensen Machine learning-based screening of complex molecules for polymer solar cells J. Chem. Phys. 2018",
 "OCELOT": "OCELOT infrastructure for data-driven research to discover and design crystalline organic semiconductors J. Chem. Phys. 2021",
}
out = {}
for k, q in Q.items():
    d = get("https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(q))
    out[k] = [summ(it) for it in d.get("message", {}).get("items", [])] or d
    time.sleep(1)
json.dump(out, open("docs/newrefs.json", "w"), ensure_ascii=False, indent=1)
