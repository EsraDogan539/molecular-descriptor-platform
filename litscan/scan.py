import json, time, urllib.request, urllib.parse
from pathlib import Path
out = Path("docs/litscan"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "ChalMolDB-litscan (mailto:esra.6615@gmail.com)"}
def get(url):
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return json.load(r)
        except Exception as e:
            err = e; time.sleep(3)
    return {"error": str(err)}
def oa_work(doi):
    return get("https://api.openalex.org/works/https://doi.org/" + doi)
def slim(w):
    return {"id": w.get("id"), "doi": w.get("doi"), "title": w.get("title"), "year": w.get("publication_year"),
            "venue": ((w.get("primary_location") or {}).get("source") or {}).get("display_name"),
            "cited_by": w.get("cited_by_count"), "type": w.get("type"),
            "authors": [a["author"]["display_name"] for a in w.get("authorships", [])][:8],
            "abstract": " ".join(k for k, _ in sorted(((k, p) for k, ps in (w.get("abstract_inverted_index") or {}).items() for p in ps), key=lambda x: x[1]))[:1500]}
res = {}
seeds = {"nanoscale": "10.1039/d4nr03702b", "erol_jcim": "10.1021/acs.jcim.5c00345", "saeki_jpcl": "10.1021/acs.jpclett.8b00635"}
for k, doi in seeds.items():
    w = oa_work(doi)
    if "id" not in w:
        res[k] = {"error": w}; continue
    wid = w["id"].split("/")[-1]
    cites = []
    cur = "*"
    while cur:
        r = get(f"https://api.openalex.org/works?filter=cites:{wid}&per-page=200&cursor={cur}")
        cites += [slim(x) for x in r.get("results", [])]
        cur = (r.get("meta") or {}).get("next_cursor")
        if len(cites) > 600: break
    res[k] = {"work": slim(w), "citing": cites}
queries = [
    "machine learning optical gap conjugated polymers dataset",
    "machine learning band gap organic semiconductors selenophene",
    "tellurophene DFT HOMO LUMO dataset",
    "selenophene donor acceptor DFT machine learning",
    "organic photovoltaic donor acceptor dataset DFT HOMO LUMO machine learning 2025",
    "chalcogen substitution band gap DFT oligomer machine learning",
    "non-fullerene acceptor dataset DFT machine learning",
    "conjugated polymer database HOMO LUMO experimental",
    "benzothiadiazole benzoselenadiazole DFT dataset",
    "optoelectronic property prediction small molecules organic solar cells dataset",
]
res["search"] = {}
for q in queries:
    r = get("https://api.openalex.org/works?search=" + urllib.parse.quote(q) + "&filter=from_publication_date:2023-01-01&sort=relevance_score:desc&per-page=40")
    res["search"][q] = [slim(x) for x in r.get("results", [])]
(out / "openalex.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print({k: (len(v["citing"]) if isinstance(v, dict) and "citing" in v else v) for k, v in res.items() if k != "search"})
