import json, time, urllib.request, urllib.parse
from pathlib import Path
out = Path("docs/litscan"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "ChalMolDB-litscan (mailto:esra.6615@gmail.com)"}
def get(url):
    err = None
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return json.load(r)
        except Exception as e:
            err = e; time.sleep(3)
    return {"error": str(err)}
def inv(w):
    idx = w.get("abstract_inverted_index") or {}
    pos = sorted(((p, k) for k, ps in idx.items() for p in ps))
    return " ".join(k for _, k in pos)
def slim(w):
    return {"doi": w.get("doi"), "title": w.get("title"), "year": w.get("publication_year"),
            "venue": ((w.get("primary_location") or {}).get("source") or {}).get("display_name"),
            "cited_by": w.get("cited_by_count"), "type": w.get("type"),
            "oa": (w.get("open_access") or {}).get("oa_status"),
            "authors": [a["author"]["display_name"] for a in w.get("authorships", [])][:6],
            "abstract": inv(w)[:2500]}
def citing(doi, since="2024-01-01"):
    w = get("https://api.openalex.org/works/https://doi.org/" + doi)
    if "id" not in w: return {"error": w}
    wid = w["id"].split("/")[-1]; res = []; cur = "*"
    while cur and len(res) < 800:
        r = get(f"https://api.openalex.org/works?filter=cites:{wid},from_publication_date:{since}&per-page=200&cursor={cur}")
        res += [slim(x) for x in r.get("results", [])]; cur = (r.get("meta") or {}).get("next_cursor")
    return {"work": slim(w), "citing": res}
res = {"citing": {}}
seeds = {"nanoscale": "10.1039/d4nr03702b", "erol_jcim": "10.1021/acs.jcim.5c00345",
         "saeki_jpcl2018": "10.1021/acs.jpclett.8b00635", "hopv15": "10.1038/sdata.2016.86",
         "deng_msde2025": "10.1039/d4me00188e", "cep_lopez2016": "10.1016/j.joule.2017.10.006"}
for k, d in seeds.items():
    res["citing"][k] = citing(d)
queries = [
 "dataset conjugated polymers optical gap experimental", "database organic semiconductors HOMO LUMO open data",
 "machine learning band gap conjugated polymers", "machine learning optical gap organic molecules dataset",
 "machine learning HOMO LUMO prediction organic photovoltaic", "graph neural network conjugated polymer electronic properties",
 "non-fullerene acceptor database energy levels", "polymer donor dataset photovoltaic machine learning",
 "selenophene machine learning", "selenium-containing organic semiconductors dataset", "tellurophene", 
 "chalcogen organic semiconductors machine learning", "benchmark dataset excited states organic chromophores",
 "TDDFT dataset optical gap organic molecules", "experimental absorption dataset dyes machine learning",
 "open dataset organic electronics Zenodo", "high-throughput virtual screening donor acceptor polymers DFT",
 "conjugated polymer bandgap prediction transfer learning", "polymer electron affinity ionization energy dataset",
 "OPV materials database HOMO LUMO bandgap machine learning 2025"]
res["search"] = {}
for q in queries:
    r = get("https://api.openalex.org/works?search=" + urllib.parse.quote(q) +
            "&filter=from_publication_date:2024-01-01,type:article&sort=relevance_score:desc&per-page=50")
    res["search"][q] = [slim(x) for x in r.get("results", [])]
(out / "openalex2.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print({k: len(v.get("citing", [])) for k, v in res["citing"].items()})
