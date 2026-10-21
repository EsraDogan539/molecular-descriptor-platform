import io, json, re, zipfile, urllib.request
from pathlib import Path
out = Path("docs/litscan"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (ChalMolDB literature check; mailto:esra.6615@gmail.com)"}
def fetch(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
            return r.read()
    except Exception as e:
        return ("__ERR__ " + str(e)).encode()
res = {"github": {}, "supp": {}}
for repo in ["Levitsiy/PolymersPropertiesPrediction", "AI4DD/nablaColors", "HMakkiMD/SCP_HT", "pkzheng94/ring", "Songyosk/UVVIS"]:
    meta = fetch(f"https://api.github.com/repos/{repo}")
    tree = fetch(f"https://api.github.com/repos/{repo}/git/trees/HEAD?recursive=1")
    try:
        m = json.loads(meta); t = json.loads(tree)
        res["github"][repo] = {"exists": "id" in m, "license": (m.get("license") or {}).get("spdx_id"), "pushed": m.get("pushed_at"),
                               "files": [(x["path"], x.get("size")) for x in t.get("tree", []) if x["type"] == "blob"][:120]}
    except Exception as e:
        res["github"][repo] = {"error": str(e), "raw": meta[:200].decode("utf-8", "replace")}
for key, pmc in {"telluro_helicenes": "PMC12223471", "se_te_bicyclic": "PMC13000587", "zhang_nfa": "PMC11077656"}.items():
    data = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles")
    info = {"bytes": len(data)}
    try:
        z = zipfile.ZipFile(io.BytesIO(data))
        info["files"] = [(i.filename, i.file_size) for i in z.infolist()]
        sub = Path(out / "supp" / key); sub.mkdir(parents=True, exist_ok=True)
        for i in z.infolist():
            if i.file_size < 15_000_000:
                (sub / Path(i.filename).name).write_bytes(z.read(i))
    except Exception as e:
        info["error"] = str(e); info["head"] = data[:200].decode("utf-8", "replace")
    res["supp"][key] = info
# NREL polymer dataset (St. John et al. 2019) record
for url in ["https://api.openalex.org/works/https://doi.org/10.1063/1.5099132"]:
    try:
        w = json.loads(fetch(url)); res["stjohn2019"] = {"title": w.get("title"), "year": w.get("publication_year"), "cited_by": w.get("cited_by_count")}
    except Exception as e:
        res["stjohn2019"] = {"error": str(e)}
(out / "availability2.json").write_text(json.dumps(res, indent=1))
print(json.dumps({k: (v if k != "github" else {r: (x.get("exists"), x.get("license"), len(x.get("files", []))) for r, x in v.items()}) for k, v in res.items()}, indent=1)[:3000])
