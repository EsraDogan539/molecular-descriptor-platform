import gzip, io, json, re, urllib.request
from pathlib import Path
out = Path("docs/litscan"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (ChalMolDB literature check; mailto:esra.6615@gmail.com)"}
def fetch(url, limit=None):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
            return r.read(limit) if limit else r.read()
    except Exception as e:
        return ("__ERR__ " + str(e)).encode()
res = {}
html = fetch("https://data.nrel.gov/submissions/236").decode("utf-8", "replace")
links = sorted(set(re.findall(r'href="([^"]+)"', html)))
res["nrel_links"] = [l for l in links if any(k in l for k in ("csv", "files", "license", "License", "download"))][:40]
res["nrel_license_text"] = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))[:3000]
for l in res["nrel_links"]:
    if "smiles_test" in l or "smiles_valid" in l:
        url = l if l.startswith("http") else "https://data.nrel.gov" + l
        data = fetch(url)
        try:
            txt = gzip.decompress(data).decode()
            lines = txt.splitlines()
            res["nrel_smiles_sample"] = {"url": url, "rows": len(lines) - 1, "header": lines[0][:200],
                "with_Se": sum(("[se]" in x or "[Se]" in x) for x in lines), "with_Te": sum(("[te]" in x or "[Te]" in x) for x in lines),
                "with_S": sum(bool(re.search(r"(?<![A-Za-z\[])[sS](?![eincbr])", x)) for x in lines), "example": lines[1][:200]}
        except Exception as e:
            res["nrel_smiles_sample"] = {"url": url, "error": str(e), "head": data[:200].decode("latin1")}
        break
z = fetch("https://zenodo.org/api/records/18061300")
try:
    zj = json.loads(z); res["nablacolors_zenodo"] = {"title": zj["metadata"]["title"], "license": zj["metadata"].get("license"),
        "files": [(f["key"], f["size"]) for f in zj.get("files", [])]}
except Exception as e:
    res["nablacolors_zenodo"] = {"error": str(e), "raw": z[:200].decode("latin1")}
for key, pmc in {"se_te_bicyclic": "PMC13000587", "telluro_helicenes": "PMC12223471"}.items():
    xml = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML").decode("utf-8", "replace")
    tabs = re.findall(r"<table-wrap[^>]*>(.*?)</table-wrap>", xml, re.S)
    res[key + "_tables"] = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " | ", t))[:4000] for t in tabs]
(out / "availability3.json").write_text(json.dumps(res, indent=1))
print({k: (len(v) if isinstance(v, list) else (v if not isinstance(v, str) else v[:200])) for k, v in res.items()})
