import json, re, time, urllib.request, urllib.parse
from pathlib import Path
out = Path("docs/litscan"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (ChalMolDB literature check; mailto:esra.6615@gmail.com)"}
def fetch(url, binary=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
            data = r.read()
            return data if binary else data.decode("utf-8", "replace")
    except Exception as e:
        return f"__ERR__ {e}"
DOIS = {
 "nanoscale_liu2025": "10.1039/d4nr03702b",
 "telluro_helicenes_jcc2025": "10.1002/jcc.70175",
 "se_te_bicyclic_polymers_omega2026": "10.1021/acsomega.5c04257",
 "nablacolors3d_commchem2026": "10.1038/s42004-026-01944-5",
 "deepacceptor_npj2024": "10.1038/s41524-024-01367-7",
 "hutchison_opv_jmca2024": "10.1039/d4ta01942c",
 "zhang_nfa_advsci2024": "10.1002/advs.202308652",
 "koskin_polymers2026": "10.3390/polym18070879",
 "meng_jcp2026": "10.1063/5.0333521",
 "saeki_dd2025": "10.1039/d5dd00418g",
 "deng_msde2025": "10.1039/d4me00188e",
 "peng_npj2024_fmo": "10.1038/s41524-024-01403-6",
 "ie_ea_aenm2025": "10.1002/aenm.202403707",
 "sizeresgraph_omega2026": "10.1021/acsomega.6c08651",
 "liu_pufp_npj2025": "10.1038/s41524-025-01608-3",
 "ringvault_chemsci2025": "10.1039/d5sc04079e",
 "ofetnow_jacs2026": "10.1021/jacs.6c11486",
 "cemp_jcheminf2026": "10.1186/s13321-026-01295-2",
 "solarrrl_544pairs_2025": "10.1002/solr.202400846",
 "cnn_absorption_jcim2024": "10.1021/acs.jcim.3c01792",
 "mapping_scp_mh2025": "10.1039/d5mh00485c",
 "polytellurophene_review_2025": "10.1039/d5py00439j",
 "kyhoiesh_jmgm2025": "10.1016/j.jmgm.2025.109154",
 "zindo_dml_jctc2026": "10.1021/acs.jctc.5c01587",
}
URLRE = re.compile(r"https?://(?:www\.)?(?:github\.com|zenodo\.org|doi\.org/10\.5281|figshare\.com|doi\.org/10\.6084|gitlab\.com|huggingface\.co|osf\.io|data\.mendeley\.com|doi\.org/10\.17632)[^\s\"'<>)\]]*", re.I)
res = {}
for key, doi in DOIS.items():
    rec = {"doi": doi}
    q = urllib.parse.quote(f'DOI:"{doi}"')
    s = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core")
    try:
        hit = json.loads(s)["resultList"]["result"][0]
        rec["title"] = hit.get("title"); rec["pmcid"] = hit.get("pmcid"); rec["license"] = hit.get("license")
        rec["abstract"] = (hit.get("abstractText") or "")[:3000]
    except Exception:
        hit = None
    text = ""
    if rec.get("pmcid"):
        xml = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{rec['pmcid']}/fullTextXML")
        if not xml.startswith("__ERR__"):
            text = xml
            secs = re.findall(r"<sec[^>]*>\s*<title>([^<]*(?:[Dd]ata|[Cc]ode|[Ss]upporting|[Aa]vailab|[Aa]ssociated)[^<]*)</title>(.*?)</sec>", xml, re.S)
            rec["availability_sections"] = [(t, re.sub(r"<[^>]+>", " ", b)[:1500]) for t, b in secs]
            notes = re.findall(r"<notes[^>]*>(.*?)</notes>", xml, re.S)
            rec["notes"] = [re.sub(r"<[^>]+>", " ", n)[:1200] for n in notes if re.search(r"[Aa]vailab", n)]
            rec["table_captions"] = [re.sub(r"<[^>]+>", " ", t)[:300] for t in re.findall(r"<caption>(.*?)</caption>", xml, re.S)][:25]
            rec["mentions_Se"] = bool(re.search(r"selen|\bSe\b", xml)); rec["mentions_Te"] = bool(re.search(r"tellur|\bTe\b", xml))
            rec["data_sentences"] = list(dict.fromkeys(m.strip() for m in re.findall(r"[^.]{0,200}(?:dataset|data set|database)[^.]{0,200}\.", re.sub(r"<[^>]+>", " ", xml))))[:15]
    if not text:
        oa = fetch("https://api.openalex.org/works/https://doi.org/" + doi)
        try:
            w = json.loads(oa)
            rec["title"] = rec.get("title") or w.get("title")
            loc = w.get("best_oa_location") or {}
            rec["oa_landing"] = loc.get("landing_page_url"); rec["oa_pdf"] = loc.get("pdf_url"); rec["oa_license"] = loc.get("license")
            page = fetch(loc.get("landing_page_url") or ("https://doi.org/" + doi))
            text = page if not page.startswith("__ERR__") else ""
            rec["landing_fetch"] = "ok" if text else page[:200]
        except Exception as e:
            rec["openalex_error"] = str(e)
    rec["repo_urls"] = sorted(set(u.rstrip(".,;") for u in URLRE.findall(text)))[:30]
    res[key] = rec
    time.sleep(1)
# Saeki Zenodo record
z = fetch("https://zenodo.org/api/records/17430843")
try:
    zj = json.loads(z); res["saeki_zenodo_17430843"] = {"title": zj["metadata"]["title"], "license": zj["metadata"].get("license"),
        "files": [(f["key"], f["size"]) for f in zj.get("files", [])], "description": re.sub(r"<[^>]+>", " ", zj["metadata"].get("description", ""))[:1500]}
except Exception as e:
    res["saeki_zenodo_17430843"] = {"error": str(e), "raw": z[:300]}
(out / "availability.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print({k: (v.get("pmcid"), len(v.get("repo_urls", []))) for k, v in res.items() if isinstance(v, dict)})
