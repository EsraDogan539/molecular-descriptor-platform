"""Literature scan: 2024-2026 chalcogen conjugated polymers with optical gaps, not in ChalMolDB."""
import json, re, time, urllib.parse, urllib.request
from pathlib import Path
OUT = Path("docs/lit"); (OUT / "xml").mkdir(parents=True, exist_ok=True)
have = set(Path("scripts/pol_dois.txt").read_text().split())
UA = {"User-Agent": "ChalMolDB-litscan/1.0 (mailto:esra.6615@gmail.com)"}
def get(url, raw=False):
    for k in range(3):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60).read()
            return r if raw else json.loads(r)
        except Exception as e:
            err = e; time.sleep(3)
    print("FAIL", url[:120], err); return None
found = {}
# ---- Europe PMC (abstract search, open-access flag)
chal = '(ABSTRACT:"selenophene" OR ABSTRACT:"benzoselenadiazole" OR ABSTRACT:"selenadiazole" OR ABSTRACT:"selenium" OR ABSTRACT:"benzothiadiazole" OR ABSTRACT:"thiophene" OR ABSTRACT:"diketopyrrolopyrrole" OR ABSTRACT:"benzodithiophene")'
gapq = '(ABSTRACT:"band gap" OR ABSTRACT:"bandgap" OR ABSTRACT:"band-gap")'
polq = '(ABSTRACT:"polymer" OR ABSTRACT:"copolymer" OR ABSTRACT:"polymers")'
q = f'{chal} AND {gapq} AND {polq} AND (ABSTRACT:"conjugated" OR ABSTRACT:"donor-acceptor" OR ABSTRACT:"D-A") AND PUB_YEAR:[2024 TO 2026]'
cursor = "*"
for page in range(12):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 1000, "resultType": "core", "cursorMark": cursor})
    d = get(u)
    if not d: break
    for r in d.get("resultList", {}).get("result", []):
        doi = (r.get("doi") or "").lower()
        if not doi: continue
        found[doi] = {"doi": doi, "title": r.get("title"), "year": r.get("pubYear"), "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"),
                      "abstract": r.get("abstractText", ""), "oa": r.get("isOpenAccess") == "Y", "pmcid": r.get("pmcid"), "src": "epmc", "in_db": doi in have}
    nc = d.get("nextCursorMark")
    if not nc or nc == cursor or not d.get("resultList", {}).get("result"): break
    cursor = nc
print("europepmc", len(found), flush=True)
# ---- Crossref (title/abstract bibliographic query)
for qq in ["selenophene conjugated polymer optical band gap", "benzoselenadiazole copolymer optical bandgap", "selenium substituted donor acceptor polymer solar cell",
           "selenophene polymer field-effect transistor", "selenophene-based polymer donor organic solar cells", "thiophene donor acceptor conjugated polymer optical bandgap synthesis",
           "benzothiadiazole copolymer optical band gap", "diketopyrrolopyrrole selenophene polymer", "electrochromic polymer selenophene", "low bandgap conjugated polymer near-infrared selenium"]:
    for off in range(0, 400, 100):
        u = "https://api.crossref.org/works?" + urllib.parse.urlencode({"query.bibliographic": qq, "filter": "from-pub-date:2024-01-01,type:journal-article", "rows": 100, "offset": off,
                                                                       "select": "DOI,title,issued,container-title,abstract"})
        d = get(u)
        if not d: break
        items = d["message"]["items"]
        for it in items:
            doi = it["DOI"].lower()
            if doi in found: continue
            found[doi] = {"doi": doi, "title": (it.get("title") or [""])[0], "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
                          "journal": (it.get("container-title") or [""])[0], "abstract": re.sub("<[^>]+>", " ", it.get("abstract", "") or ""), "oa": None, "pmcid": None, "src": "crossref", "in_db": doi in have}
        if len(items) < 100: break
        time.sleep(1)
print("total", len(found), flush=True)
# relevance filter
kw = re.compile(r"polymer|copolymer", re.I); ch = re.compile(r"selen|thiophen|thiadiazol|chalcogen|sulfur|dithiophene|thieno", re.I)
gap = re.compile(r"band ?-?gap|E\s*g|optical gap", re.I)
for v in found.values():
    t = (v["title"] or "") + " " + (v["abstract"] or "")
    v["rel"] = bool(kw.search(t) and ch.search(t)); v["has_gap_in_abs"] = bool(gap.search(v["abstract"] or ""))
    v["se"] = bool(re.search(r"selen", t, re.I))
    m = re.findall(r"(\d\.\d{1,2})\s*eV", v["abstract"] or ""); v["ev_values"] = m[:8]
json.dump(sorted(found.values(), key=lambda v: (not v["rel"], not v["se"], str(v["year"]))), open(OUT / "candidates.json", "w"), indent=0)
# full texts for OA PMC records that are relevant
n = 0
for v in found.values():
    if v["rel"] and v["pmcid"] and v["oa"] and not v["in_db"]:
        x = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{v['pmcid']}/fullTextXML", raw=True)
        if x: (OUT / "xml" / f"{v['pmcid']}.xml").write_bytes(x); n += 1
        time.sleep(0.3)
print("fulltexts", n)
