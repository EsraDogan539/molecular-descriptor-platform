import json, re, time, urllib.parse, urllib.request
refs = json.load(open("refcheck/refs.json"))
UA = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:hkayi@ankara.edu.tr)"}
def get(url):
    for _ in range(3):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30))
        except Exception as e:
            err = str(e); time.sleep(3)
    return {"error": err}
out = []
for i, r in enumerate(refs, 1):
    if len(r) == 5:
        a, j, y, v, p = r
        q = f"{a} {j} {y} {v} {p}"
    else:
        q = re.sub(r"[*]", "", r[0]); y = None
    url = "https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(q)
    d = get(url)
    hits = []
    for it in d.get("message", {}).get("items", [])[:3]:
        hits.append({
            "doi": it.get("DOI"), "title": (it.get("title") or [""])[0][:160],
            "journal": (it.get("container-title") or [""])[0], "short": (it.get("short-container-title") or [""])[:1],
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "volume": it.get("volume"), "issue": it.get("issue"), "page": it.get("page"), "article": it.get("article-number"),
            "authors": [f"{x.get('given','')} {x.get('family','')}".strip() for x in it.get("author", [])][:10],
            "n_authors": len(it.get("author", [])), "score": it.get("score")})
    out.append({"n": i, "ref": r, "hits": hits, "error": d.get("error")})
    time.sleep(1)
json.dump(out, open("docs/refcheck.json", "w"), ensure_ascii=False, indent=1)
