import json, time, urllib.request
UA = {"User-Agent": "ChalMolDB-refcheck/1.0 (mailto:hkayi@ankara.edu.tr)"}
DOIS = ["10.1103/PhysRevLett.39.1098", "10.1038/347539a0", "10.1126/science.270.5243.1789", "10.1021/cr050149z",
        "10.1021/cm102419z", "10.1002/adma.200501717", "10.1021/cr950257t", "10.3390/polym9020039",
        "10.1007/s00894-014-2269-6", "10.1016/j.comptc.2014.12.011", "10.1016/j.comptc.2016.11.014", None,
        "10.1007/s00894-019-4043-2", "10.1007/s00894-024-05985-2", "10.1038/sdata.2016.18", "10.1007/s10822-011-9440-2",
        "10.1093/nar/gkr320", "10.1021/ci049714+", "10.1093/nar/gkac956", "10.1038/sdata.2014.22",
        "10.1021/acscentsci.9b00476", "10.1039/D3PY01372C", "10.1039/D4DD00266K", "10.1039/D4DD00196F",
        "10.1021/acs.jpcc.8b02913", "10.1021/acs.jcim.0c00726", "10.1038/s41524-017-0056-5", "10.1016/j.patter.2021.100238",
        "10.1038/s41467-023-39868-6", "10.1038/s41524-025-01823-y", "10.1021/acs.jcim.5c00345", "10.1186/s13321-015-0068-4",
        None, "10.1021/ci100050t", "10.1186/s13321-015-0069-3", "10.1023/A:1010933404324", "10.1007/BF00994018",
        "10.1145/2939672.2939785", None, "10.1021/acs.jcim.7b00403", "10.1021/ci400084k", "10.1039/C7SC02664A",
        "10.1093/bioinformatics/btu829", "10.1186/1758-2946-3-S1-P3"]
out = []
for doi in DOIS:
    if not doi: out.append(None); continue
    for _ in range(3):
        try:
            it = json.load(urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + doi, headers=UA), timeout=30))["message"]; break
        except Exception as e:
            it = {"error": str(e)}; time.sleep(3)
    out.append({"doi": doi, "title": (it.get("title") or [""])[0], "authors": [[x.get("given", ""), x.get("family", ""), x.get("suffix", "")] for x in it.get("author", [])],
                "journal": (it.get("container-title") or [""])[0], "short": (it.get("short-container-title") or [""])[0],
                "issued": (it.get("issued", {}).get("date-parts") or [[None]])[0], "published_print": (it.get("published-print", {}).get("date-parts") or [[None]])[0],
                "volume": it.get("volume"), "issue": it.get("issue"), "page": it.get("page"), "article": it.get("article-number"), "error": it.get("error")})
    time.sleep(0.5)
json.dump(out, open("docs/refmeta.json", "w"), ensure_ascii=False, indent=1)
