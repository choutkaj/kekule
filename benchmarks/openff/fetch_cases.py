"""Fetch supplementary external PubChem cases; preserve response bytes and provenance."""
import hashlib,json,urllib.request,urllib.parse
from pathlib import Path

names=["cyclopropane","cyclobutane","cubane","norbornane","naphthalene","cis-2-butene","trans-2-butene","D-glucose","L-cysteine","dimethyl sulfoxide"]
here=Path(__file__).resolve().parent
destination=here/"fixtures/pubchem"
destination.mkdir(exist_ok=True)
sources=[]
for i,name in enumerate(names):
    url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+urllib.parse.quote(name,safe="")+"/property/IsomericSMILES/JSON"
    request=urllib.request.Request(url,headers={"User-Agent":"kekule-openff-scientific-reference/0.1"})
    with urllib.request.urlopen(request,timeout=30) as response: data=response.read()
    path=destination/f"{i}.json";path.write_bytes(data)
    sources.append(dict(name=name,path=str(path.relative_to(here)).replace("\\","/"),url=url,sha256=hashlib.sha256(data).hexdigest()))
(here/"supplementary-sources.lock.json").write_text(json.dumps(dict(schema=1,sources=sources),indent=2)+"\n",encoding="utf-8",newline="\n")
