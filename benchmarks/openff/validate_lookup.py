"""Exercise fixed-H InChI generation against every entry of the pinned Ash table."""
import argparse
import json
import subprocess
from pathlib import Path
from audit import MODEL_HASH, digest

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    model=json.loads((args.model/"model.json").read_text())
    assert model["checkpoint_sha256"]==MODEL_HASH
    entries=model["lookup_tables"]["am1bcc_charges"]
    process=subprocess.Popen([str(args.binary.resolve()),str(args.model.resolve())],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    results=[]
    for i,entry in enumerate(entries):
        process.stdin.write(json.dumps(dict(smiles=entry["mapped_smiles"],mode="identity"))+"\n")
        process.stdin.flush()
        result=json.loads(process.stdout.readline())
        passed=result.get("status")=="ok" and result["fixed_h_inchi"]==entry["inchi"]
        if not passed:
            results.append(dict(index=i,smiles=entry["mapped_smiles"],expected=entry["inchi"],actual=result))
        if i%1000==0: print(i,"checked;",len(results),"disagreements",flush=True)
    process.stdin.close()
    assert process.wait()==0
    report=dict(schema=1,checkpoint_sha256=MODEL_HASH,model_manifest_sha256=digest(args.model/"model.json"),
                binary_sha256=digest(args.binary),total=len(entries),passed=len(entries)-len(results),disagreements=results)
    args.output.write_text(json.dumps(report,allow_nan=False)+"\n",encoding="utf-8")
    print({k:v for k,v in report.items() if k!="disagreements"})
    raise SystemExit(bool(results))
