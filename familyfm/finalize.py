"""Make a trained directory into a servable bundle: target index, predict.py,
reference predictions. Run once per model directory."""
import argparse, hashlib, json, os, shutil, sys
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import embed as emb

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--bundle",required=True)
    ap.add_argument("--measurements",required=True)
    a=ap.parse_args()
    d=pd.read_csv(a.measurements,low_memory=False,
                  usecols=["uniprot","sequence","family"]).drop_duplicates("uniprot")
    idx=json.load(open(os.path.join(a.bundle,"sequence_index.json")))["seq_key_to_row"]
    acc2key={}; fam={}
    for u,s,f in zip(d.uniprot,d.sequence,d.family):
        k=emb.sequence_key(s)
        if k in idx:
            acc2key[u]=k; fam[u]=f
    json.dump({"accession_to_seq_key":acc2key,"family_of":fam},
              open(os.path.join(a.bundle,"target_index.json"),"w"))
    shutil.copy2(os.path.join(HERE,"predict.py"),os.path.join(a.bundle,"predict.py"))
    print("targets resolvable in this bundle: %d"%len(acc2key))

    spec_path=os.path.join(a.bundle,"predict.py")
    import importlib.util
    spec=importlib.util.spec_from_file_location("fp",spec_path)
    fp=importlib.util.module_from_spec(spec); sys.modules["fp"]=fp
    sys.path.insert(0,a.bundle); spec.loader.exec_module(fp)
    m=fp.load(a.bundle)
    SMI=["CC(=O)Oc1ccccc1C(=O)O",
         "CN1CCN(c2ccc(Nc3ncc(F)c(-c4cc(F)c5nc(C)n(C(C)C)c5c4)n3)cc2)CC1",
         "Cc1ccc(cc1Nc1nccc(n1)-c1cccnc1)C(=O)Nc1ccc(C)c(c1)N1CCN(C)CC1",
         "CCOC(=O)c1ccccc1O","c1ccc2c(c1)[nH]c1ccccc12"]
    avail=[x for x in ["P00533","P08684","P31645","Q12809","P11712","P28223","P43220"]
           if x in acc2key]
    ref=[]
    for s in SMI:
        for i in range(len(avail)-1):
            a1,b1=avail[i],avail[i+1]
            ref.append({"smiles":s,"target_a":a1,"target_b":b1,
                        "p":round(fp.compare_targets(m,s,a1,b1),6)})
    selfc=fp.compare_targets(m,SMI[0],avail[0],avail[0])
    anti=fp.compare_targets(m,SMI[0],avail[0],avail[1])+fp.compare_targets(m,SMI[0],avail[1],avail[0])
    sha=hashlib.sha256(open(os.path.join(a.bundle,"familyfm_selectivity.joblib"),"rb").read()).hexdigest()
    json.dump({"tolerance":1e-6,"self_comparison":round(selfc,6),
               "antisymmetry_sum":round(anti,12),"model_sha256":sha,
               "predictions":ref},
              open(os.path.join(a.bundle,"reference_predictions.json"),"w"),indent=1)
    man=json.load(open(os.path.join(a.bundle,"MANIFEST.json")))
    man["model_sha256"]=sha
    man["targets_servable"]=len(acc2key)
    json.dump(man,open(os.path.join(a.bundle,"MANIFEST.json"),"w"),indent=1)
    print("self-comparison %.6f (must be 0.500000)"%selfc)
    print("antisymmetry sum %.12f (must be 1.0)"%anti)
    print("reference predictions: %d | sha256 %s"%(len(ref),sha[:16]))
if __name__=="__main__": main()
