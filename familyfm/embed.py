"""ESM2 sequence embedding: the 480-dimension sequence block.

Checkpoint facebook/esm2_t12_35M_UR50D, sequences truncated at 3,000 residues,
mean-pooled over residues with BOS and EOS excluded, first 480 dimensions.

Every line is load-bearing and getting any of it wrong is SILENT. In particular
the tokenizer has no lowercase vocabulary: an un-uppercased sequence collapses
to a single <unk> and returns a meaningless vector with no error at all.
"""
import hashlib
import re

import numpy as np

CHECKPOINT = "facebook/esm2_t12_35M_UR50D"
MAX_RESIDUES = 3000
MAX_LEN = 3002
DIMS = 480

_WS = re.compile(r"\s+")
_model = None
_tok = None
_device = None


def normalize(sequence):
    return _WS.sub("", sequence).upper()


def sequence_key(sequence):
    return hashlib.md5(normalize(sequence).encode()).hexdigest()


def _load(device=None):
    global _model, _tok, _device
    if _model is None:
        import torch
        from transformers import AutoModel, AutoTokenizer
        if device is None:
            device = "mps" if torch.backends.mps.is_available() else "cpu"
        _device = device
        _tok = AutoTokenizer.from_pretrained(CHECKPOINT)
        _model = AutoModel.from_pretrained(CHECKPOINT).eval().to(device)
        for p in _model.parameters():
            assert p.dtype == torch.float32, "expected float32 weights"
    return _tok, _model


def embed_residues(sequence, device=None):
    import torch
    tok, model = _load(device)
    seq = normalize(sequence)
    if not seq:
        raise ValueError("empty sequence after normalization")
    enc = tok(seq, return_tensors="pt", truncation=True, max_length=MAX_LEN)
    enc = {k: v.to(_device) for k, v in enc.items()}
    with torch.no_grad():
        h = model(**enc).last_hidden_state[0][1:-1].float().cpu().numpy()
    if h.shape[0] != min(len(seq), MAX_RESIDUES):
        raise RuntimeError(
            f"residue count {h.shape[0]} does not match sequence length "
            f"{len(seq)}; the tokenizer dropped or merged characters.")
    return h


def embed_sequence(sequence, device=None):
    v = embed_residues(sequence, device).mean(axis=0)
    assert v.shape == (DIMS,), v.shape
    return v.astype(np.float32)
