"""CisLink: CRISPR-trained models that link noncoding variants and cis-regulatory elements to their target genes.
CisLink-Extended is the full model; CisLink-Base needs only accessibility and RNA."""
from .model import MODELS, CisLinkModel, load_model, score

__version__ = "1.5.0"
__all__ = ["MODELS", "CisLinkModel", "load_model", "score", "__version__"]
