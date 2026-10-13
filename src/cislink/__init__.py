"""CisLink: general models that link noncoding variants (CisLink-V2G) and cis-regulatory elements
(CisLink-E2G) to their target genes."""
from .model import MODELS, CisLinkModel, load_model, score

__version__ = "1.0.0"
__all__ = ["MODELS", "CisLinkModel", "load_model", "score", "__version__"]
