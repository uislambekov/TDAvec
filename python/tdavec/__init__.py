print("Hello from tdavec")

from .vectorizer import TDAvectorizer
from .vectorizer import createEllipse, pmin, pmax, test_package

# Optional: if you're sure tdavec_core.so is built and available
try:
    from .tdavec_core import compute_diagram
except ImportError:
    compute_diagram = None  # Or raise a custom warning

__all__ = ["TDAvectorizer", "pmin", "pmax", "createEllipse", "test_package"]
