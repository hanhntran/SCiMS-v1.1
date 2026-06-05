# __init__.py: initialize SCiMS package

from .utils import (
    read_metadata,
    normalize_colname,
    find_sample_id_column,
    extract_sample_id,
)

from .analytical import (
    classify_sample,
    compute_chromosome_probs,
    _empty_result,
)

from .process_input_file import process_input_file

# Package version
from ._version import __version__
 