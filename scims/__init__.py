# __init__.py: initialize SCiMS package

from .utils import (
    read_metadata,
    normalize_colname,
    find_sample_id_column,
    read_master_file,
    extract_sample_id
)

from .helpers import (
    load_training_data
)

from .classification import (
    process_sample_xy,
    process_sample_zw
)

from .__main__ import main

__all__ = [
    'main',
    'read_metadata',
    'normalize_colname',
    'find_sample_id_column',
    'read_master_file',
    'extract_sample_id',
    'load_training_data',
    'process_sample_xy',
    'process_sample_zw'
]

# Package version
__version__ = '1.1.0'