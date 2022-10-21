# scaffold lengths
scafLen = {'m13mp18': 7249, 'phix174': 5386, 'p8064': 8064, 'p7308': 7308, 'p7560': 7560}

# set allowed extensions
ALLOWED_EXTENSIONS = set(['stl'])

from .connections import *
from .draw import *
from .pathway import *
from .stl import *
from .submit_sequences import *
from .verify import *
