# Single source of truth for scaffold lengths and allowed file extensions.
# Merges the previously divergent dicts from views/__init__.py and forms/__init__.py.

SCAFFOLD_LENGTHS = {
    'm13mp18': 7249,
    'phix174': 5386,
    'p8064': 8064,
    'p7308': 7308,
    'p7560': 7560,
    'custom': 0,
    'debug': 50000,
}

ALLOWED_EXTENSIONS = {'stl'}
