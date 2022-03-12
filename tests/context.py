
"""
See:
    https://docs.python-guide.org/writing/structure/
"""
from pathlib import Path
import inspect
import sys

# print(__file__, __name__)
current_file = inspect.getfile(inspect.currentframe())
current_dir = Path(current_file).parent.absolute()
root_dir = current_dir.parent
sys.path.insert(0, str(root_dir))

print(f"current_file: {current_file}")
print(f"current_dir: {current_dir}")
print(f"root_dir: {root_dir}")
print()
[print(p) for p in sys.path]

# import solarenergy
