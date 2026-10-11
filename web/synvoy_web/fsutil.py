"""File operations that several parts of the service need to get right."""

import os
import shutil
import stat
import sys
from pathlib import Path


def remove_tree(path, root):
    """Delete ``path``, which must lie inside ``root``. A symbolic link is removed as a
    link; what it points to is never touched."""
    path, root = Path(path), Path(root)
    if path.is_symlink():
        path.unlink()
        return
    if not path.exists():
        return
    resolved, resolved_root = path.resolve(), root.resolve()
    if resolved == resolved_root or resolved_root not in resolved.parents:
        raise RuntimeError("refusing to delete %s: it is not inside %s" % (path, root))

    def retry(function, failed_path, *_):
        # Read-only files and folders (task folders can contain them).
        try:
            os.chmod(os.path.dirname(failed_path), stat.S_IRWXU)
            if not os.path.islink(failed_path):
                os.chmod(failed_path, stat.S_IRWXU)
            function(failed_path)
        except OSError:
            pass

    if sys.version_info >= (3, 12):
        shutil.rmtree(str(path), onexc=retry)
    else:
        shutil.rmtree(str(path), onerror=retry)
    if path.exists():
        raise OSError("could not delete %s completely" % path)


def tree_bytes(path):
    """Bytes of the files under ``path``, or of the file itself. Links are counted as
    links, not followed; what disappears while it is counted is left out."""
    path = Path(path)
    try:
        if path.is_symlink() or not path.is_dir():
            return path.lstat().st_size
    except OSError:
        return 0
    total = 0
    for folder, _, names in os.walk(str(path)):
        for name in names:
            try:
                total += os.lstat(os.path.join(folder, name)).st_size
            except OSError:
                pass
    return total


def write_atomic(path, text):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    os.replace(str(temporary), str(path))
