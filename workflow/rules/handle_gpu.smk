# Source: https://github.com/snakemake/snakemake/issues/281

from contextlib import contextmanager
from itertools import cycle
import GPUtil
from random import randint
from time import sleep


@contextmanager
def get_gpu_id(get_gid=True, n_gpus=2):
    """get next available GPU and lock it

    Note that the lock file created will be at
    /tmp/LCK_gpu_{allocated_gid}.lock

    This is based on the solution proposed at
    https://github.com/snakemake/snakemake/issues/281#issuecomment-610796104
    and then modified slightly

    Parameters
    ----------
    get_gid : bool, optional
        if True, return the ID of the first available GPU. If False,
        return None. This weirdness is to allow us to still use this
        contextmanager when we don't actually want to create a lockfile
    n_gpus : int, optional
        number of GPUs on this device

    Returns
    -------
    allocated_gid : int
        the ID of the GPU to use

    """
    allocated_gid = None
    sleep(randint(0, 5))
    avail_gpus = GPUtil.getAvailable(
        order="random", maxLoad=0.1, maxMemory=0.1, includeNan=False, limit=n_gpus
    )
    if not get_gid:
        avail_gpus = []
    for gid in cycle(avail_gpus):
        # then we've successfully created the lockfile
        if os.system(f"dotlockfile -r 1 /tmp/LCK_gpu_{gid}.lock") == 0:
            allocated_gid = gid
            break
    try:
        yield allocated_gid
    finally:
        os.system(f"dotlockfile -u /tmp/LCK_gpu_{allocated_gid}.lock")
