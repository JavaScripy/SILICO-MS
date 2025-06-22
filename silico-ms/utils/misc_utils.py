import numpy as np 
from pathos import multiprocessing as mp
from tqdm import tqdm


def softmin(x: np.ndarray) -> np.ndarray:
    return np.exp(-x) / sum(np.exp(-x))

def softmax(x: np.ndarray) -> np.ndarray: 
    return np.exp(x) / sum(np.exp(x))


def simple_parallel(
    input_list, function, n_threads=12, timeout=4000, max_retries=3
):
    cpus = min(mp.cpu_count(), n_threads)
    pool = mp.ProcessingPool(nodes=cpus)
    results = list(tqdm(pool.imap(function, input_list), total=len(input_list)))
    return results


def chunked_parallel(
    input_list, function, chunck_size = 10000, n_threads=12, timeout=4000, max_retries=3
):
    def batch_func(list_inputs):
        outputs = []
        for i in list_inputs:
            outputs.append(function(i))
        return outputs
    
    chunked_list = [
        input_list[i : i + chunck_size] for i in range(0, len(input_list), chunck_size)
    ]
    
    list_outputs = simple_parallel(
        chunked_list,
        batch_func,
        n_threads=n_threads,
        timeout=timeout,
        max_retries=max_retries,
    )
    # Unroll
    full_output = [j for i in list_outputs for j in i]

    return full_output
