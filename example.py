fib_stop = 10
storage = [None for _ in range(fib_stop + 1)]


def fib_recursive_with_storage(n):
    global storage
    if n == 0:  # Base case
        storage[n] = 1
        return 1
    if n == 1:  # Base case
        storage[n] = 1
        return 1
    # Check storage
    if storage[n-1] is not None and storage[n-2] is not None:
        storage[n] = storage[n-1] + storage[n-2]
        return storage[n]
    # Fill storage for n-2
    if storage[n-1] is not None:
        storage[n-2] = fib_recursive_with_storage(n-2)
        storage[n] = storage[n-1] + storage[n-2]
        return storage[n]
    # Fill storage for n-1
    if storage[n-2] is not None:
        storage[n-1] = fib_recursive_with_storage(n-1)
        storage[n] = storage[n-1] + storage[n-2]
        return storage[n]
    storage[n - 1] = fib_recursive_with_storage(n - 1)
    storage[n - 2] = fib_recursive_with_storage(n - 2)
    storage[n] = storage[n - 1] + storage[n - 2]
    return storage[n]
    pass


def fib_dynamic(n):
    storage[0] = 1
    storage[1] = 1
    for i in range(2, n+1):
        storage[i] = storage[i-1] + storage[i-2]
        pass
    return storage[n]


def main():
    global storage
    print(fib_recursive_with_storage(fib_stop))
    print(storage)
    storage = [None for _ in range(fib_stop + 1)]
    print(fib_dynamic(fib_stop))
    print(storage)
    pass


main()