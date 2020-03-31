import random


def simulate(n):
    container = [0 for i in range(20)]
    for num in range(n):
        container[int(random.random() * 20)] += 1
    if all(m > 0 for m in container):
        return 1
    else:
        return 0

if __name__ == '__main__':
    n = 160
    count = 0
    for i in range(10000):
        count += simulate(n)
    print(count/10000)
