import time

class Numbers:
    def __init__(self, a, b, c) -> None:
        self.a = a
        self.b = b
        self.c = c

    def permuteNumbers(self):
        return Numbers(self.b, self.c, self.a)
    
class Name:
    def __init__(self, name) -> None:
        self.name = name

    def reverseName(self):
        reverse = self.name[::-1]
        return Name(reverse)

class Total:
    def __init__(self, name: Name, numbers: Numbers) -> None:
        self.name = name
        self.numbers = numbers
    
    def changeTotal(self):
        return Total(self.name.reverseName(), self.numbers.permuteNumbers())

class OneTotal:
    def __init__(self, name: str, a: int, b: int, c: int):
        self.name = name
        self.a = a
        self.b = b
        self.c = c

    def changeTotal(self):
        return OneTotal(self.name[::-1], self.b, self.c, self.a)

total = Total(Name('JohnDoe'), Numbers(2,4,0))
oneTotal = OneTotal('JohnDoe', 2, 4, 0)
start_time = time.time()
for _ in range(2000000):
    total = total.changeTotal()
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()
for _ in range(2000000):
    oneTotal = oneTotal.changeTotal()
print("--- %s seconds ---" % (time.time() - start_time))

# Output:
# --- 1.512286901473999 seconds ---
# --- 0.7540619373321533 seconds ---