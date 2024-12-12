import random

def generate_random_color():
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return f'#{r:02X}{g:02X}{b:02X}'

random_color_table = [generate_random_color() for _ in range(200)]
print(random_color_table)